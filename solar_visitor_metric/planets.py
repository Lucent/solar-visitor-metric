"""Planet score: how close a visitor came to the planets, against visitors of the same speed and perihelion distance in every orientation, arriving at any time.

Closeness is geometry alone. A visitor is its heliocentric two-body hyperbola; the planets move on their osculating orbits at the visitor's epoch. A planet's pass is found where the path crosses the planet's distance from the Sun, refined by straight-line relative motion, which is exact for the close passes that matter. Units: AU, day.

The null keeps what discovery fixes, the speed far from the Sun and the perihelion distance, and randomizes what a planner would choose: the orientation of the path (uniform over all rotations, so a path in the planets' plane is as likely as it is for chance and no likelier) and the time of perihelion, over SPAN years of the planets' motion on one clock. A path aimed into the ecliptic and timed to meet several planets therefore scores for both.

Per planet, p is the chance a null visitor passes that planet at least this close (`Null.p`). Across the planets, the statistic is the product of the three smallest p's, the chance of the visitor's three closest passes taken together: one very close pass dominates it, and two or three moderate ones compound. Its own null distribution, from an independent sample of null visitors, turns it into the score (`Null.bits`), with the planets' dependence on a shared path built in.
"""
import math
from dataclasses import dataclass

import numpy as np
from scipy.stats import beta

MU = 2.959122082855911e-4   # GM_sun, AU^3 / day^2
NAMES = ("Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune")
RADIUS_KM = (2440, 6052, 6371, 3390, 69911, 58232, 25362, 24622)   # mean radii, to read a pass in planet radii
SPAN = 1e4 * 365.25          # days of perihelion times in the null
N_CURVE = 1_000_000          # null visitors for each planet's chance curve
N_STAT = 1_000_000           # independent null visitors for the statistic's distribution; planned tours reach deeper (real_visitors.MISSION_STAT)
FIT = 100                    # closest null passes per planet that fix its small-distance law
CHUNK = 100_000


def rotation(inc, Om, om):
	"""Perifocal to ecliptic: Rz(Om) Rx(inc) Rz(om), broadcast over leading axes."""
	cO, sO, ci, si, cw, sw = np.cos(Om), np.sin(Om), np.cos(inc), np.sin(inc), np.cos(om), np.sin(om)
	return np.stack([np.stack([cO * cw - sO * sw * ci, -cO * sw - sO * cw * ci, sO * si], -1),
					 np.stack([sO * cw + cO * sw * ci, -sO * sw + cO * cw * ci, -cO * si], -1),
					 np.stack([sw * si, cw * si, ci], -1)], -2)


@dataclass
class Orbits:
	"""The planets' osculating heliocentric orbits at epoch t0 (arrays of 8)."""
	a: np.ndarray
	e: np.ndarray
	R: np.ndarray    # (8, 3, 3) perifocal to ecliptic
	M0: np.ndarray
	t0: float

	@classmethod
	def at(cls, sim, t0):
		"""From a REBOUND simulation holding the Sun and the eight planets at JD t0 (real.planets_at, whose own clock starts at zero)."""
		o = [sim.particles[j].orbit(primary=sim.particles[0]) for j in range(1, 9)]
		f = lambda k: np.array([getattr(x, k) for x in o])
		return cls(f("a"), f("e"), rotation(f("inc"), f("Omega"), f("omega")), f("M"), t0)

	def state(self, t):
		"""Positions and velocities at times t of shape (..., 8)."""
		n = np.sqrt(MU / self.a ** 3)
		M = self.M0 + n * (t - self.t0)
		E = M.copy()
		for _ in range(12):
			E -= (E - self.e * np.sin(E) - M) / (1 - self.e * np.cos(E))
		b = self.a * np.sqrt(1 - self.e ** 2)
		x, y = self.a * (np.cos(E) - self.e), b * np.sin(E)
		Edot = n / (1 - self.e * np.cos(E))
		vx, vy = -self.a * np.sin(E) * Edot, b * np.cos(E) * Edot
		return (np.einsum("jab,...jb->...ja", self.R[..., :2], np.stack([x, y], -1)),
				np.einsum("jab,...jb->...ja", self.R[..., :2], np.stack([vx, vy], -1)))


def random_rotations(rng, n):
	"""n rotations uniform over SO(3), from uniform unit quaternions."""
	u = rng.random((3, n))
	q = np.stack([np.sqrt(1 - u[0]) * np.sin(2 * np.pi * u[1]), np.sqrt(1 - u[0]) * np.cos(2 * np.pi * u[1]),
				  np.sqrt(u[0]) * np.sin(2 * np.pi * u[2]), np.sqrt(u[0]) * np.cos(2 * np.pi * u[2])], -1)
	x, y, z, w = q.T
	return np.stack([np.stack([1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)], -1),
					 np.stack([2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)], -1),
					 np.stack([2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)], -1)], -2)


def conic(q, v_inf, R, tp, t):
	"""Positions and velocities at times t (N, 8) on hyperbolas of perihelion q, speed v_inf, orientations R (N, 3, 3), perihelion times tp (N)."""
	e = 1 + q * v_inf ** 2 / MU
	A = MU / v_inf ** 2                 # |a|
	M = v_inf ** 3 / MU * (t - tp[:, None])
	F = np.arcsinh(M / e)
	for _ in range(30):
		F -= (e * np.sinh(F) - F - M) / (e * np.cosh(F) - 1)
	Fdot = v_inf ** 3 / MU / (e * np.cosh(F) - 1)
	s = math.sqrt(e * e - 1)
	xy = np.stack([A * (e - np.cosh(F)), A * s * np.sinh(F)], -1)
	vxy = np.stack([-A * np.sinh(F) * Fdot, A * s * np.cosh(F) * Fdot], -1)
	return np.einsum("nab,njb->nja", R[..., :2], xy), np.einsum("nab,njb->nja", R[..., :2], vxy)


def passes(q, v_inf, R, tp, orbits):
	"""Closest approach (AU) and its time to each planet for hyperbolas of perihelion q (AU), speed v_inf (AU/day), orientations R (N, 3, 3) and perihelion times tp (N): arrays (N, 8). Each pass starts where the path crosses the planet's distance from the Sun, inbound and outbound, or at perihelion if the path stays outside it, and is walked to the nearest separation on the true curves."""
	e = 1 + q * v_inf ** 2 / MU
	p = q * (1 + e)
	N = len(tp)
	best_d = np.full((N, 8), np.inf)
	best_t = np.zeros((N, 8))
	for side in (-1.0, 1.0):
		r = np.broadcast_to(orbits.a, (N, 8)).copy()
		for _ in range(4):
			cos_nu = np.clip((p / r - 1) / e, -1.0, 1.0)
			nu = side * np.arccos(np.where(r > q, cos_nu, 1.0))
			F = 2 * np.arctanh(math.sqrt((e - 1) / (e + 1)) * np.tan(nu / 2))
			t = tp[:, None] + (e * np.sinh(F) - F) * MU / v_inf ** 3
			r = np.linalg.norm(orbits.state(t)[0], axis=-1)
		# walk to the nearest separation: each step is the straight-line correction at the current time
		for _ in range(12):
			X, V = conic(q, v_inf, R, tp, t)
			P, VP = orbits.state(t)
			rel, vrel = X - P, V - VP
			t = t - np.einsum("nja,nja->nj", rel, vrel) / np.einsum("nja,nja->nj", vrel, vrel)
		X, _ = conic(q, v_inf, R, tp, t)
		d = np.linalg.norm(X - orbits.state(t)[0], axis=-1)
		better = d < best_d
		best_d = np.where(better, d, best_d)
		best_t = np.where(better, t, best_t)
	return best_d, best_t


def _null_chunk(args):
	q, v_inf, orbits, m, seed = args
	rng = np.random.default_rng(seed)
	tp = orbits.t0 + (rng.random(m) - 0.5) * SPAN
	return passes(q, v_inf, random_rotations(rng, m), tp, orbits)[0]


def null_passes(q, v_inf, orbits, n, seed, mapper=map):
	"""Closest approaches (n, 8) of n null visitors: random orientation, perihelion at a random time within SPAN of the epoch. `seed` is a tuple of ints, extended per chunk; `mapper` spreads the chunks, e.g. a Pool's map."""
	chunks = [(q, v_inf, orbits, min(CHUNK, n - k), (*seed, k)) for k in range(0, n, CHUNK)]
	return np.vstack(list(mapper(_null_chunk, chunks)))


TOP = 3   # passes in the statistic: of the measures tried, this one detects random planned tours best


def statistic(ps):
	"""The product of the TOP smallest p's, rows of ps; smaller is more surprising."""
	return np.prod(np.sort(np.atleast_2d(ps), axis=-1)[..., :TOP], axis=-1)


class Null:
	"""The planet null for one speed and perihelion distance."""

	def __init__(self, q, v_inf, orbits, mapper=map, seed=0, n_stat=None):
		"""`mapper` may be lazy (a Pool's imap): the statistic's sample is reduced chunk by chunk, so tens of millions of null visitors fit in memory."""
		self.curve = np.sort(null_passes(q, v_inf, orbits, N_CURVE, (seed, 0), mapper), axis=0)
		n = n_stat or N_STAT
		chunks = [(q, v_inf, orbits, min(CHUNK, n - k), (seed, 1, k)) for k in range(0, n, CHUNK)]
		self.T = np.sort(np.concatenate([statistic(self.p(d)) for d in mapper(_null_chunk, chunks)]))

	def p(self, d):
		"""Per planet, the chance a null visitor passes at least as close as d (rows of 8). Below its FIT-th closest null pass a planet's chance follows the area law of a straight pass, d^2."""
		d = np.atleast_2d(d)
		out = np.empty(d.shape)
		for j in range(8):
			c = self.curve[:, j]
			count = np.searchsorted(c, d[:, j], side="right")
			law = FIT / len(c) * (d[:, j] / c[FIT - 1]) ** 2
			out[:, j] = np.where(count >= FIT, count / len(c), law)
		return out

	def bits(self, d):
		"""Score of closest approaches d (8): -log2 of the null chance of a statistic this small, counted among the null visitors. Few null visitors reach a deep statistic, so the chance is the count's 95% upper confidence limit, which never overstates; where none reach it, the score is the most the sample can show, a floor on the truth."""
		t = float(statistic(self.p(d))[0])
		count = int(np.searchsorted(self.T, t, side="right"))
		n = len(self.T)
		return max(0.0, -math.log2(min(1.0, float(beta.ppf(0.95, count + 1, max(n - count, 1))))))


def euler(R):
	"""(inc, Om, om) of a rotation built by `rotation`."""
	return np.arccos(np.clip(R[2, 2], -1, 1)), np.arctan2(R[0, 2], -R[1, 2]), np.arctan2(R[2, 0], R[2, 1])


def aimed(V_in, V_out):
	"""The orientation of the hyperbola whose asymptotic directions of motion are V_in (arriving) and V_out (leaving), both ecliptic: perihelion lies along -(out - in), the plane holds both."""
	u_in, u_out = V_in / np.linalg.norm(V_in), V_out / np.linalg.norm(V_out)
	x = -(u_out - u_in) / np.linalg.norm(u_out - u_in)
	y = (u_out + u_in) / np.linalg.norm(u_out + u_in)
	return np.stack([x, y, np.cross(x, y)], -1)


def plan(q, v_inf, orbits, targets, R=None, seed=0, n=200_000, starts=20):
	"""A planned path: orientation (unless R fixes it) and perihelion time whose passes come as close as possible to `targets`, {planet index: distance AU}. Returns (R, tp, closest approaches)."""
	from scipy.optimize import minimize
	rng = np.random.default_rng(seed)
	idx = np.array(list(targets))
	goal = np.array([targets[j] for j in idx])
	miss = lambda d: np.sum(np.log(np.maximum(d[..., idx], 1e-6) / goal) ** 2, axis=-1)
	tp = orbits.t0 + (rng.random(n) - 0.5) * SPAN
	Rs = random_rotations(rng, n) if R is None else np.broadcast_to(R, (n, 3, 3))
	first = np.argsort(miss(passes(q, v_inf, Rs, tp, orbits)[0]))[:starts]
	scale = 365.25

	def unpack(x):
		Rx = rotation(*x[:3]) if R is None else R
		return Rx[None], np.array([x[-1] * scale])

	best = None
	for k in first:
		x0 = np.array([*euler(Rs[k]), tp[k] / scale]) if R is None else np.array([tp[k] / scale])
		res = minimize(lambda x: float(miss(passes(q, v_inf, *unpack(x), orbits)[0])[0]), x0, method="Nelder-Mead",
					   options={"xatol": 1e-9, "fatol": 1e-12, "maxiter": 4000})
		if best is None or res.fun < best.fun:
			best = res
	Rb, tb = unpack(best.x)
	return Rb[0], float(tb[0]), passes(q, v_inf, Rb, tb, orbits)[0][0]
