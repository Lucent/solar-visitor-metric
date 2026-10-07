"""Planar patched-conic flyby physics on circular, coplanar planet orbits.

A heliocentric conic is described timing-free by its specific energy E and signed angular momentum h (h > 0 is prograde, same sense as the planets). Because planet phases are treated as independent and uniformly random, the only thing that matters at each ring crossing is the local velocity there, so no Kepler timing is needed.

At a crossing of ring a with local velocity v = (v_r, v_t), the planet's along-track offset s maps linearly onto the planetocentric impact parameter b_p = s * v_r / |w|, with w = v - v_planet. A uniformly random phase is therefore a uniformly random b_p with density |w| / (2 pi a |v_r|).
"""
import math
from dataclasses import dataclass

from .constants import PLANETS


def ecc(E, h):
	return math.sqrt(max(0.0, 1 + 2 * E * h * h))


def perihelion(E, h):
	return h * h / (1 + ecc(E, h))


def aphelion(E, h):
	e = ecc(E, h)
	return math.inf if e >= 1 else h * h / (1 - e)


def local_velocity(E, h, r, sign):
	"""(v_r, v_t) on the conic at radius r; sign of v_r given (-1 inbound)."""
	vr2 = 2 * E + 2 / r - (h / r) ** 2
	return sign * math.sqrt(max(vr2, 0.0)), h / r


def crossings(E, h, r0, sign, skip=()):
	"""Ring crossings in time order, starting at radius r0 moving `sign`.

	Returns (list of (planet_index, v_r sign), bound). A bound orbit is followed for one radial period and flagged so callers can treat capture separately. Planets in `skip` are omitted (already encountered: a return to the same planet is phase-correlated, not independent).
	"""
	q, Q = perihelion(E, h), aphelion(E, h)
	idx = [i for i, p in enumerate(PLANETS) if i not in skip and q < p.a < Q]
	inner = sorted((i for i in idx if PLANETS[i].a < r0), key=lambda i: -PLANETS[i].a)
	outer = sorted((i for i in idx if PLANETS[i].a > r0), key=lambda i: PLANETS[i].a)
	bound = Q < math.inf
	if sign < 0:
		seq = [(i, -1) for i in inner] + [(i, +1) for i in inner[::-1]] + [(i, +1) for i in outer]
		if bound:
			seq += [(i, -1) for i in outer[::-1]]
	else:
		seq = [(i, +1) for i in outer]
		if bound:
			seq += [(i, -1) for i in outer[::-1]] + [(i, -1) for i in inner] + [(i, +1) for i in inner[::-1]]
	return seq, bound


@dataclass
class Crossing:
	planet: int
	v_r: float
	v_t: float

	@property
	def p(self):
		return PLANETS[self.planet]

	@property
	def w(self):
		"""Planetocentric velocity vector (radial, tangential)."""
		return self.v_r, self.v_t - self.p.v_circ

	@property
	def w_mag(self):
		return math.hypot(*self.w)

	@property
	def v_mag(self):
		return math.hypot(self.v_r, self.v_t)

	@property
	def density(self):
		"""Probability per unit b_p (both signs) for a random planet phase."""
		return self.w_mag / (2 * math.pi * self.p.a * abs(self.v_r))

	def b_surface(self):
		"""Smallest |b_p| that misses the planet's surface."""
		R, w2 = self.p.radius, self.w_mag ** 2
		return R * math.sqrt(1 + 2 * self.p.mu / (R * w2))

	def b_for_rp(self, rp):
		"""|b_p| giving closest approach rp."""
		return rp * math.sqrt(1 + 2 * self.p.mu / (rp * self.w_mag ** 2))

	def b_for_dv(self, frac):
		"""Largest |b_p| whose flyby changes heliocentric v by >= frac*|v|.

		Returns 0 if the planet cannot deliver that much even when grazing.
		"""
		s = frac * self.v_mag / (2 * self.w_mag)
		if s >= 1:
			return 0.0
		half = math.asin(s)
		return self.p.mu / (self.w_mag ** 2 * math.tan(half))

	def deflect(self, b_p):
		"""Apply a flyby with signed impact parameter b_p.

		Returns (v_r', v_t', |dv|). Positive b_p means positive planetocentric angular momentum, which turns w counterclockwise.
		"""
		wr, wt = self.w
		wm = self.w_mag
		delta = 2 * math.atan(self.p.mu / (abs(b_p) * wm * wm))
		if b_p < 0:
			delta = -delta
		c, s = math.cos(delta), math.sin(delta)
		# rotate in the (r, t) plane, r x t = +z
		wr2, wt2 = c * wr - s * wt, s * wr + c * wt
		dv = 2 * wm * abs(math.sin(delta / 2))
		return wr2, wt2 + self.p.v_circ, dv


def state_after(cr, b_p):
	"""(E, h, sign, |dv|) of the heliocentric conic after a flyby."""
	v_r, v_t, dv = cr.deflect(b_p)
	a = cr.p.a
	E = 0.5 * (v_r * v_r + v_t * v_t) - 1 / a
	return E, a * v_t, (1 if v_r >= 0 else -1), dv
