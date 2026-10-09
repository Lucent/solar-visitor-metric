"""Score the real interstellar objects in 3D against real planet positions.

Orbits come from JPL SBDB (osculating heliocentric ecliptic elements); planets from JPL Horizons via REBOUND, cached in data/. Units: AU, day.

Visit score: for each planet, the actual closest approach d, and p = the fraction of the planet's orbital phase (its real orbit, phase shifted uniformly) that would have brought it at least this close to the visitor's actual path. Combined across the 8 planets with `order_bits`.

Assist: the visitor is traced back to 60 yr before perihelion, then flown forward in three universes (see `assist`), separating the planets' direct pull from the Sun's reflex wobble, which shifts every visitor's perihelion slightly and is not something a trajectory could "use". Non-gravitational accelerations are ignored.
"""
import json
import math
import os

import numpy as np
import rebound

from .scores import order_bits

HORIZONS_IDS = ["10", "199", "299", "399", "499", "5", "6", "7", "8"]
DATA = os.path.join(os.path.dirname(__file__), "..", "data")


SBDB_URL = "https://ssd-api.jpl.nasa.gov/sbdb.api?sstr={}&full-prec=true"


def elements(tag):
	"""Osculating elements from JPL SBDB, fetched once and cached in data/."""
	path = os.path.join(DATA, f"sbdb_{tag}.json")
	if not os.path.exists(path):
		from urllib.parse import quote
		from urllib.request import urlopen
		with urlopen(SBDB_URL.format(quote(tag)), timeout=60) as r:
			body = r.read()
		if b'"orbit"' not in body:
			raise ValueError(f"SBDB has no orbit for {tag!r}")
		os.makedirs(DATA, exist_ok=True)
		with open(path, "wb") as f:
			f.write(body)
	with open(path) as f:
		d = json.load(f)
	orb = d["orbit"]
	el = {e["name"]: float(e["value"]) for e in orb["elements"] if e["value"] is not None}
	el["epoch"] = float(orb["epoch"])
	o = d["object"]
	el["name"] = o["fullname"].strip()
	# the name without its provisional designation: 'Oumuamua (A/2017 U1) -> 'Oumuamua, C/2019 Q4 (Borisov) -> Borisov
	el["short"] = el["name"].replace(f'{o["prefix"]}/{o["des"]}', "").strip(" ()")
	return el


def planets_at(jd):
	"""Sun + planets at JD from Horizons, cached as a REBOUND binary."""
	path = os.path.join(DATA, f"planets_{jd:.1f}.bin")
	if os.path.exists(path):
		return rebound.Simulation(path)
	sim = rebound.Simulation()
	sim.units = ("day", "AU", "Msun")
	for hid in HORIZONS_IDS:
		sim.add(hid, date=f"JD{jd:.5f}")
	sim.save_to_file(path)
	return sim


def build(tag):
	el = elements(tag)
	sim = planets_at(el["epoch"])
	sim.integrator = "ias15"
	sun = sim.particles[0]
	A = abs(el["a"])
	n = math.sqrt(sim.G / A ** 3)
	M = n * (el["epoch"] - el["tp"])
	sim.add(primary=sun, a=el["a"], e=el["e"], inc=math.radians(el["i"]),
			Omega=math.radians(el["om"]), omega=math.radians(el["w"]), M=M)
	sim.N_active = len(HORIZONS_IDS)
	sim.testparticle_type = 0
	sim.t = el["epoch"]
	return sim, el


def _helio(sim):
	ps = sim.particles
	s = np.array([ps[0].x, ps[0].y, ps[0].z])
	return np.array([[p.x, p.y, p.z] for p in ps]) - s, \
		np.array([[p.vx, p.vy, p.vz] for p in ps]) - np.array([ps[0].vx, ps[0].vy, ps[0].vz])


def track(tag, years=12.0, dt=0.1):
	"""Heliocentric positions of planets and visitor on a uniform grid."""
	sim, el = build(tag)
	t0 = el["tp"]
	times = np.arange(t0 - years * 365.25, t0 + years * 365.25, dt)
	out = np.empty((len(times), len(HORIZONS_IDS) + 1, 3))
	vel = np.empty_like(out)
	# backwards from epoch, then forwards
	back = times[times <= sim.t][::-1]
	fwd = times[times > sim.t]
	state = sim.copy()
	for seq, s in ((back, sim), (fwd, state)):
		for t in seq:
			s.integrate(t, exact_finish_time=1)
			k = int(round((t - times[0]) / dt))
			out[k], vel[k] = _helio(s)
	return times, out, vel, el


def closest(times, pos, vel):
	"""Per planet: (d_min AU, time JD, relative speed AU/day)."""
	iso = pos[:, -1]
	res = []
	for j in range(1, len(HORIZONS_IDS)):
		d = np.linalg.norm(iso - pos[:, j], axis=1)
		k = int(np.argmin(d))
		res.append((float(d[k]), float(times[k]),
					float(np.linalg.norm(vel[k, -1] - vel[k, j]))))
	return res


def kepler_positions(sim_el, M_offsets, times, t_ref, mu):
	"""Positions of a Keplerian orbit for many mean-anomaly offsets."""
	a, e, inc, Om, om, M0 = sim_el
	n = math.sqrt(mu / a ** 3)
	M = M0 + M_offsets[:, None] + n * (times[None, :] - t_ref)
	E = M.copy()
	for _ in range(30):
		E -= (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
	x = a * (np.cos(E) - e)
	y = a * math.sqrt(1 - e * e) * np.sin(E)
	cO, sO, ci, si, cw, sw = math.cos(Om), math.sin(Om), math.cos(inc), math.sin(inc), \
		math.cos(om), math.sin(om)
	X = (cO * cw - sO * sw * ci) * x + (-cO * sw - sO * cw * ci) * y
	Y = (sO * cw + cO * sw * ci) * x + (-sO * sw + cO * cw * ci) * y
	Z = (sw * si) * x + (cw * si) * y
	return np.stack([X, Y, Z], axis=-1)


def phase_p(tag, times, pos, vel, j, d_obs, n_phase=20000):
	"""Fraction of planet j's phase that comes within d_obs of the path."""
	sim, el = build(tag)
	o = sim.particles[j].orbit(primary=sim.particles[0])
	orb = (o.a, o.e, o.inc, o.Omega, o.omega, o.M)
	iso = pos[:, -1]
	r = np.linalg.norm(iso, axis=1)
	q_p, Q_p = o.a * (1 - o.e), o.a * (1 + o.e)
	sel = (r > q_p - d_obs) & (r < Q_p + d_obs)
	if not sel.any():
		return 1.0
	# time step so the visitor moves << d_obs between samples
	v_max = float(np.max(np.linalg.norm(vel[sel, -1], axis=1))) + o.n * o.a
	stride = max(1, int(0.2 * d_obs / (v_max * (times[1] - times[0]))))
	idx = np.flatnonzero(sel)[::stride]
	offs = np.linspace(0, 2 * math.pi, n_phase, endpoint=False)
	hit = np.zeros(n_phase, bool)
	for chunk in np.array_split(idx, max(1, len(idx) // 200)):
		P = kepler_positions(orb, offs, times[chunk], sim.t, sim.G)
		dmin = np.linalg.norm(P - iso[chunk][None, :, :], axis=-1).min(axis=1)
		hit |= dmin <= d_obs
	return max(hit.mean(), 1 / n_phase)


OBLIQUITY = math.radians(23.4392911)
AU_DAY_KMS = 1.495978707e8 / 86400


def ecl_to_icrs(v):
	c, s = math.cos(OBLIQUITY), math.sin(OBLIQUITY)
	return np.array([v[0], c * v[1] - s * v[2], s * v[1] + c * v[2]])


def icrs_to_ecl(v):
	"""Rows of ICRS vectors to the ecliptic frame (inverse of ecl_to_icrs)."""
	c, s = math.cos(OBLIQUITY), math.sin(OBLIQUITY)
	return np.stack([v[:, 0], c * v[:, 1] + s * v[:, 2], -s * v[:, 1] + c * v[:, 2]], 1)


def _two_body(s):
	"""The visitor about one body holding all the mass at the barycenter."""
	com = s.com(last=len(HORIZONS_IDS))
	p = s.particles[-1]
	two = rebound.Simulation()
	two.G = s.G
	two.t = s.t
	two.add(m=com.m, x=com.x, y=com.y, z=com.z, vx=com.vx, vy=com.vy, vz=com.vz)
	two.add(x=p.x, y=p.y, z=p.z, vx=p.vx, vy=p.vy, vz=p.vz)
	return two


def path(tag, times, years=12.0):
	"""Visitor positions (AU, heliocentric ecliptic) at JD `times`.

	N-body among the planets within `years` of perihelion; beyond, the two-body hyperbola about the barycenter, where the planets no longer matter and IAS15 would otherwise step at Mercury's pace for a million years.
	"""
	sim, el = build(tag)
	tp = el["tp"]
	sim.integrate(tp, exact_finish_time=1)
	out = np.empty((len(times), 3))
	for sign in (-1, 1):
		s = sim.copy()
		two = None
		side = np.flatnonzero(sign * (times - tp) >= 0)
		for k in side[np.argsort(np.abs(times[side] - tp))]:
			if abs(times[k] - tp) <= years * 365.25:
				s.integrate(times[k], exact_finish_time=1)
				out[k] = _helio(s)[0][-1]
				continue
			if two is None:
				s.integrate(tp + sign * years * 365.25, exact_finish_time=1)
				two = _two_body(s)
			two.integrate(times[k], exact_finish_time=1)
			c, p = two.particles
			out[k] = (p.x - c.x, p.y - c.y, p.z - c.z)
	return out


def _inbound(s):
	"""Inbound asymptotic velocity (motion at t -> -inf), barycentric."""
	p = s.particles[-1]
	com = s.com(last=len(HORIZONS_IDS))
	r = np.array([p.x - com.x, p.y - com.y, p.z - com.z])
	v = np.array([p.vx - com.vx, p.vy - com.vy, p.vz - com.vz])
	mu = s.G * com.m
	h = np.cross(r, v)
	evec = np.cross(v, h) / mu - r / np.linalg.norm(r)
	e = np.linalg.norm(evec)
	f_inf = math.acos(-1 / e)
	ph = evec / e
	qh = np.cross(h / np.linalg.norm(h), ph)
	vinf = math.sqrt(max(np.dot(v, v) - 2 * mu / np.linalg.norm(r), 0.0))
	return -vinf * (math.cos(f_inf) * ph - math.sin(f_inf) * qh)


def asymptotes(tag, years=60.0):
	"""(V_in, V_out) in km/s, ICRS, from states 60 yr before/after perihelion."""
	sim, el = build(tag)
	s = sim.copy()
	s.integrate(el["tp"] - years * 365.25)
	V_in = _inbound(s)
	s = sim.copy()
	s.integrate(el["tp"] + years * 365.25)
	V_out = _outbound(s)
	return ecl_to_icrs(V_in) * AU_DAY_KMS, ecl_to_icrs(V_out) * AU_DAY_KMS


def _outbound(s):
	"""Outbound asymptotic velocity of the visitor (last particle), barycentric."""
	p = s.particles[-1]
	com = s.com(last=len(HORIZONS_IDS))
	r = np.array([p.x - com.x, p.y - com.y, p.z - com.z])
	v = np.array([p.vx - com.vx, p.vy - com.vy, p.vz - com.vz])
	mu = s.G * com.m
	h = np.cross(r, v)
	evec = np.cross(v, h) / mu - r / np.linalg.norm(r)
	e = np.linalg.norm(evec)
	f_inf = math.acos(-1 / e)
	ph = evec / e
	qh = np.cross(h / np.linalg.norm(h), ph)
	vinf = math.sqrt(max(np.dot(v, v) - 2 * mu / np.linalg.norm(r), 0.0))
	return vinf * (math.cos(f_inf) * ph + math.sin(f_inf) * qh)


def _cancel_direct(simp):
	"""Additional force removing the planets' pull on the visitor."""
	s = simp.contents
	ps = s.particles
	v = ps[s.N - 1]
	G = s.G
	for k in range(1, len(HORIZONS_IDS)):
		pk = ps[k]
		dx, dy, dz = pk.x - v.x, pk.y - v.y, pk.z - v.z
		r3 = (dx * dx + dy * dy + dz * dz) ** 1.5
		f = G * pk.m / r3
		v.ax -= f * dx
		v.ay -= f * dy
		v.az -= f * dz


def assist(tag, years_back=60.0):
	"""Planets' change to the outbound asymptote, split into parts.

	Three universes from the same state 60 yr before perihelion: real; real but with the planets' direct pull on the visitor canceled (the Sun still wobbles); and Sun-only (all mass at the barycenter). Returns (total, direct, reflex) as fractions of v_inf, and v_inf (AU/day).
	"""
	sim, el = build(tag)
	sim.integrate(el["tp"] - years_back * 365.25)
	t_end = el["tp"] + years_back * 365.25
	outs = {}
	for mode in ("real", "no_direct", "sun_only"):
		s = sim.copy()
		if mode == "no_direct":
			s.additional_forces = _cancel_direct
			s.force_is_velocity_dependent = 0
		if mode == "sun_only":
			# one body of the total mass at the barycenter (else the Sun drifts off its wobble and shifts the perihelion)
			com = s.com(last=len(HORIZONS_IDS))
			sun = s.particles[0]
			sun.m = com.m
			sun.x, sun.y, sun.z = com.x, com.y, com.z
			sun.vx, sun.vy, sun.vz = com.vx, com.vy, com.vz
			for k in range(1, len(HORIZONS_IDS)):
				s.particles[k].m = 0.0
		s.integrate(t_end)
		outs[mode] = _outbound(s)
	v0 = float(np.linalg.norm(outs["sun_only"]))
	frac = lambda a, b: float(np.linalg.norm(outs[a] - outs[b]) / v0)
	return {"total": frac("real", "sun_only"), "direct": frac("real", "no_direct"),
			"reflex": frac("no_direct", "sun_only")}, v0


def score(tag):
	times, pos, vel, el = track(tag)
	close = closest(times, pos, vel)
	ps = [phase_p(tag, times, pos, vel, j, d) for j, (d, _, _) in enumerate(close, start=1)]
	a, vinf = assist(tag)
	return {"name": el["name"], "close": close, "p": ps, "visit_bits": order_bits(ps),
			"assist": a, "v_inf_kms": vinf * 1.495978707e8 / 86400, "times": times,
			"pos": pos}
