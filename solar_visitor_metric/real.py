"""Score the real interstellar objects in 3D against real planet positions.

Orbits come from JPL SBDB (osculating heliocentric ecliptic elements); planets from JPL Horizons via REBOUND, cached in data/. Units: AU, day.

Planet passes are scored by `planets` on the visitor's two-body hyperbola from these elements; the N-body path here is drawn by the site, and its asymptotes feed the route. Non-gravitational accelerations are ignored.
"""
import json
import math
import os

import numpy as np
import rebound

from . import planets

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


COMETS_URL = "https://ssd-api.jpl.nasa.gov/sbdb_query.api?fields=full_name,e,i,q,first_obs&sb-kind=c&sb-cdata=%7B%22AND%22%3A%5B%22e%7CGT%7C0.9%22%5D%7D"
ERA = 5   # years either side of a visitor's perihelion whose comet discoveries measure its surveys


def survey_excess(year):
	"""How much more often than random the surveys of a visitor's era found paths at each tilt from the ecliptic: per planets.TILT_EDGES band, the share of long-period comets first seen within ERA years of `year`, over the share random orientations give. Long-period comets arrive from nearly every direction, so what the surveys found shows what they favor. Sungrazers (q < 0.1 AU, found by coronagraphs) and interstellar objects (e > 1.5) are left out. Returns (excess, comets counted)."""
	path = os.path.join(DATA, "comets.json")
	if not os.path.exists(path):
		from urllib.request import urlopen
		with urlopen(COMETS_URL, timeout=120) as r:
			body = r.read()
		if b'"data"' not in body:
			raise ValueError("SBDB returned no comets")
		with open(path, "wb") as f:
			f.write(body)
	with open(path) as f:
		rows = json.load(f)["data"]
	inc = np.radians([float(i) for name, e, i, q, first in rows
					  if name.strip().startswith("C/") and first and float(q) > 0.1 and float(e) < 1.5 and abs(int(first[:4]) - year) <= ERA])
	share = np.histogram(np.arccos(np.abs(np.cos(inc))), planets.TILT_EDGES)[0] / len(inc)
	return share / np.diff(1 - np.cos(planets.TILT_EDGES)), len(inc)


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


def score(tag, mapper=map):
	"""The visitor's passes of the eight planets and their score against the planet null at its speed and perihelion distance."""
	el = elements(tag)
	orbits = planets.Orbits.at(planets_at(el["epoch"]), el["epoch"])
	v = math.sqrt(planets.MU / abs(el["a"]))
	R = planets.rotation(*(np.array([math.radians(el[k])]) for k in ("i", "om", "w")))
	d, t = planets.passes(el["q"], v, R, np.array([el["tp"]]), orbits)
	null = planets.Null(el["q"], v, orbits, mapper)
	excess, n_comets = survey_excess(round(2000 + (el["tp"] - 2451545.0) / 365.25))
	return {"name": el["name"], "short": el["short"], "v_inf_kms": v * AU_DAY_KMS, "q": el["q"], "tp": el["tp"],
			"d": d[0], "jd": t[0], "p": null.p(d[0])[0], "bits": null.bits(d[0]),
			"segments": null.segments(el["q"], v, orbits, R[0], d[0], excess, mapper), "excess": excess, "comets": n_comets,
			"tilt": math.degrees(float(planets.tilt(R[0])))}
