"""Score interstellar objects: planet passes, and the route between stars (Gaia DR3 and Hipparcos).

	.venv/bin/python real_visitors.py > out/real.txt         # 1I, 2I, 3I, and the site's data
	.venv/bin/python real_visitors.py 4I                      # any SBDB designation

Needs network on first run (JPL SBDB/Horizons, Gaia archive); results are cached in data/.
"""
import datetime
import itertools
import json
import math
import sys
from multiprocessing import Pool

import numpy as np

from solar_visitor_metric import planets, real, route, simbad
from solar_visitor_metric.constants import AU_KM

OUT = "out"
SITE = "site/data.json"
NAMES = "data/simbad_names.csv"
TAGS = ("1I", "2I", "3I")
N_NULL = 36_000
# The site's time axis: t - t_peri = TAU sinh(u asinh(SPAN / TAU)) for u in [-1, 1], days through the perihelion passage and millions of years out to the stars on one slider.
TAU = 10.0
SPAN = 3e6 * 365.25
N_PATH = 801
GAIA_EPOCH = 2457389.0  # J2016.0, JD
PC_AU = 206264.806
KMS_AU_DAY = 86400 / 1.495978707e8
SHOWN = 5               # best stars per end on the site
# Published searches for a visitor's birthplace, shown beside its origin.
SEARCHES = {"3I": [("Pérez-Couto et al. 2025", "https://arxiv.org/abs/2509.03361"),
				   ("Guo et al. 2025", "https://arxiv.org/abs/2509.07678")]}
# The site's simulated missions, at 3I's speed among the planets as they were for 3I: a tour of three planets, a route between two stars everyone knows, and that route's path timed to pass two planets.
MISSION_TAG = "3I"
MISSION_SPEED = 58.0
MISSION_Q = 1.0                               # AU, the tour's perihelion
MISSION_TOUR = {2: 0.05, 3: 0.05, 4: 0.05}    # Earth, Mars and Jupiter, each within 0.05 AU
MISSION_ROUTE = ("HIP 32349", "HIP 71683")   # Sirius, then α Cen A
MISSION_TIMED = {2: 0.01, 3: 0.01}            # the route's path, timed for Earth and Mars
MISSION_STAT = 50_000_000                     # null visitors for a planned tour's statistic, which lies far in the tail


def date(jd):
	return (datetime.datetime(2000, 1, 1, 12)
			+ datetime.timedelta(days=jd - 2451545.0)).strftime("%Y-%m-%d")


def main():
	field, errors, n_untraced = route.Field.from_catalogs("data/gaia_reach.csv", "data/hipparcos_reach.csv")
	gaia = np.char.startswith(field.ids.astype(str), "Gaia")
	print(f"Stars with radial velocities that could meet a visitor faster than {route.V_FLOOR:g} km/s within "
		  f"{route.HOP:g} pc of the Sun: {len(field.D)}, {gaia.sum()} from Gaia DR3 and {(~gaia).sum()} from Hipparcos "
		  f"({(field.D < route.HOP).sum()} within {route.HOP:g} pc today; {n_untraced} more there lack "
		  f"a usable radial velocity and are excluded); median velocity error "
		  f"{np.median(errors.star_v[gaia]):.2f} km/s per axis for Gaia, {np.median(errors.star_v[~gaia]):.2f} for Hipparcos\n")
	site = []
	pool = Pool()
	for tag in (sys.argv[1:] or TAGS):
		r = real.score(tag, pool.map)
		V_in, V_out = real.asymptotes(tag)
		print(f"{r['name']}  v_inf {r['v_inf_kms']:.2f} km/s, perihelion {r['q']:.3f} AU")
		print(f"  {'planet':<9}{'closest':>10}{'radii':>10}  {'date':<11}{'p':>8}")
		for name, radius, d, t, p in zip(planets.NAMES, planets.RADIUS_KM, r["d"], r["jd"], r["p"]):
			print(f"  {name:<9}{d:>8.3f}AU{d * AU_KM / radius:>10,.0f}  {date(t):<11}{p:>8.3f}")
		print(f"  planet score: {r['bits']:.2f} bits, against the same speed and perihelion in every orientation at any time")
		s = r["segments"]
		print(f"  where they come from ({s['total']:.2f} bits on their own sample): survey bias {s['bias']:.2f} (from {r['comets']} comets, "
			  f"excess by tilt {np.round(r['excess'], 2).tolist()}), lying in the plane {s['plane']:.2f}, direction within it {s['direction']:.2f}, timing {s['timing']:.2f}")
		ends = []
		shown = {}
		turn = math.degrees(math.acos(V_in @ V_out / (np.linalg.norm(V_in) * np.linalg.norm(V_out))))
		print(f"  the Sun turned it through {turn:.1f}°; the departure is scored given the arrival, over orientations of that turn")
		for V, sign, lab, (p, k, info) in ((V_in, -1, "origin", route.end_score(field, V_in, -1, errors)),
										   (V_out, +1, "destination", route.turn_score(field, V_in, V_out, errors))):
			ends.append(p)
			shown[lab] = encounters(field, V, info)
			cand = np.flatnonzero(info["weights"])
			j = cand[np.argmin(info["d"][cand])]
			hop = np.linalg.norm(V) * abs(info["t"][j]) * route.KMS_PC_MYR
			print(f"  route {lab}: {max(0.0, -math.log2(p)):.2f} bits ({len(cand)} candidate stars); "
				  f"nearest pass: {field.ids[j]}, {field.D[j]:.1f} pc away now, missed by "
				  f"{info['d'][j]:.2f} ± {info['sigma'][j]:.2f} pc {hop:.1f} pc from the Sun, "
				  f"{abs(info['t'][j]):.2f} Myr {'ago' if sign < 0 else 'ahead'}")
		rb = max(0.0, -math.log2(route.combine(*ends)))
		plan = planned(field, float(np.linalg.norm(V_in)), errors)
		print(f"  a route planned at this speed from the nearest reachable star ({field.ids[plan['origin']]}) "
			  f"to the next ({field.ids[plan['destination']]}), perihelion {plan['q']:.2f} AU: {plan['bits']:.1f} bits")
		nr = route.null_route_parallel(field, float(np.linalg.norm(V_in)), N_NULL, errors, pool, seed=7)
		frac = (nr[:, 2] <= route.combine(*ends)).mean()
		print(f"  route combined: {rb:.2f} bits; {frac:.1%} of {N_NULL:,} random visitors "
			  f"at this speed score at least as high on this catalog\n")
		el = real.elements(tag)
		u = np.linspace(-1, 1, N_PATH)
		times = el["tp"] + TAU * np.sinh(u * math.asinh(SPAN / TAU))
		site.append({
			"tag": tag, "name": r["name"], "short": el["short"], "vinf": r["v_inf_kms"], "tp": el["tp"], "q": el["q"],
			"path": real.path(tag, times),
			"planets": {"bits": r["bits"], "segments": r["segments"], "comets": r["comets"], "tilt": r["tilt"], "passes": [{"name": name, "d": d, "radii": d * AU_KM / radius, "jd": t, "p": p}
				for name, radius, d, t, p in zip(planets.NAMES, planets.RADIUS_KM, r["d"], r["jd"], r["p"])]},
			"together": together(r["bits"], rb), "searches": SEARCHES.get(tag, []),
			"route": {"bits": rb, "beaten": frac, "turn": turn, "planned": plan, **shown}})
	if not sys.argv[1:]:
		write_site(field, site, missions(field, errors, pool.imap))
	pool.close()


def together(planet_bits, route_bits):
	"""Bits for planets and route as one test: their p's combined exactly, so looking in two places is paid for (calibrated jointly by together_check.py)."""
	return max(0.0, -math.log2(route.combine(2 ** -planet_bits, 2 ** -route_bits)))


def missions(field, errors, mapper):
	"""Three simulated missions at MISSION_SPEED among the planets as they were for MISSION_TAG: a planned tour of three planets, the star-to-star route, and the route's own path timed to pass two planets. A part a mission does not aim at scores nothing."""
	el = real.elements(MISSION_TAG)
	orbits = planets.Orbits.at(real.planets_at(el["epoch"]), el["epoch"])
	v = MISSION_SPEED * KMS_AU_DAY
	_, _, tour = planets.plan(MISSION_Q, v, orbits, MISSION_TOUR)
	tour_bits = planets.Null(MISSION_Q, v, orbits, mapper, n_stat=MISSION_STAT).bits(tour)
	i, j = (int(np.flatnonzero(field.ids == d)[0]) for d in MISSION_ROUTE)
	V_in, V_out, q = route.plan_route(field, i, j, MISSION_SPEED)
	route_bits = float(route.planted_bits(V_in, V_out, field, errors)[2])
	_, _, timed = planets.plan(q, v, orbits, MISSION_TIMED, R=planets.aimed(*real.icrs_to_ecl(np.array([V_in, V_out]))))
	timed_bits = planets.Null(q, v, orbits, mapper, n_stat=MISSION_STAT).bits(timed)
	passes = lambda d, targets: {planets.NAMES[k]: float(d[k]) for k in targets}
	rows = [{"kind": "tour", "passes": passes(tour, MISSION_TOUR), "planets": tour_bits, "route": 0.0},
			{"kind": "route", "passes": {}, "planets": 0.0, "route": route_bits},
			{"kind": "both", "passes": passes(timed, MISSION_TIMED), "planets": timed_bits, "route": route_bits}]
	for row in rows:
		row["together"] = together(row["planets"], row["route"])
		print(f"Mission ({row['kind']}): passes {row['passes']}; planets {row['planets']:.1f} bits, route {row['route']:.1f}, together {row['together']:.1f}")
	return {"speed": MISSION_SPEED, "q": MISSION_Q, "route_q": q, "origin": i, "destination": j, "rows": rows}


def planned(field, v, errors):
	"""The reference itinerary at this speed on the same catalog: from the nearest star it can meet to the next one a Sun flyby can turn it to, scored as the catalog would see it."""
	(i, _), = route.reachable(field, (1,), v, -1)
	for rank in itertools.count(2):
		(j, _), = route.reachable(field, (rank,), v, +1, {i})
		plan = route.plan_route(field, i, j, v)
		if plan:
			V_in, V_out, q = plan
			return {"bits": float(route.planted_bits(V_in, V_out, field, errors)[2]),
					"origin": int(i), "destination": int(j), "q": q}


def encounters(field, V, info):
	"""One end's candidate count and its SHOWN best stars by look-elsewhere-corrected bits."""
	w = info["weights"]
	cand = np.flatnonzero(w)
	ratio = info["p"][cand] / w[cand]
	best = cand[np.argsort(ratio)[:SHOWN]]
	v = np.linalg.norm(V)
	def star(j):
		hop = v * abs(info["t"][j]) * route.KMS_PC_MYR
		# the pass as an aim: the miss and the catalog's blur on it, as angles seen from the Sun at the meeting
		return {"i": int(j), "bits": max(0.0, -math.log2(min(1.0, info["p"][j] / w[j]))),
				"miss": info["d"][j], "sigma": info["sigma"][j], "hop": hop,
				"aim": math.degrees(math.atan2(info["d"][j], hop)), "blur": math.degrees(math.atan2(info["sigma"][j], hop)),
				"jd": GAIA_EPOCH + info["t"][j] * 1e6 * 365.25}
	return {"candidates": len(cand), "stars": [star(j) for j in best]}


def write_site(field, site, missions):
	"""The site's data: every star any end could have met, ecliptic AU and AU/day from the Gaia epoch, with each visitor's paths and scores indexing into them."""
	keep = sum(route.end_weights(field, v["vinf"], sign) for v in site for sign in (-1, 1)) > 0
	keep[[missions["origin"], missions["destination"]]] = True
	keep = np.flatnonzero(keep)
	index = {j: k for k, j in enumerate(keep)}
	for v in site:
		for end in ("origin", "destination"):
			for s in v["route"][end]["stars"]:
				s["i"] = index[s["i"]]
			v["route"]["planned"][end] = index[v["route"]["planned"][end]]
	for end in ("origin", "destination"):
		missions[end] = index[missions[end]]
	data = {"time": {"tau": TAU, "span": SPAN}, "epoch": GAIA_EPOCH,
			"stars": {"ids": list(field.ids[keep]),
					  "names": [n if n != i else "" for i, n in simbad.names(field.ids[keep], NAMES).items()],
					  "pos": real.icrs_to_ecl(field.pos[keep] * PC_AU),
					  "vel": real.icrs_to_ecl(field.vel[keep] * KMS_AU_DAY)},
			"visitors": site, "missions": missions}
	with open(SITE, "w") as f:
		json.dump(data, f, default=compact, separators=(",", ":"))


def compact(a):
	"""Arrays as flat lists at five significant digits: the site draws, it does not recompute."""
	return [float(f"{v:.5g}") for v in a.ravel()]


if __name__ == "__main__":
	sys.exit(main())
