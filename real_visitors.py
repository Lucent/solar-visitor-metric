"""Score interstellar objects: planet visits (3D), assist, route (Gaia DR3).

	.venv/bin/python real_visitors.py > out/real.txt         # 1I, 2I, 3I
	.venv/bin/python real_visitors.py 4I                      # any SBDB designation

Needs network on first run (JPL SBDB/Horizons, Gaia archive); results are cached in data/.
"""
import datetime
import math
import sys
from multiprocessing import Pool

import numpy as np

from solar_visitor_metric import real, route, svg
from solar_visitor_metric.constants import PLANETS

OUT = "out"
TAGS = ("1I", "2I", "3I")
N_NULL = 36_000


def date(jd):
	return (datetime.datetime(2000, 1, 1, 12)
			+ datetime.timedelta(days=jd - 2451545.0)).strftime("%Y-%m-%d")


def main():
	field, errors, n_no_rv = route.Field.from_gaia("data/gaia_reach.csv")
	print(f"Gaia DR3: {len(field.D)} stars with radial velocities that could meet a visitor faster "
		  f"than {route.V_FLOOR:g} km/s within {route.HOP:g} pc of the Sun "
		  f"({(field.D < route.HOP).sum()} within {route.HOP:g} pc today; {n_no_rv} more there lack "
		  f"a radial velocity and are excluded); median velocity error "
		  f"{np.median(errors.star_v):.2f} km/s per axis\n")
	figs = []
	for tag in (sys.argv[1:] or TAGS):
		r = real.score(tag)
		V_in, V_out = real.asymptotes(tag)
		print(f"{r['name']}  v_inf {r['v_inf_kms']:.2f} km/s")
		print(f"  {'planet':<9}{'closest':>10}{'R_Hill':>9}  {'date':<11}{'p':>8}")
		for p, (d, t, _), pv in zip(PLANETS, r["close"], r["p"]):
			print(f"  {p.name:<9}{d:>8.3f}AU{d / p.hill:>9.1f}  {date(t):<11}{pv:>8.3f}")
		a = r["assist"]
		print(f"  visit score: {r['visit_bits']:.2f} bits (8 planets, look-elsewhere corrected)")
		print(f"  assist: total {a['total']:.2e} of v_inf; parts (vectors, may cancel): "
			  f"planets' direct pull {a['direct']:.2e}, Sun's reflex wobble {a['reflex']:.2e}")
		ends = []
		turn = math.degrees(math.acos(V_in @ V_out / (np.linalg.norm(V_in) * np.linalg.norm(V_out))))
		print(f"  the Sun turned it through {turn:.1f}°; the departure is scored given the arrival, over orientations of that turn")
		for V, sign, lab, (p, k, info) in ((V_in, -1, "from", route.end_score(field, V_in, -1, errors)),
										   (V_out, +1, "to", route.turn_score(field, V_in, V_out, errors))):
			ends.append(p)
			cand = np.flatnonzero(info["weights"])
			j = cand[np.argmin(info["d"][cand])]
			hop = np.linalg.norm(V) * abs(info["t"][j]) * route.KMS_PC_MYR
			print(f"  route {lab}: {max(0.0, -math.log2(p)):.2f} bits ({len(cand)} candidate stars); "
				  f"nearest pass: Gaia DR3 {field.ids[j]}, {field.D[j]:.1f} pc away now, missed by "
				  f"{info['d'][j]:.2f} ± {info['sigma'][j]:.2f} pc {hop:.1f} pc from the Sun, "
				  f"{abs(info['t'][j]):.2f} Myr {'ago' if sign < 0 else 'ahead'}")
		rb = max(0.0, -math.log2(route.combine(*ends)))
		with Pool() as pool:
			nr = route.null_route_parallel(field, float(np.linalg.norm(V_in)), N_NULL,
										   errors, pool, seed=7)
		frac = (nr[:, 2] <= route.combine(*ends)).mean()
		print(f"  route combined: {rb:.2f} bits; {frac:.1%} of {N_NULL:,} random visitors "
			  f"at this speed score at least as high on the Gaia field\n")
		r["lines"] = [f"v∞ {r['v_inf_kms']:.1f} km/s; visit score {r['visit_bits']:.1f} bits",
					  f"assist {a['total']:.1e} of v∞ (direct {a['direct']:.1e})",
					  f"route {rb:.1f} bits (Gaia DR3, hops ≤ {route.HOP:g} pc)"]
		figs.append(r)
	if not sys.argv[1:]:
		svg.real_figure(f"{OUT}/real_visitors.svg", figs)


if __name__ == "__main__":
	sys.exit(main())
