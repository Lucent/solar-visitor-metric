"""Route-score calibration under arrivals that move as stars do, against isotropic arrivals, on the real catalog.

	.venv/bin/python kinematic_check.py > out/kinematic.txt      # about 45 minutes on 18 cores

The route p's assume isotropic arrival directions. Interstellar objects, natural or not, inherit stellar motions, so this draws each null visitor's direction from the stars moving at its speed (`route.kinematic_directions`) and checks P(bits >= x) * 2^x stays at or below 1.
"""
import sys
from multiprocessing import Pool

import numpy as np

from solar_visitor_metric import route
from solar_visitor_metric.constants import VISITORS

N_NULL = 72_000


def main():
	field, errors, _ = route.Field.from_catalogs("data/gaia_reach.csv", "data/hipparcos_reach.csv")
	with Pool() as pool:
		for name, v in VISITORS.items():
			d = route.kinematic_directions(field, v)
			for null, dirs in (("isotropic", None), ("kinematic", d)):
				nr = route.null_route_parallel(field, v, N_NULL, errors, pool, seed=900, directions=dirs)
				b = -np.log2(np.maximum(nr[:, 2], 1e-300))
				print(f"{name} {v:g} km/s, {null:<9} ({len(d)} stars in the speed band): P(bits >= x) * 2^x "
					  + "  ".join(f"{x}: {(b >= x).mean() * 2 ** x:.2f}" for x in (1, 2, 4, 6, 8)), flush=True)


if __name__ == "__main__":
	sys.exit(main())
