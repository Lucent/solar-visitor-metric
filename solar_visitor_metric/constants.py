"""Units and the known visitors.

Velocities in the route's Sun-only turn are in units of Earth's mean orbital speed (~29.78 km/s), so GM_sun = 1 with lengths in AU.
"""
import math

AU_KM = 1.495978707e8
V_UNIT_KMS = 2 * math.pi * AU_KM / (365.25 * 86400)  # ~29.785 km/s

# v_inf of the known interstellar objects, km/s
VISITORS = {"1I/Oumuamua": 26.3, "2I/Borisov": 32.3, "3I/ATLAS": 58.0}
