"""Units and planet data.

Internal units: length = AU, time = yr/(2*pi), so GM_sun = 1 and the velocity unit is Earth's mean orbital speed (~29.78 km/s).
"""
from dataclasses import dataclass
import math

AU_KM = 1.495978707e8
V_UNIT_KMS = 2 * math.pi * AU_KM / (365.25 * 86400)  # ~29.785 km/s


def kms(v):
	"""km/s -> internal velocity units."""
	return v / V_UNIT_KMS


@dataclass(frozen=True)
class Planet:
	name: str
	a: float       # semimajor axis, AU (orbits treated as circular)
	mu: float      # mass / M_sun (= GM in internal units)
	radius: float  # physical radius, AU

	@property
	def v_circ(self):
		return self.a ** -0.5

	@property
	def hill(self):
		return self.a * (self.mu / 3) ** (1 / 3)

	@property
	def soi(self):
		"""Laplace sphere of influence."""
		return self.a * self.mu ** 0.4


PLANETS = (
	Planet("Mercury", 0.38710, 1.6601e-7, 2440 / AU_KM),
	Planet("Venus",   0.72333, 2.4478e-6, 6052 / AU_KM),
	Planet("Earth",   1.00000, 3.0035e-6, 6371 / AU_KM),
	Planet("Mars",    1.52366, 3.2272e-7, 3390 / AU_KM),
	Planet("Jupiter", 5.20289, 9.5479e-4, 69911 / AU_KM),
	Planet("Saturn",  9.53668, 2.8589e-4, 58232 / AU_KM),
	Planet("Uranus", 19.18916, 4.3662e-5, 25362 / AU_KM),
	Planet("Neptune", 30.06992, 5.1514e-5, 24622 / AU_KM),
)

# v_inf of the known interstellar objects, km/s
VISITORS = {"1I/Oumuamua": 26.3, "2I/Borisov": 32.3, "3I/ATLAS": 58.0}
