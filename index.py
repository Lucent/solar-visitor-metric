# sweep_visits_ephemeris.py
# Minimal sweep using real ephemerides (poliastro + Astropy), no manual positions,
# no error handling. Measure-correct sampling on the orbital plane.
#
# What it does:
# - Samples inbound states on log-spaced rings (R_min..R_max) with uniform angles (left→right).
# - Propagates each test particle under Sun-only gravity (two-body) from epoch0 forward.
# - Counts "visits" when closest approach to a planet ≤ 5% of that planet's Hill radius.
# - Conditions the probability on actually "grazing the system": only trials that cross the 30 AU gate
#   are included in the statistics (both numerator and denominator).
# - Reports hit histogram and a simple "1 in x" rarity for ≥ Y hits.
#
# Install:
#   pip install poliastro astropy numpy

import numpy as np
from astropy import units as u, constants as const
from astropy.time import Time, TimeDelta
from poliastro.bodies import Sun, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune
from poliastro.twobody import Orbit, propagation
from functools import lru_cache

# ------------------------
# Config (edit as needed)
# ------------------------
epoch0 = Time("2025-11-05 00:00:00", scale="tdb")  # start epoch
R_MIN  = 30.0 * u.AU            # inner "encircling" ring ~ just outside Neptune
R_MAX  = 3000.0 * u.AU          # "effectively infinity"
NR     = 4                      # number of rings (log-spaced)
NTHETA = 36                     # headings (left→right) per ring (uniform)
NPERI  = 3                      # number of target pericenters (log-spaced)
PERI_MIN = 0.2 * u.AU           # innermost pericenter to sample
PERI_MAX = 2.5 * u.AU           # outer pericenter (still inside the belt)
V_INF  = 30.0 * (u.km/u.s)      # inbound hyperbolic excess speed (at infinity)
T_MAX  = 8.0 * u.year           # integrate each particle this long
DT     = 20.0 * u.day           # step
ALPHA  = 0.05                   # visit threshold = 5% Hill radius
R_OUTER_GATE = 30.0 * u.AU      # start counting once inside this radius
R_BELT_GATE  = 3.0 * u.AU       # require trajectories to penetrate the asteroid belt

# Scoring (explodes with hits + closeness)
GAMMA_PER_PLANET = 10.0
K_POWER = 2.0

# Targets for rarity estimate
Y_MIN_HITS = 2                   # e.g., rarity of "≥ 2 planets visited"

# ------------------------
# Bodies & constants
# ------------------------
BODIES = [Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune]
BODY_BY_NAME = {b.name: b for b in BODIES}
MU_SUN = Sun.k  # GM of the Sun

# ------------------------
# Helpers
# ------------------------
def hill_radius(body, epoch):
    """Instantaneous Hill radius using osculating a and mass ratio at given epoch."""
    orb = Orbit.from_body_ephem(body, epoch)
    a = orb.a.to(u.AU)
    m_body = (body.k / const.G).to(u.kg)
    m_sun  = (Sun.k  / const.G).to(u.kg)
    RH = a * (m_body/(3*m_sun))**(1/3)
    return RH.to(u.AU)

@lru_cache(maxsize=None)
def _planet_r_cached(body_name, jd_tdb):
    body = BODY_BY_NAME[body_name]
    return Orbit.from_body_ephem(body, Time(jd_tdb, format="jd", scale="tdb")).r.to(u.AU)


def planet_r(body, epoch):
    """Heliocentric position vector at epoch (AU)."""
    return _planet_r_cached(body.name, epoch.tdb.jd)

def inbound_state(R, theta, target_peri, epoch):
    """Construct inbound state targeting a specific pericenter inside the system."""
    R = R.to(u.AU)
    rp = target_peri.to(u.AU)

    mu = MU_SUN.to(u.AU**3 / u.day**2)
    v_inf = V_INF.to(u.AU / u.day)

    # Impact parameter producing the requested pericenter for a hyperbolic flyby
    b = np.sqrt(rp**2 + (2 * mu * rp) / v_inf**2)
    h = b * v_inf  # specific angular momentum
    r_hat = np.array([np.cos(theta), np.sin(theta), 0.0])
    t_hat = np.array([-np.sin(theta), np.cos(theta), 0.0])

    # Position vector on the sampling ring
    r_vec = u.Quantity(r_hat, unit=u.one) * R

    # Speed decomposition at the sampling radius
    v_total = np.sqrt(v_inf**2 + 2 * mu / R)
    v_tan = (h / R).to(u.AU / u.day)
    term = (v_total**2 - v_tan**2).to(u.AU**2 / u.day**2)
    if term.value < 0:
        term = 0 * term.unit  # guard numeric noise; will surface if configuration invalid
    v_rad = -np.sqrt(term)

    v_vec = (u.Quantity(r_hat, unit=u.one) * v_rad +
             u.Quantity(t_hat, unit=u.one) * v_tan).to(u.AU / u.day)

    return r_vec, v_vec

def advance_two_body_sun(r0, v0, epoch, steps):
    """Yield (epoch_i, r_i) along two-body (Sun) propagation with fixed DT steps."""
    orb = Orbit.from_vectors(Sun, r0, v0, epoch)
    te = epoch
    yield te, orb.r.to(u.AU)
    for _ in range(steps):
        te = te + TimeDelta(DT)
        orb = propagation.propagate(orb, DT, method=propagation.cowell)
        yield te, orb.r.to(u.AU)

# ------------------------
# Precompute Hill radii at epoch0
# ------------------------
HILLS = {b.name: hill_radius(b, epoch0) for b in BODIES}

# ------------------------
# Sampling grid (measure-correct):
# - Uniform theta per ring
# - Log-spaced rings
#   => approximates uniform sampling in impact parameter at the 30 AU gate
# ------------------------
rings   = u.AU * np.geomspace(R_MIN.value, R_MAX.value, NR)
thetas  = np.linspace(-np.pi, np.pi, NTHETA, endpoint=False)
periaps = u.AU * np.geomspace(PERI_MIN.value, PERI_MAX.value, NPERI)


def run_sweep():
    counted_hits  = []
    counted_score = []

    belt_crossers = 0
    total_trials = len(rings) * len(periaps) * len(thetas)

    for R in rings:
        for rp in periaps:
            for th in thetas:
                r0, v0 = inbound_state(R, th, rp, epoch0)

                # Track min distances after crossing the outer gate
                dmin = {b.name: np.inf * u.AU for b in BODIES}
                crossed_outer = False
                crossed_belt = False

                steps = int((T_MAX / DT).decompose())  # number of integration steps
                for te, r_vec in advance_two_body_sun(r0, v0, epoch0, steps):
                    r_now = np.linalg.norm(r_vec.value) * u.AU

                    # Start counting only after entering the outer system
                    if (not crossed_outer) and (r_now <= R_OUTER_GATE):
                        crossed_outer = True

                    if (not crossed_belt) and (r_now <= R_BELT_GATE):
                        crossed_belt = True

                    if crossed_outer:
                        # Update closest approach to each planet at this time
                        for b in BODIES:
                            rb = planet_r(b, te)
                            d = np.linalg.norm((r_vec - rb).value) * u.AU
                            if d < dmin[b.name]:
                                dmin[b.name] = d

                    # Optional early exit: once we have gone back outside the outer gate
                    if crossed_outer and (r_now > (R_OUTER_GATE + 2 * u.AU)):
                        break

                if not crossed_belt:
                    # Missed the asteroid belt entirely: EXCLUDE from denominator
                    continue

                belt_crossers += 1

                # Visits + simple score for belt-crossers only
                hits = 0
                score = 1.0
                for b in BODIES:
                    y = ALPHA * HILLS[b.name]
                    if dmin[b.name] <= y:
                        hits += 1
                        score *= (y / dmin[b.name]).decompose().value ** K_POWER
                        score *= GAMMA_PER_PLANET  # ×10 per visited planet (default)

                counted_hits.append(hits)
                counted_score.append(score)

    counted_hits_arr = np.array(counted_hits)
    counted_score_arr = np.array(counted_score)

    # ------------------------
    # Rarity ("1 in x") conditioned on crossing the gate
    # ------------------------
    S_obs_min = GAMMA_PER_PLANET ** Y_MIN_HITS  # minimal score if Y hits happen right at the threshold
    mask_Y = counted_score_arr >= S_obs_min     # "≥ Y hits" in this simple scoring
    p = np.mean(mask_Y) if counted_score_arr.size > 0 else np.nan
    rarity = (np.inf if (p == 0) else (1.0 / p)) if np.isfinite(p) else np.nan

    print(f"Belt-crossing trials used: {counted_score_arr.size} of total {belt_crossers} (from {total_trials} samples)")
    print(f"P(≥ {Y_MIN_HITS} hits | crosses {R_BELT_GATE}) ≈ {p:.3e}" if np.isfinite(p) else "No belt-crossers sampled.")
    print(f"→ Rarity ≈ 1 in {rarity:.2g}" if np.isfinite(rarity) else "")

    # Hit histogram
    if counted_hits_arr.size > 0:
        unique, counts = np.unique(counted_hits_arr, return_counts=True)
        print("Hit histogram (hits: count):", dict(zip(unique.tolist(), counts.tolist())))

    return counted_hits_arr, counted_score_arr


if __name__ == "__main__":
    run_sweep()
