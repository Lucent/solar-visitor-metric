# solar-visitor-metric — agent & architecture guide

The craft doctrine is the one that governs Ptable and Ptree (`../Ptable/STYLE.md`, `../Ptable/AGENTS.md`, `../Ptree/AGENTS.md`): single-developer craft, design by elimination, fix causes not cases, trust internal invariants and let impossible states crash loudly. This project is Python and science rather than a browser app, so the DOM/CSS rules do not apply; everything about honesty, deletion, and root causes does.

## Mission

Quantify how surprising a visitor's path is, as bits against an explicit null, so a planned path (gravity assists, close passes, a stellar itinerary) would stand out and a natural one would not. A score is only as good as its null: **every score must be a calibrated p-value.** Under the null, P(score ≥ x) must not exceed 2⁻ˣ; conservative is acceptable, overstating is a bug. Any change to a score or its null reruns the calibration checks (tests and the scripts' calibration figures) before the change is trusted.

## Map

- `constants.py` — units (AU, yr/2π, GM☉ = 1) and planets on circular orbits.
- `flyby.py` — planar patched conics: a conic is (E, h); at each ring crossing a random planet phase is a uniform planetocentric impact parameter b_p, and a flyby rotates the planet-relative velocity.
- `chain.py` — semi-analytic chance of 1–3 flybys that each change heliocentric velocity by ≥ 10%, recursing over where in each window the flyby happened.
- `scores.py` — one path walker (`_walk`) with two nulls: independent random phase per crossing (`fly`) and real J2000 phases on a common clock swept over heading and epoch (`fly_timed`). Visit score, assist, and their Monte Carlo nulls.
- `route.py` — 3D far-field route score: straight asymptotes against linearly moving stars, synthetic or Gaia DR3, each field complete for every star that could meet a visitor within the hop limit.
- `nbody.py` — REBOUND checks of the patched conic: single flybys, targeted chains with the removal test, planned multi-flyby paths.
- `real.py` — the real visitors in 3D: SBDB orbits, Horizons planets, per-planet phase p-values, assist split into direct pull and solar reflex.
- `svg.py` — plain-SVG figures, light colors inline with a dark-mode stylesheet.

## Decisions

- **Planar first, 3D as a factor.** The planar null is the working model; the extra 3D chance per flyby (crossing within the target distance of the planet's plane) is measured separately by `factor3d.py` and applied later. The real visitors are scored fully in 3D.
- **Use is not proximity.** A qualifying flyby changes heliocentric velocity by ≥ 10%. The sphere of influence and Hill fractions only measure closeness; at interstellar speeds they say little about use.
- **Random phases are the null, checked against a common clock.** Independent random planet phases make the null timing-free and fast; `fly_timed` with real phases verifies it. Where they differ (a planet's inbound and outbound crossings share one phase), the score is built to be right under both.
- **Visit score: one p per planet, order statistic across planets.** p is the chance a random phase passes at least this close (by impact parameter), folded over a planet's crossings as 1 − (1 − p_min)^n. The score is n · min_k P(at least k of n planets this close), so one very close pass and several fairly close ones both count, without Fisher's dilution of a single strong pass.
- **Assist: change to the outbound asymptote against a Sun-only universe.** The counterfactual puts all mass at the barycenter; letting the Sun drift off its wobble fakes a perihelion shift. For real visitors the planets' direct pull and the Sun's reflex wobble are reported separately.
- **Visit and assist add only over independent planets.** Combined in-system bits are assist bits plus the visit score of planets that did not bend the path.
- **Route: exact per-star p, noise convolved, short hops cheap.** Per star, the exact chance an isotropic visitor at this speed passes this close, convolved with the Rice distribution of the measured miss (otherwise noise below a barely reachable star's floor scores as impossible). Look-elsewhere is a weighted Bonferroni, 1/rank by hop among stars reachable at this speed; reachability and hop depend only on speed, so restricting to them keeps the test valid.
- **A route is one itinerary: origin, Sun, destination.** The arrival end is scored against isotropic arrivals; the departure end given the arrival, where the Sun's turn angle is taken as observed (it is set by the perihelion, and discovery favors close perihelia) and only its orientation about the arrival direction is random. The two p's are then independent by construction and combine exactly. Scoring the ends separately counted a pair of stars on one line through the Sun twice (a corridor; it overstated 3I-speed nulls by 1.3× in the tail); conditioned, a path the Sun barely bends earns nothing for a destination straight ahead.
- **A hop is where the meeting happens, not where the star is today.** Candidates are the stars a visitor at v∞ could meet within 20 pc of the Sun, ranked by that distance. Stars outrun slow visitors, so most candidates are farther than 20 pc now; a field is complete for visitors faster than 20 km/s when it holds every star within 20 pc × (1 + |w| / 20 km/s) today. Convolution is spent only on stars whose lower bound on p is below their weight, since no other star can set the score.
- **A bad radial velocity fakes an encounter.** Real proper motion plus a spurious radial speed of hundreds of km/s puts a star's track almost through the Sun. Gaia velocities count only with `rv_expected_sig_to_noise` ≥ 5 (Katz et al. 2023) and a Galactic rest-frame speed under 600 km/s, above the Galaxy's escape speed.
- **N-body is the referee, not the engine.** The patched conic carries the statistics; REBOUND confirms deflections, targeted chains, and asymptotes, and the visitor starts relative to the barycenter because far out it orbits that, not the wobbling Sun.

## Conventions

- **Tabs**, per `.editorconfig`. Comments explain why. **No hard line breaks in prose** — docstrings, comments, and docs are one line per paragraph or bullet.
- **No speculative error handling.** Validate only external data at its boundary (SBDB, Horizons, Gaia responses); trust internal invariants.
- **Generated output is regenerated, never hand-edited.** `out/` belongs to the scripts; `data/` is a fetch cache. Neither is committed.
- **AGENTS.md records decisions and operating instructions, not results; the README presents the findings.** When a number in the README changes, rerun the script that produced it and update the README in the same commit; when a doc and the code disagree, fix the doc.
- Python packages only in the project `.venv`. Parallel scripts keep their work under `if __name__ == "__main__"` (Python 3.14 starts workers with forkserver).
- Commits are atomic: one logical change with its tests and doc lines.
