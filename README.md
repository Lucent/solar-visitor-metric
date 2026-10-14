# solar-visitor-metric

**If something built a probe and sent it through our solar system, would its path give it away?**

Three interstellar objects have crossed the solar system in under a decade: 1I/ʻOumuamua in 2017, 2I/Borisov in 2019, and 3I/ATLAS in 2025. Avi Loeb has argued that ʻOumuamua might be an artificial light sail, and that 3I/ATLAS's path is suspiciously well arranged: a retrograde orbit lying almost exactly in the ecliptic, with close passes by Venus, Mars and Jupiter. Arguments like these rest on the word *suspicious*, and that word needs a number: how often would chance produce the same thing?

solar-visitor-metric provides that number. It borrows the idea behind assembly theory (an object is evidence of design when chance cannot plausibly build it) and turns every claim about a trajectory into **bits of surprise**: −log₂ of the probability under an explicit, calibrated null model. Ten bits means one in a thousand. Twenty means one in a million.

**[Explore the visitors →](https://hyperbolic.lucent.tools/)** Each path in 3D among the planets and the stars it could have come from, from a million years ago to a million years ahead, with every score behind it.

## Findings

**One visitor passed the planets more closely than chance usually allows. None followed a route between stars.**

| | v∞ | Perihelion | Planets | Route between stars | Together |
|---|---|---|---|---|---|
| 1I/ʻOumuamua | 26.4 km/s | 0.26 AU | 2.8 bits | 0.5 bits | 1.6 bits |
| 2I/Borisov | 32.3 km/s | 2.01 AU | 0.05 bits | 0.2 bits | 0 bits |
| 3I/ATLAS | 58.0 km/s | 1.36 AU | **8.2 bits** | 0 bits | 5.5 bits |

**3I/ATLAS threaded three planets: 8.2 bits, about one chance in 300.** It passed Mars at 0.194 AU on 3 October 2025, Venus at 0.649 AU on 3 November 2025 and Jupiter at 0.358 AU on 16 March 2026. Among visitors arriving at its speed and perihelion distance in any orientation and at any time, 2.1% pass Mars that close, 1.0% Venus and 0.3% Jupiter, and the three together happen about once in 300 arrivals. What makes it so is the plane of its orbit, within 5° of the ecliptic, where the planets are: holding its orientation as observed and randomizing only its timing leaves 2.4 bits. Drawing arrivals from the directions the catalog's stars actually move in, rather than evenly from the whole sky, makes the passes rarer still: 8.7 bits. A planned tour at its speed sets the scale: passing Earth, Mars and Jupiter each within 0.05 AU scores 16.8 bits, one chance in 115,000.

**ʻOumuamua's best is Earth, Borisov's nothing.** 1I passed Earth at 0.161 AU on 14 October 2017, five days before it was discovered; 1.7% of random arrivals at its speed and perihelion come that close, and with its other passes it scores 2.8 bits. 2I passed no planet more closely than a random arrival typically does.

**None came from or is heading for a star.** Traced back and forward against the 79,001 stars of Gaia DR3 and Hipparcos that could have met a visitor within 20 parsecs of the Sun, the best route is ʻOumuamua's, 0.5 bits; half of all random arrivals at its speed do as well. 3I's closest origin, the young flare-star pair AT Mic B, which it missed by 0.39 parsecs (1.8°) 210,000 years ago, is a one-in-3,800 pass on its own; but the best of about 2,000 stars tried typically reaches one in 2,000 by chance, so it scores nothing. Two independent searches of Gaia for 3I's birthplace found dozens of sub-parsec encounters and no plausible parent ([Pérez-Couto et al. 2025](https://arxiv.org/abs/2509.03361); [Guo et al. 2025](https://arxiv.org/abs/2509.07678)).

**A planned path would stand out.** Three missions flown at 3I's speed among the planets as they were for 3I, each scored by the same tests:

| Simulated mission at 58 km/s | Planets | Route | Together |
|---|---|---|---|
| A tour passing Earth, Mars and Jupiter, each within 0.05 AU | 16.8 bits | 0 | 13.1 bits |
| From Sirius around the Sun to α Centauri A | 0 | 13.4 bits | 10.1 bits |
| That route's own path, timed to pass Earth at 0.26 AU and Mars at 0.54 AU | 2.9 bits | 13.4 bits | 12.7 bits |

Across 60 random planned tours, one to four planets each passed within 0.03–0.5 AU, the planet score averages 13.4 bits. On the real catalog, with its own errors, a route planted from the nearest star a visitor can reach to the next scores 7.5 bits at ʻOumuamua's speed, 8.2 at Borisov's and 11.8 at 3I's.

**Three closest passes, taken together.** Each planet's pass is a chance in its own right, how often a random arrival comes that close to that planet; the score is the chance of the three smallest, taken together and calibrated against random arrivals. Of ten ways tried to combine the eight planets, all calibrated alike, it catches random planned tours best: 13.4 bits on average, against 12.0 for Fisher's method over all eight and 8.6 for the nearest pass measured in planet radii. A single very close pass dominates it, and two or three moderate ones compound. The statistic was chosen on planned tours before the visitors were scored with it.

**A stellar itinerary is detectable, and precision is the lever.** On a synthetic neighborhood with realistic stellar motions, a visitor planted to come from one nearby star and turn at the Sun for another scores 11.6–15.1 bits; one threading far stars scores almost nothing, because the stars' own velocity errors smear where they were a million years ago. Every threefold improvement in stellar velocities buys about 4 bits: 3.3 bits at 3 km/s errors, 11.6 at 0.3 km/s, 20.2 at 0.03 km/s.

**A traveler favors near stops, but not one star.** Which stars a visitor would plausibly aim for is a weight, separate from how precisely it aimed. Weights falling as the fourth power of the meeting distance, each doubling of distance getting half the share of the one before, with meetings nearer than 2 parsecs counted as 2, detect planted routes best of 24 falloffs and floors tried: 8.1 bits on average, against 3.7 with every star weighted equally. Without the floor, the one or two stars met nearest the Sun take most of the budget. No choice among the 24 gives any visitor a route: the best, ʻOumuamua's, is 0.9 bits before paying for having tried 24.

**The famous neighbors change nothing.** Gaia saturates on the brightest stars and measures radial velocities only for stars bright and cool enough, so Sirius, Vega, Procyon and α Centauri A are missing from it. Hipparcos with ground-based radial velocities adds 949 such stars, about one candidate in twenty at each end of a route, without any of them scoring.

**The stars that matter are mostly not near us now.** A hop is measured where the visitor and the star meet, not where the star is today, and stars outrun slow visitors: of the 2,089 stars ʻOumuamua could have come from within 20 parsecs of the Sun, four in five are farther than that now, some nearly 200 parsecs away. A search confined to today's neighborhood misses most of the candidates.

**A route has to be scored as one trip.** Scoring where a visitor came from and where it is heading as two independent tests overstates chance routes: a pair of stars on one line through the Sun is hit at both ends by any visitor aimed down that corridor, and counted twice. The score treats a route as origin, Sun, destination: the destination is judged given the arrival and the angle the Sun turned it through, so only the orientation of the turn is free. A path the Sun barely bends earns nothing for a destination straight ahead; a sharp turn that lands on a star is rare, and scores.

**A mismeasured star nearly gave ʻOumuamua a home.** On the first pass through Gaia, 1I scored 4.3 bits for its origin. The star responsible, Gaia DR3 1928116685826146304, lies 42 parsecs away and appeared to have passed 1.3 parsecs from the Sun 63,000 years ago, just where ʻOumuamua then was. Its cataloged radial velocity is 654 ± 7 km/s, against 21 km/s across the sky. All twelve of its Gaia spectra were deblended from a brighter companion 1.1″ away at the same distance, the configuration Gaia's own validation flags for spurious high velocities. Taken at face value the pair would be crossing the Galaxy at 880 km/s, well beyond escape speed and faster than any binary could be thrown without being torn apart. That such a star topped the list is no coincidence: a spurious radial velocity, added to a genuine motion across the sky, points a star's track almost through the Sun, which manufactures exactly the short, cheap hop the score looks at first. The statistics held (the star wins for 0.29% of random visitors at 1I's speed, and 4.3 bits was a one-in-twenty result, reported as such), but a route built on it would have been an artifact. Radial velocities now count only at Gaia's recommended signal-to-noise (≥ 5, which sets aside 6,886) and for stars the Galaxy can hold (under 600 km/s about its center, which removes six, this one among them).

## How it works

**Planets.** A visitor's path is its two-body hyperbola about the Sun, and the planets move on their real orbits; planetary gravity plays no part, because closeness is geometry. For every planet, the chance that a random arrival passes it at least as closely is measured against a million random arrivals with the visitor's speed and perihelion distance, which discovery fixes, in uniformly random orientations and at random times across 10,000 years of the planets' motion, which a planner would choose. A path in the planets' plane and timed to meet several of them therefore scores for both. The three smallest chances are multiplied, and that product is calibrated against another million random arrivals, so planets passed together on one path are allowed for. Earth counts like any other planet: visitors are found because they pass near us, which favors a close Earth pass, and a visitor interested in life would favor one too.

**Route.** Inbound and outbound asymptotes are straight lines through a field of moving stars. For the origin, the exact chance that an isotropic arrival at this speed passes each star this closely is computed by integrating over the sphere of possible relative velocities; it is the solid angle of the miss, so it measures the aim in degrees, the same for any star. For the destination, the chance is taken over the circle of departures the Sun's turn allows, given the arrival. Both are convolved with the measurement noise of the miss. Which star was worth aiming for is a separate weight, falling as the fourth power of the distance at which the visitor could meet it, nearer than 2 parsecs counted as 2: the preference for near stops of the [traveling salesman problem](https://en.wikipedia.org/wiki/Traveling_salesman_problem), the score's main assumption about what a planned route looks like. The two ends do not say the same thing. A high origin score means a recent stellar origin, which a natural object ejected by a nearby star within the last few million years would also show; only the destination, a next star the path is aimed at, has no natural reason. The chances assume arrivals from every direction alike, while real interstellar objects, natural or not, move as stars do; drawn instead from the motions of the catalog's own stars at each visitor's speed, random arrivals reach a route score no more often than its bits allow, P(bits ≥ x) · 2ˣ between 0.67 and 0.96 from 1 to 8 bits, against 0.72 to 0.89 for directions at random, so the score reports one number.

**Together.** Planets and route combine as two chances, p₁p₂(1 − ln p₁p₂), so looking in both places is paid for and a strong result in one still counts. Both scores randomize the arrival direction, so they could move together under chance; scored for both on the same 36,000 random paths at each visitor's speed and perihelion, they are uncorrelated (the correlation of their logarithms stays within ±0.05) and the combination stays calibrated, P(bits ≥ x) · 2ˣ between 0.68 and 0.97 from 1 to 8 bits.

**What it does not test.** The planned paths are close passes and stellar itineraries, what a probe surveying the planets or stopping at stars would leave behind. A large telescope could survey Earth from several AU, and a probe could go somewhere without passing anything, so 0 bits rules out a visit of this kind, not design in general. Propulsion is out of scope too: orbits are taken as JPL fitted them, and ʻOumuamua's non-gravitational acceleration, the anomaly the debate turns on, is not scored.

**The null is the product.** Every score is −log₂ of a chance under an explicit model of random arrivals, and every one is checked against chance arrivals to confirm that P(score ≥ x) never exceeds 2⁻ˣ: a score is allowed to understate evidence, never to overstate it.

Orbits come from JPL's Small-Body Database, planet states from JPL Horizons, and stars from Gaia DR3, with Hipparcos (van Leeuwen 2007) and Pulkovo radial velocities (Gontcharov 2006) for those Gaia cannot measure, named from SIMBAD, all fetched on first use. When a fourth interstellar object is discovered, `real_visitors.py 4I` scores it.

## Run

```
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt && mkdir -p out
.venv/bin/python -m pytest -q
.venv/bin/python real_visitors.py > out/real.txt          # 1I, 2I, 3I, the missions, and the site's data
.venv/bin/python statistic_check.py > out/statistic.txt   # how to combine the planets, by power on planned tours
.venv/bin/python route_check.py > out/route.txt           # the route score on a synthetic star field
.venv/bin/python falloff_check.py > out/falloff.txt       # how steeply the route favors near stars, by power on planted routes
.venv/bin/python kinematic_check.py > out/kinematic.txt   # route calibration under arrivals moving as stars do
.venv/bin/python together_check.py > out/together.txt     # planets and route calibrated together
.venv/bin/python -m http.server -d site                   # the explorer, locally
```

Figures are written as SVG to `out/`; the explorer in `site/` is published to GitHub Pages on every push that changes it. Architecture, decisions and conventions are in [AGENTS.md](AGENTS.md).
