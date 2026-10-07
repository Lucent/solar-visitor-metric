# solar-visitor-metric

**If something built a probe and sent it through our solar system, would its path give it away?**

Three interstellar objects have crossed the solar system in under a decade: 1I/ʻOumuamua in 2017, 2I/Borisov in 2019, and 3I/ATLAS in 2025. Avi Loeb has argued that ʻOumuamua might be an artificial light sail, and that 3I/ATLAS's path is suspiciously well arranged: a retrograde orbit lying almost exactly in the ecliptic, with close passes by Venus, Mars and Jupiter. Arguments like these rest on the word *suspicious*, and that word needs a number: how often would chance produce the same thing?

solar-visitor-metric provides that number. It borrows the idea behind assembly theory (an object is evidence of design when chance cannot plausibly build it) and turns every claim about a trajectory into **bits of surprise**: −log₂ of the probability under an explicit, calibrated null model. Ten bits means one in a thousand. Twenty means one in a million.

**[Explore the visitors →](https://hyperbolic.lucent.tools/)** Each path in 3D among the planets and the stars it could have come from, from a million years ago to a million years ahead, with every score behind it.

## Findings

**None of the three known visitors shows any trace of a planned path.**

| | v∞ | Planet visits | Planets' pull on its outbound path | Star-to-star route |
|---|---|---|---|---|
| 1I/ʻOumuamua | 26.4 km/s | 0.2 bits | 0.13% of v∞ | 0 bits |
| 2I/Borisov | 32.3 km/s | 0 bits | 0.02% | 0 bits |
| 3I/ATLAS | 58.0 km/s | 0 bits | 0.12% | 0 bits |

3I/ATLAS passed Jupiter at 0.358 AU on 16 March 2026, almost exactly one Hill radius, and Mars at 0.194 AU on 3 October 2025. That looks pointed until you ask how often a randomly placed Jupiter comes that close: **3% of the time**, and Mars 19%. Allow for the fact that eight planets were available to be "visited" and the surprise drops to zero. The Jupiter pass changed 3I's velocity by about a thousandth, a hundred times short of anything a mission designer would call a gravity assist. Traced back and forward against the 78,052 Gaia DR3 stars that could have met a visitor within 20 parsecs of the Sun, none of the three came from or is heading for any of them more closely than chance routinely manages: four random visitors in five at 3I's speed score at least as well.

The one genuinely unusual fact about 3I is the one Loeb leads with. Its orbit lies within 5° of the ecliptic, which an isotropic arrival does about once in 275 tries: **8 bits**, a real coincidence, but a single one among the many features that could have been checked.

**A planned path would stand out, and by a wide margin.** Chance gives two velocity-changing flybys (≥ 10% each) to roughly one in three million in-plane visitors at ʻOumuamua's speed (**22 bits**), and one in a hundred million at 3I's (**27 bits**). A chance chain falls off roughly as the fourth power of speed, because fast objects barely bend. Real paths flown in a full N-body simulation score accordingly:

| Planned path at 26 km/s | Visits | Assist | Combined |
|---|---|---|---|
| Earth and Mars at 1.5 million km, then a Jupiter assist | 13.0 bits | 10.3 bits | **17 bits** |
| Earth and Mars at 150,000 km, then a Jupiter assist | 21.1 bits | 10.3 bits | **24 bits** |
| Jupiter assist alone | 4.0 bits | 10.3 bits | 10 bits |

A targeted Jupiter → Saturn chain flown in REBOUND misses Saturn by 5.2 AU when Jupiter's mass is switched off: every significant chain is causally dependent on its first flyby, so the second encounter is no coincidence. And all of these are in-plane numbers. In three dimensions each flyby must also hit the planet's orbital plane to within its target distance, a further factor of about 2 × 10⁻³ per flyby by direct Monte Carlo, so a chance two-flyby chain sits beyond **40 bits**.

**The solar system is, for this purpose, two planets.** At interstellar speeds Mercury, Venus, Earth and Mars cannot change a visitor's velocity by 10% even on a grazing pass, and Jupiter ↔ Saturn accounts for nearly nine chance chains in ten. Merely entering a planet's sphere of influence, by contrast, is common: one in-plane visitor in twenty enters at least one, which is why "it came within the Hill radius" carries almost no information.

**A stellar itinerary is detectable, and precision is the lever.** On a synthetic neighborhood with realistic stellar motions, a visitor planted to come from one nearby star and turn at the Sun for another scores 12–15 bits; one threading far stars scores almost nothing, because the stars' own velocity errors smear where they were a million years ago. Every threefold improvement in stellar velocities buys about 4 bits.

**The stars that matter are mostly not near us now.** A hop is measured where the visitor and the star meet, not where the star is today, and stars outrun slow visitors: of the 2,031 Gaia stars ʻOumuamua could have come from within 20 parsecs of the Sun, four in five are farther than that now, some nearly 200 parsecs away. A search confined to today's neighborhood misses most of the candidates.

**A route has to be scored as one trip.** Scoring where a visitor came from and where it is heading as two independent tests overstates chance routes: a pair of stars on one line through the Sun is hit at both ends by any visitor aimed down that corridor, and counted twice. The score treats a route as origin, Sun, destination: the destination is judged given the arrival and the angle the Sun turned it through, so only the orientation of the turn is free. A path the Sun barely bends earns nothing for a destination straight ahead; a sharp turn that lands on a star is rare, and scores.

**A mismeasured star nearly gave ʻOumuamua a home.** On the first pass through Gaia, 1I scored 4.3 bits for its origin. The star responsible, Gaia DR3 1928116685826146304, lies 42 parsecs away and appeared to have passed 1.3 parsecs from the Sun 63,000 years ago, just where ʻOumuamua then was. Its cataloged radial velocity is 654 ± 7 km/s, against 21 km/s across the sky. All twelve of its Gaia spectra were deblended from a brighter companion 1.1″ away at the same distance, the configuration Gaia's own validation flags for spurious high velocities. Taken at face value the pair would be crossing the Galaxy at 880 km/s, well beyond escape speed and faster than any binary could be thrown without being torn apart. That such a star topped the list is no coincidence: a spurious radial velocity, added to a genuine motion across the sky, points a star's track almost through the Sun, which manufactures exactly the short, cheap hop the score looks at first. The statistics held (the star wins for 0.29% of random visitors at 1I's speed, and 4.3 bits was a one-in-twenty result, reported as such), but a route built on it would have been an artifact. Radial velocities now count only at Gaia's recommended signal-to-noise (≥ 5, which sets aside 6,886) and for stars the Galaxy can hold (under 600 km/s about its center, which removes six, this one among them).

## How it works

**Visit.** For every planet, the score asks what fraction of the planet's possible positions along its orbit would have brought it at least as close as it actually came. That is a p-value. Across the planets, an order statistic asks how surprising it is that *k* planets were all this close, for the best *k*, and corrects for having tried every *k*: one very close pass and several fairly close ones both count, and a near-miss is never diluted by the planets that were nowhere near.

**Assist.** The visitor is flown twice, once among the real planets and once in a universe where all their mass sits in the Sun at the barycenter, and the two outbound asymptotes are compared. For the real objects the change is split into the planets' direct pull and the Sun's own reflex wobble, which shifts every visitor's perihelion a little and is nothing a trajectory could exploit.

**Route.** Inbound and outbound asymptotes are straight lines through a field of moving stars. For the origin, the exact chance that an isotropic arrival at this speed passes each star this closely is computed by integrating over the sphere of possible relative velocities. For the destination, the chance is taken over the circle of departures the Sun's turn allows, given the arrival. Both are convolved with the measurement noise of the miss distance. The look-elsewhere correction weights stars by 1/rank of the hop, the distance from the Sun at which the visitor could meet them, among those it could physically reach, encoding a traveler's preference for short hops.

**The null is the product.** Planets are given random phases; arrivals are swept over impact parameter, and in a second, time-progressed version over heading and epoch with the real planets advancing on one clock. The two agree. Flybys are patched conics, checked against REBOUND N-body integration to within half a percent. Every score is tested against millions of chance visitors to confirm that P(score ≥ x) never exceeds 2⁻ˣ: a score is allowed to understate evidence, never to overstate it.

Orbits come from JPL's Small-Body Database, planet states from JPL Horizons, and stars from Gaia DR3, all fetched on first use. When a fourth interstellar object is discovered, `real_visitors.py 4I` scores it.

## Run

```
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt && mkdir -p out
.venv/bin/python -m pytest -q
.venv/bin/python run.py > out/report.txt          # chance of 1–3 velocity-changing flybys per perihelion ring
.venv/bin/python score.py > out/scores.txt        # visit and assist scores of planned paths
.venv/bin/python route_check.py > out/route.txt   # route score on a synthetic star field
.venv/bin/python real_visitors.py > out/real.txt  # 1I, 2I, 3I against JPL planets and Gaia DR3, and the site's data
.venv/bin/python -m http.server -d site           # the explorer, locally
.venv/bin/python factor3d.py                      # in-plane to 3D factor for one Jupiter flyby
```

Figures are written as SVG to `out/`; the explorer in `site/` is published to GitHub Pages on every push that changes it. Architecture, decisions and conventions are in [AGENTS.md](AGENTS.md).
