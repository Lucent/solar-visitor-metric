# solar-visitor-metric

**If something built a probe and sent it through our solar system, would its path give it away?**

Three interstellar objects have crossed the solar system in under a decade: 1I/ʻOumuamua in 2017, 2I/Borisov in 2019, and 3I/ATLAS in 2025. Avi Loeb has argued that ʻOumuamua might be an artificial light sail, and that 3I/ATLAS's path is suspiciously well arranged: a retrograde orbit lying almost exactly in the ecliptic, with close passes by Venus, Mars and Jupiter. Arguments like these rest on the word *suspicious*, and that word needs a number: how often would chance produce the same thing?

solar-visitor-metric provides that number. It borrows the idea behind assembly theory (an object is evidence of design when chance cannot plausibly build it) and turns every claim about a trajectory into **bits of surprise**: −log₂ of the probability under an explicit, calibrated null model. Ten bits means one in a thousand. Twenty means one in a million.

## Findings

**None of the three known visitors shows any trace of a planned path.**

| | v∞ | Planet visits | Planets' pull on its outbound path | Star-to-star route |
|---|---|---|---|---|
| 1I/ʻOumuamua | 26.4 km/s | 0.2 bits | 0.13% of v∞ | 0 bits |
| 2I/Borisov | 32.3 km/s | 0 bits | 0.02% | 0.8 bits |
| 3I/ATLAS | 58.0 km/s | 0 bits | 0.12% | 2.6 bits |

3I/ATLAS passed Jupiter at 0.358 AU on 16 March 2026, almost exactly one Hill radius, and Mars at 0.194 AU on 3 October 2025. That looks pointed until you ask how often a randomly placed Jupiter comes that close: **3% of the time**, and Mars 19%. Allow for the fact that eight planets were available to be "visited" and the surprise drops to zero. The Jupiter pass changed 3I's velocity by about a thousandth, a hundred times short of anything a mission designer would call a gravity assist. Traced back and forward against 1,850 Gaia DR3 stars within 20 parsecs, none of the three visitors came from or is heading for any of them more closely than chance routinely manages: one random visitor in ten at 3I's speed scores a better route.

The one genuinely unusual fact about 3I is the one Loeb leads with. Its orbit lies within 5° of the ecliptic, which an isotropic arrival does about once in 275 tries: **8 bits**, a real coincidence, but a single one among the many features that could have been checked.

**A planned path would stand out, and by a wide margin.** Chance gives two velocity-changing flybys (≥ 10% each) to roughly one in three million in-plane visitors at ʻOumuamua's speed (**22 bits**), and one in a hundred million at 3I's (**27 bits**). A chance chain falls off roughly as the fourth power of speed, because fast objects barely bend. Real paths flown in a full N-body simulation score accordingly:

| Planned path at 26 km/s | Visits | Assist | Combined |
|---|---|---|---|
| Earth and Mars at 1.5 million km, then a Jupiter assist | 13.0 bits | 10.3 bits | **17 bits** |
| Earth and Mars at 150,000 km, then a Jupiter assist | 21.1 bits | 10.3 bits | **24 bits** |
| Jupiter assist alone | 4.0 bits | 10.3 bits | 10 bits |

A targeted Jupiter → Saturn chain flown in REBOUND misses Saturn by 5.2 AU when Jupiter's mass is switched off: every significant chain is causally dependent on its first flyby, so the second encounter is no coincidence. And all of these are in-plane numbers. In three dimensions each flyby must also hit the planet's orbital plane to within its target distance, a further factor of about 2 × 10⁻³ per flyby by direct Monte Carlo, so a chance two-flyby chain sits beyond **40 bits**.

**The solar system is, for this purpose, two planets.** At interstellar speeds Mercury, Venus, Earth and Mars cannot change a visitor's velocity by 10% even on a grazing pass, and Jupiter ↔ Saturn accounts for nearly nine chance chains in ten. Merely entering a planet's sphere of influence, by contrast, is common: one in-plane visitor in twenty enters at least one, which is why "it came within the Hill radius" carries almost no information.

**A stellar itinerary is detectable, and precision is the lever.** On a synthetic neighborhood with realistic stellar motions, a visitor planted to pass exactly through two nearby stars scores 17–21 bits; one threading stars 15 parsecs out scores almost nothing, because the stars' own velocity errors smear where they were a million years ago. Every threefold improvement in stellar radial velocities buys about 6 bits. Most stars near the Sun move faster than ʻOumuamua does relative to it, so a slow visitor physically cannot have come from most of them, and the score counts only the stars it could have reached.

## How it works

**Visit.** For every planet, the score asks what fraction of the planet's possible positions along its orbit would have brought it at least as close as it actually came. That is a p-value. Across the planets, an order statistic asks how surprising it is that *k* planets were all this close, for the best *k*, and corrects for having tried every *k*: one very close pass and several fairly close ones both count, and a near-miss is never diluted by the planets that were nowhere near.

**Assist.** The visitor is flown twice, once among the real planets and once in a universe where all their mass sits in the Sun at the barycenter, and the two outbound asymptotes are compared. For the real objects the change is split into the planets' direct pull and the Sun's own reflex wobble, which shifts every visitor's perihelion a little and is nothing a trajectory could exploit.

**Route.** Inbound and outbound asymptotes are straight lines through a field of moving stars. For each star, the exact chance that an isotropic visitor at this speed passes it this closely is computed by integrating over the sphere of possible relative velocities, then convolved with the measurement noise of the miss distance. The look-elsewhere correction weights stars by 1/rank of distance among those the visitor could physically reach, encoding a traveler's preference for short hops.

**The null is the product.** Planets are given random phases; arrivals are swept over impact parameter, and in a second, time-progressed version over heading and epoch with the real planets advancing on one clock. The two agree. Flybys are patched conics, checked against REBOUND N-body integration to within half a percent. Every score is tested against millions of chance visitors to confirm that P(score ≥ x) never exceeds 2⁻ˣ: a score is allowed to understate evidence, never to overstate it.

Orbits come from JPL's Small-Body Database, planet states from JPL Horizons, and stars from Gaia DR3, all fetched on first use. When a fourth interstellar object is discovered, `real_visitors.py 4I` scores it.

## Run

```
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt && mkdir -p out
.venv/bin/python -m pytest -q
.venv/bin/python run.py > out/report.txt          # chance of 1–3 velocity-changing flybys per perihelion ring
.venv/bin/python score.py > out/scores.txt        # visit and assist scores of planned paths
.venv/bin/python route_check.py > out/route.txt   # route score on a synthetic star field
.venv/bin/python real_visitors.py > out/real.txt  # 1I, 2I, 3I against JPL planets and Gaia DR3
.venv/bin/python factor3d.py                      # in-plane to 3D factor for one Jupiter flyby
```

Figures are written as SVG to `out/`. Architecture, decisions and conventions are in [AGENTS.md](AGENTS.md).
