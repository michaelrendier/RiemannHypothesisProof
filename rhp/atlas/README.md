# The atlas — every picture beside the code that draws it

Each generator is a module in this directory; running it rewrites its SVG here. The maths it draws lives in `rhp/` (named in each docstring), the tests in `tests/`.
Older figures (`../../figures/`) predate the atlas and are left where notebooks 01–03 expect them; new pictures go here.

| SVG | Generator | Maths drawn | Notebook |
|---|---|---|---|
| `equator.svg` | `python -m rhp.atlas.equator` | `rhp/sphere.py`: the Cayley map w = (s−1)/s, the critical line as the equator, poles s = 1 and s = 0, π factored out of θ(t), the winding the flat circle forgets | 04 §1 |
| `polarisation.svg` | `python -m rhp.atlas.polarisation` | `rhp/sphere.py`: forward and mirror wave of each term; ellipse by latitude; the standing-wave cross (σ = ½ and n = √(t/2π)) | 04 §2 |
| `wrap.svg` | `python -m rhp.atlas.wrap` | `rhp/wrap.py`: integers wrapped on 6 and 30, twin centres, the thinning 0 spot, twin-wheel beat, the springs' tension | 04 §3 |
| `wheels.svg`, `wheels_animated.svg` | `python -m rhp.atlas.wheels` (≈ 1 min) | `rhp/wrap.py`: how far K prime wheels reproduce |ζ| by latitude; the first six wheels in motion (SMIL; open in a browser) | 04 §4 |

Conventions: matplotlib with text kept as text; the animated one is hand-written SVG with plain classes and a dark-mode media query. Labels follow the paper: ESTABLISHED / OURS / THEORETICAL / CONJECTURE / OPEN, with `:CALCULATED`.
A picture here shows a space and a machinery; none of it tests RH.
