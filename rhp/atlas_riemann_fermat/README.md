# The Riemann–Fermat Atlas — every picture beside the code that draws it

**Which atlas this is.** The Riemann–Fermat Atlas: the pictures and generators for the UNS space that contains The Two Trees (Telperion/Laurelin) and the ζ pathway. It is one part of what 0_RB covers, pointed at one space. It is NOT the Box-Kite (sedenion, 42-vertex) Atlas of ValaQuenta/VAPMIP, and NOT the three Abrikosov Atlases in AbrikosovTree/render/. Registry of all atlases: `ContextPlease/claude/ATLAS_REGISTRY.md`. Renamed 2026-11-09 from `rhp/atlas/`.

Each generator is a module in this directory; running it rewrites its SVG here. The maths it draws lives in `rhp/` (named in each docstring), the tests in `tests/`.
Older figures (`../../figures/`) predate the Riemann–Fermat Atlas and are left where notebooks 01–03 expect them; new pictures go here.

| SVG | Generator | Maths drawn | Notebook |
|---|---|---|---|
| `equator.svg` | `python -m rhp.atlas_riemann_fermat.equator` | `rhp/sphere.py`: the Cayley map w = (s−1)/s, the critical line as the equator, poles s = 1 and s = 0, π factored out of θ(t), the winding the flat circle forgets | 04 §1 |
| `polarisation.svg` | `python -m rhp.atlas_riemann_fermat.polarisation` | `rhp/sphere.py`: forward and mirror wave of each term; ellipse by latitude; the standing-wave cross (σ = ½ and n = √(t/2π)) | 04 §2 |
| `wrap.svg` | `python -m rhp.atlas_riemann_fermat.wrap` | `rhp/wrap.py`: integers wrapped on 6 and 30, twin centres, the thinning 0 spot, twin-wheel beat, the springs' tension | 04 §3 |
| `wheels.svg`, `wheels_animated.svg` | `python -m rhp.atlas_riemann_fermat.wheels` (≈ 1 min) | `rhp/wrap.py`: how far K prime wheels reproduce |ζ| by latitude; the first six wheels in motion (SMIL; open in a browser) | 04 §4 |
| `twins.svg` | `python -m rhp.atlas_riemann_fermat.twins` | `rhp/wrap.py`: sieve × un-sieve plane (π(spf), π(gpf)); the circles' twin constant C₂; Hardy–Littlewood count; twin-gap lengths; the clock–wheel resonance at t = 2πn² | 04 §8 |
| `flow.svg` | `python -m rhp.atlas_riemann_fermat.flow` | `rhp/extinction.py`: sieve survivors and smooth numbers as flows in τ = ln y against Buchstab ω and Dickman ρ; the Laplace dual and the reciprocity s′ρ̂·(1+ω̂) = 1 | 05 §1–2 |
| `flow_player.html` | `python -m rhp.atlas_riemann_fermat.flow_player` | the same flow, interactive: drag τ, watch each n ≤ 2500 die (Telperion) and be born (Laurelin); open in a browser | 05 §1 |
| `net.svg` | `python -m rhp.atlas_riemann_fermat.net` | the six faces (u,τ) (u,σ′) (u,t) (τ,σ′) (τ,t) (σ′,t) of the four coordinates, unfolded like a hypercube net | 05 §4 |
| `nodal.svg` | `python -m rhp.atlas_riemann_fermat.nodal` (≈ 2 min) | `rhp/nodal.py`: Re ξ = 0 and Im ξ = 0 on both hemispheres (the ladder), saddles, heartbeat F(t), the weight |χ|, trivial zeros of ζ | 05 §5 |
| `diverge.svg` | `python -m rhp.atlas_riemann_fermat.diverge` | `rhp/diverge.py`: the divergence detector — control on real ψ, injected off-line quadruple, recovered growth slope, detection limit by reach, frequency scan | 05 §8 |
| `uns_map.svg` | `python -m rhp.atlas_riemann_fermat.uns_map` | `rhp/unsmap.py`: number-free Cayley+stereographic map; lines/circles → circles; true cusp survives; edge-on rotation vs cusp; spiral-law identifier | 06 §1 |
| `chladni.svg` | `python -m rhp.atlas_riemann_fermat.chladni` (≈ 20 s) | `rhp/chladni.py`: even/odd parities, flux ∂σ ln|ξ| (motion), nodal rungs on the sphere, off-line flux pocket, lobe angle | 06 §2 |
| `spirals.svg` | `python -m rhp.atlas_riemann_fermat.spirals` | `rhp/coil.py`: ζ in the cylindrical UNS chart, partial-sum cone/cylinder/funnel, conical helices of the explicit formula, the POE pancake coil with zeros as beads | 06 §3 |
| `heartbeat.svg` | `python -m rhp.atlas_riemann_fermat.heartbeat` | `rhp/heartbeat.py`: the explicit formula's standing wave as counter-rotating circles; area budget of the unit circle; who holds the zero-circle area; mean power on exact ψ; ADD/SCALE/SIGN as spiral similarity and conjugation | 06 / docs |

| `wave.svg` | `python -m rhp.atlas_riemann_fermat.wave` | `rhp/wave.py`: the Heartbeat waveform from the two sets (Telperion/Laurelin mirror), zero-tone rebuild, wind envelope, a factor step against the wind, bandwidth to resolve a prime | docs |

| `primepath_flat.svg`, `primepath_3axis.svg`, `primepath_uns.svg` | `python -m rhp.atlas_riemann_fermat.primepath` | `rhp/primepath.py`: log of the partial Euler product at s = ½ + it for 65 / 10⁴ / 10⁵ / 10⁶ primes, flat, isometric + three orthographic views, and on the Riemann sphere | docs |

Image catalog (every image, panel by panel): `docs/IMAGE_CATALOG.md`; PDF `docs/RHP_Image_Catalog_2026-11-08.pdf` (mirror of `~/Desktop/SpaceClaude/`; rebuild with `python3 docs/build_image_catalog.py`).

Conventions: matplotlib with text kept as text; the animated one is hand-written SVG with plain classes and a dark-mode media query. Labels follow the paper: ESTABLISHED / OURS / THEORETICAL / CONJECTURE / OPEN, with `:CALCULATED`.
A picture here shows a space and a machinery; none of it tests RH.

| `gemini_worm.svg` | `gemini_worm.py` | the Gemini event-horizon sims' curves beside the ζ path (2026-11-09; see `docs/two_trees_in_zeta_and_gemini_worm_2026-11-09.md`) |
