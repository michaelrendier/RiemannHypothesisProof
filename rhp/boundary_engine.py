# SPDX-License-Identifier: GPL-3.0-only
"""Boundary Condition Engine — a viewer for every pattern in the five Gemini event-horizon files, the transitions between them, and the ζ pathway.

THREE LAYERS (the nes-viewport model — the menu is a layer OVER the viewport, never part of it):
  1  MAP layer      the world: a GRID of pattern tiles (rows = source files, cols = phases), content cw × ch, scrolled by a Viewport (vw × vh tiles).
                    It draws pattern points from rhp.boundary and records exactly what it drew (`App.frame`).
  2  FRAME          the Viewport (rhp/viewport.py, from the skill): offset clamped after every move; the cursor tile is kept inside it (screen-snap archetype).
  3  OVERLAY layer  scrollbars (DERIVED from the viewport: vbar/hbar), the MENU bar (key → action, built from MENU), the HUD (diagnostics).  It holds NO state of its own:
                    `App.overlay()` is a pure function of the App state, so the bar can never disagree with the map.
Modes: MAP (browse tiles) → STAGE (one pattern big, in motion, parameters live) → TRANSITION (A → B at α with a chosen mapping).
Keys:  arrows move the cursor (the viewport follows) · PgUp/PgDn/Home/End page the viewport · Enter open / run · Esc back · space play/pause · , . step τ (or α) ·
       [ ] choose parameter · - = change it · Tab next target pattern · m next transition mapping · z toggle the ζ ghost on the same clock · s save a PNG.
Run:   python -m rhp.boundary_engine                 interactive window (matplotlib; needs a display)
       python -m rhp.boundary_engine --export DIR    static poster + strips + transition strips + GIF motion (headless)
       python -m rhp.boundary_engine --selftest      drive the menu by key presses and verify the map layer moved as the menu said (headless)
"""
from __future__ import annotations
import os, sys
import numpy as np

from . import boundary as B
from .viewport import Viewport, assert_viewport_sound

BG, INK, ACC, ACC2, DIM, GOLDC = "#0b0b14", "#d8d8e8", "#00e5ff", "#ff4d6d", "#4a4a66", "#c9a84c"
TAU_MIN, TAU_STEP = 0.02, 0.02
MENU = [  # (keys, label, action) — the menu layer is generated from this list; nothing else draws menu items
    ("arrows", "move", None), ("PgUp/PgDn", "page", None), ("Enter", "open/run", "enter"), ("Esc", "back", "back"), ("space", "play/pause", "play"),
    (", .", "step τ/α", None), ("[ ]", "param", None), ("- =", "change", None), ("Tab", "target", "target"), ("m", "mapping", "mapping"), ("z", "ζ ghost", "ghost"),
]
KEYMAP = {"up": "up", "down": "down", "left": "left", "right": "right", "pageup": "pgup", "pagedown": "pgdn", "home": "home", "end": "end", "enter": "enter",
          "escape": "back", "backspace": "back", " ": "play", ",": "tau-", ".": "tau+", "[": "param-", "]": "param+", "-": "val-", "=": "val+", "+": "val+",
          "tab": "target", "m": "mapping", "z": "ghost", "s": "save"}


class App:
    def __init__(self, vw=3, vh=2):
        self.grid = B.GRID; self.rows, self.cols = len(self.grid), max(len(r) for r in self.grid)
        self.vp = Viewport(cw=self.cols, ch=self.rows, vw=vw, vh=vh, page_overlap=0)          # screen-snap: pages are whole screens
        self.cur = (0, 0); self.mode = "map"; self.selected = "O1"; self.tau = None; self.playing = False
        self.params = {}; self.pidx = 0; self.target = "O3"; self.tmode = "polar"; self.alpha = 0.0; self.ghost = False; self.frame = {}; self.dt = 0.02
        self._ensure_visible(); self.refresh()

    # ── the map layer's world ────────────────────────────────────────────────────────────────
    def cell(self, r, c):
        return self.grid[r][c] if 0 <= r < self.rows and 0 <= c < len(self.grid[r]) else None

    def _ensure_visible(self):
        r, c = self.cur
        self.vp.set_offset(x=min(max(self.vp.ox, c - self.vp.vw + 1), c), y=min(max(self.vp.oy, r - self.vp.vh + 1), r))

    def _move(self, dr, dc):
        r, c = self.cur[0] + dr, self.cur[1] + dc
        while 0 <= r < self.rows and 0 <= c < self.cols and self.cell(r, c) is None and (dr or dc):    # skip empty cells in the direction of travel
            r += dr; c += dc
        if self.cell(r, c): self.cur = (r, c); self._ensure_visible()

    def visible_tiles(self):
        out = []
        for r in self.vp.visible_rows():
            for c in range(self.vp.ox, min(self.vp.cw, self.vp.ox + self.vp.vw)):
                out.append(((r, c), self.cell(r, c)))
        return out

    # ── the menu layer: key → action → state change ───────────────────────────────────────────
    def key(self, k):
        a = KEYMAP.get(k, k); m = self.mode
        if a == "save": return a
        if m == "map":
            if a in ("up", "down", "left", "right"): self._move({"up": -1, "down": 1}.get(a, 0), {"left": -1, "right": 1}.get(a, 0))
            elif a == "pgdn": self.vp.page_down(); self._cursor_into_view()
            elif a == "pgup": self.vp.page_up(); self._cursor_into_view()
            elif a == "home": self.cur = (0, 0); self.vp.move_to_top(); self.vp.move_to_left(); self._ensure_visible()
            elif a == "end":
                self.cur = max((r, c) for r in range(self.rows) for c in range(self.cols) if self.cell(r, c)); self.vp.move_to_bottom(); self.vp.move_to_right(); self._ensure_visible()
            elif a == "enter": self.open(self.cell(*self.cur))
        elif m == "stage":
            q = B.PATTERNS[self.selected]; names = list(self._pnames())
            if a == "back": self.mode = "map"; self.playing = False
            elif a == "play": self.playing = not self.playing
            elif a == "tau-": self.tau = max(TAU_MIN, round(self.tau - TAU_STEP, 10))
            elif a == "tau+": self.tau = min(1.0, round(self.tau + TAU_STEP, 10))
            elif a in ("param-", "param+") and names: self.pidx = (self.pidx + (1 if a == "param+" else -1)) % len(names)
            elif a in ("val-", "val+") and names: self._bump(names[self.pidx], 1 if a == "val+" else -1)
            elif a == "target": self._next_target()
            elif a == "mapping": self.tmode = B.TRANSITION_MODES[(B.TRANSITION_MODES.index(self.tmode) + 1) % len(B.TRANSITION_MODES)]
            elif a == "ghost": self.ghost = not self.ghost
            elif a == "enter": self.mode = "transition"; self.alpha = 0.0; self.playing = True
        elif m == "transition":
            if a == "back": self.mode = "stage"; self.playing = False
            elif a == "play": self.playing = not self.playing
            elif a == "tau-": self.alpha = max(0.0, round(self.alpha - .05, 10))
            elif a == "tau+": self.alpha = min(1.0, round(self.alpha + .05, 10))
            elif a == "mapping": self.tmode = B.TRANSITION_MODES[(B.TRANSITION_MODES.index(self.tmode) + 1) % len(B.TRANSITION_MODES)]
            elif a == "target": self._next_target()
            elif a == "enter": self.alpha = 0.0; self.playing = True
        self.refresh(); return a

    def _cursor_into_view(self):
        r = min(max(self.cur[0], self.vp.oy), self.vp.oy + self.vp.vh - 1); c = min(max(self.cur[1], self.vp.ox), self.vp.ox + self.vp.vw - 1)
        if not self.cell(r, c):
            for cc in range(c, -1, -1):
                if self.cell(r, cc): c = cc; break
        self.cur = (r, c)

    def open(self, pid):
        if not pid: return
        q = B.PATTERNS[pid]; self.selected = pid; self.mode = "stage"; self.tau = q.tau0; self.params = dict(q.params); self.pidx = 0; self.playing = False
        self.alpha = 0.0
        if self.target == pid: self._next_target()

    def _pnames(self):
        return [k for k, v in self.params.items() if isinstance(v, (int, float))]

    def _bump(self, name, sgn):
        v = self.params[name]; step = abs(v) * .1 if v else .1
        if name in ("K", "N"): step = max(100, v * .25)
        self.params[name] = v + sgn * step

    def _next_target(self):
        ids = [p for row in self.grid for p in row if p and p != self.selected]
        i = ids.index(self.target) if self.target in ids else -1; self.target = ids[(i + 1) % len(ids)]

    def tick(self, dt=None):
        """Advance the clock one animation step: τ in STAGE, α in TRANSITION.  Returns True if something moved."""
        if not self.playing: return False
        d = self.dt if dt is None else dt
        if self.mode == "stage":
            self.tau += d; self.tau = TAU_MIN if self.tau > 1.0 else self.tau
        elif self.mode == "transition":
            self.alpha = min(1.0, self.alpha + d * 2.5)
            if self.alpha >= 1.0: self.playing = False
        else: return False
        self.refresh(); return True

    # ── what the map layer draws (recorded, so tests can check it) ───────────────────────────────
    def refresh(self):
        if self.mode == "map":
            self.frame = dict(mode="map", panels=[dict(rc=rc, pid=pid, tau=B.PATTERNS[pid].tau0, z=_thumb(pid)) for rc, pid in self.visible_tiles() if pid])
        elif self.mode == "stage":
            z = B.points(self.selected, self.tau, self.params)
            self.frame = dict(mode="stage", panels=[dict(pid=self.selected, tau=self.tau, params=dict(self.params), z=z)], stats=B.turn_stats(z))
            if self.ghost: self.frame["ghost"] = B.points("Z1", self.tau)
        else:
            q = B.PATTERNS[self.selected]
            z, used = B.transition(self.selected, self.target, self.alpha, self.tmode, tau_a=self.tau, tau_b=None)
            self.frame = dict(mode="transition", panels=[dict(pid=f"{self.selected}→{self.target}", alpha=self.alpha, mapping=used, z=z)], stats=B.turn_stats(z))
        return self.frame

    # ── the overlay layer: a pure function of the state ──────────────────────────────────────────
    def overlay(self):
        vb = self.vp.vbar(track_len=20); hb = self.vp.hbar(track_len=20)
        return dict(mode=self.mode, vbar=(vb.pos, vb.length, vb.visible), hbar=(hb.pos, hb.length, hb.visible), offset=(self.vp.ox, self.vp.oy),
                    cursor=self.cur, selected=self.selected, target=self.target, mapping=self.tmode, tau=self.tau, alpha=self.alpha, playing=self.playing,
                    param=(self._pnames()[self.pidx] if self._pnames() else None), ghost=self.ghost,
                    menu=[(k, lab) for k, lab, _ in MENU], row_names=[B.ROW_NAMES[r] for r in self.vp.visible_rows()])


_THUMBS = {}


def _thumb(pid):
    if pid not in _THUMBS:
        z = B.points(pid); _THUMBS[pid] = z if len(z) <= 3000 else z[np.linspace(0, len(z) - 1, 3000).astype(int)]
    return _THUMBS[pid]


# ── rendering (matplotlib).  Layer order: map → overlay ───────────────────────────────────────────────
def render(app, fig):
    fig.clf(); fig.set_facecolor(BG); ov = app.overlay(); fr = app.frame
    fig.text(.02, .965, "BOUNDARY CONDITION ENGINE — the Gemini event-horizon patterns and the ζ pathway", color=GOLDC, fontsize=11, va="top")
    if fr["mode"] == "map":
        W, H = app.vp.vw, app.vp.vh; x0, y0, w, h = .04, .20, .86, .70
        for i, ((r, c), pid) in enumerate(app.visible_tiles()):
            if not pid: continue
            ax = fig.add_axes([x0 + (c - app.vp.ox) * w / W + .006, y0 + (H - 1 - (r - app.vp.oy)) * h / H + .01, w / W - .014, h / H - .075]); _style(ax)
            z = next(p["z"] for p in fr["panels"] if p["rc"] == (r, c)); ax.plot(z.real, z.imag, lw=.4, color=ACC if pid[0] != "Z" else ACC2); ax.set_aspect("equal", "datalim")
            q = B.PATTERNS[pid]; ax.set_title(f"{pid}  {q.label}", color=INK, fontsize=7.5, loc="left")
            if (r, c) == app.cur:
                for sp in ax.spines.values(): sp.set_color(GOLDC); sp.set_linewidth(2.2)
        for j, nm in enumerate(ov["row_names"]):
            fig.text(.04, y0 + (H - 1 - j) * h / H + h / H - .028, nm, color=DIM, fontsize=6.5, va="bottom")
        _bars(fig, ov, x0, y0, w, h)
    else:
        ax = fig.add_axes([.04, .20, .60, .72]); _style(ax); p = fr["panels"][0]; z = p["z"]
        ax.plot(z.real, z.imag, lw=.5, color=ACC); ax.plot(z.real[-1:], z.imag[-1:], "o", color=ACC2, ms=4); ax.set_aspect("equal", "datalim")
        if fr["mode"] == "stage":
            ax.set_title(f"{app.selected}  {B.PATTERNS[app.selected].label}   τ = {app.tau:.3f}", color=INK, fontsize=9, loc="left")
            if app.ghost:
                g = fr["ghost"]; k = np.abs(z).max() / (np.abs(g).max() + 1e-12); ax.plot((g * k).real, (g * k).imag, lw=.4, color=ACC2, alpha=.55)
        else:
            ax.set_title(f"{app.selected} → {app.target}   α = {app.alpha:.2f}   mapping: {p['mapping']}", color=INK, fontsize=9, loc="left")
        a2 = fig.add_axes([.69, .52, .29, .38]); _style(a2)
        if fr["mode"] == "transition" or True:
            al, cur = B.transition_curve(app.selected, app.target, app.tmode if app.tmode != "wheel" or (B.PATTERNS[app.selected].state and B.PATTERNS[app.target].state) else "polar", m=21, n=1500)
            for key, col in (("angular_reversals", ACC2), ("curvature_reversals", ACC), ("cusps", GOLDC)): a2.plot(al, cur[key], color=col, lw=1, label=key)
            a2.axvline(app.alpha, color=INK, lw=.6); a2.legend(fontsize=6, facecolor=BG, labelcolor=INK); a2.set_title(f"turn-arounds along {app.selected}→{app.target} ({app.tmode})", color=INK, fontsize=7, loc="left")
        st = fr["stats"]
        fig.text(.69, .46, f"angular reversals   {st['angular_reversals']}\ncurvature reversals {st['curvature_reversals']}\ncusps               {st['cusps']}\nwinding             {st['winding']:.1f}\n"
                 f"steps               {st['steps']}", color=INK, fontsize=8, family="monospace", va="top")
        if app.selected == "O2" and app.mode == "stage":
            fig.text(.69, .30, f"theory: retrograde iff β·ω > 1 (ω = τ, β = 5 → τ > 0.2)\nfraction of run going backwards = {B.retro_fraction(5, app.tau):.3f}", color=GOLDC, fontsize=7.5, family="monospace", va="top")
    mt = "   ".join(f"[{k}] {lab}" for k, lab in ov["menu"])
    fig.text(.02, .10, mt, color=INK, fontsize=7.5, family="monospace")
    fig.text(.02, .06, f"mode {ov['mode'].upper()}   cursor {ov['cursor']}   offset {ov['offset']}   τ {ov['tau'] if ov['tau'] is None else round(ov['tau'], 3)}   α {ov['alpha']:.2f}   "
             f"target {ov['target']}   mapping {ov['mapping']}   param {ov['param']}   {'▶' if ov['playing'] else '⏸'}   ζ-ghost {'on' if ov['ghost'] else 'off'}",
             color=ACC, fontsize=7.5, family="monospace")
    return fig


def _style(ax):
    ax.set_facecolor(BG); ax.tick_params(colors=DIM, labelsize=6)
    for sp in ax.spines.values(): sp.set_color(DIM)


def _bars(fig, ov, x0, y0, w, h):
    tr = 20
    for (pos, ln, vis), horiz in ((ov["vbar"], False), (ov["hbar"], True)):
        ax = fig.add_axes([x0, y0 - .045, w, .012] if horiz else [x0 + w + .012, y0, .012, h]); ax.set_facecolor(BG); ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values(): sp.set_color(DIM)
        if not vis: ax.text(.5, .5, "", transform=ax.transAxes); continue
        a, b = pos / tr, (pos + ln) / tr                                         # derived from the viewport, never stored
        (ax.axvspan if horiz else ax.axhspan)(a if horiz else 1 - b, b if horiz else 1 - a, color=GOLDC); ax.set_xlim(0, 1); ax.set_ylim(0, 1)


# ── headless drivers ─────────────────────────────────────────────────────────────────────────────
def _agg():
    import matplotlib; matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt; return plt


def sound(vp):
    """assert_viewport_sound sweeps the viewport (it is destructive); run it and put the offset back."""
    ox, oy = vp.ox, vp.oy; assert_viewport_sound(vp); vp.set_offset(x=ox, y=oy)


def selftest(verbose=True):
    """Drive the menu with key presses and verify that the MAP layer moved exactly as the menu said.  Raises AssertionError on the first defect."""
    log = []
    def chk(ok, msg):
        log.append(("PASS" if ok else "FAIL", msg)); assert ok, msg
    a = App(); chk(a.mode == "map" and a.frame["panels"], "starts in MAP with tiles in view")
    # viewport soundness + cursor always visible after every navigation key
    for k in ["right", "right", "right", "down", "down", "down", "left", "pagedown", "pageup", "end", "home", "down", "right"]:
        a.key(k); sound(a.vp)
        r, c = a.cur; chk(a.vp.oy <= r < a.vp.oy + a.vp.vh and a.vp.ox <= c < a.vp.ox + a.vp.vw, f"cursor {a.cur} visible after '{k}' (offset {a.vp.ox},{a.vp.oy})")
    # the overlay is derived: bar thumb must follow the offset
    a.key("home"); t0 = a.overlay()["vbar"]; a.key("end"); t1 = a.overlay()["vbar"]; chk(t0[0] == 0 and t1[0] == 20 - t1[1], "scrollbar thumb follows the viewport (top → bottom)")
    chk(a.overlay() == a.overlay(), "overlay is a pure function of state")
    # map tiles drawn == the maths
    a.key("home")
    for p in a.frame["panels"]: chk(np.allclose(p["z"], B.points(p["pid"])[np.linspace(0, len(B.points(p["pid"])) - 1, min(3000, len(B.points(p["pid"])))).astype(int)]), f"tile {p['pid']} shows its pattern")
    # open the worm via the menu, play, check the motion
    a.key("home"); a.key("right"); chk(B.GRID[a.cur[0]][a.cur[1]] == "O2", "cursor reached O2 by arrows"); a.key("enter")
    chk(a.mode == "stage" and a.frame["panels"][0]["pid"] == "O2", "Enter opened O2 in STAGE")
    a.key("tau-"); a.key("tau-"); t_a = a.tau; z_a = a.frame["panels"][0]["z"]; chk(abs(t_a - (B.PATTERNS["O2"].tau0 - .04)) < 1e-9, f"',' stepped τ to {t_a:.3f}")
    a.key("play"); prev = z_a; moved = []; taus = []
    for _ in range(6):
        a.tick(); z = a.frame["panels"][0]["z"]; taus.append(a.tau); moved.append(float(np.mean(np.abs(z - prev)))); prev = z
        chk(np.allclose(z, B.points("O2", a.tau, a.params)), f"frame at τ={a.tau:.3f} equals the worm formula")
    chk(all(m > 1e-6 for m in moved) and all(np.diff(taus) > 0), f"play moves the map every tick (mean displacement {min(moved):.2e}…{max(moved):.2e}) and τ rises")
    a.key("play"); chk(not a.tick(), "pause freezes the clock")
    # parameter menu changes the map
    z0 = a.frame["panels"][0]["z"]; a.key("param+"); a.key("val+"); chk(not np.allclose(a.frame["panels"][0]["z"], z0), f"[ ] and = changed '{a._pnames()[a.pidx]}' and the map followed")
    # the exact turn-around law at the worm
    a.key("back"); a.key("enter"); a.tau = .15; a.refresh(); lo = a.frame["stats"]["angular_reversals"]; a.tau = .6; a.refresh(); hi = a.frame["stats"]["angular_reversals"]
    chk(lo == 0 and hi == 6, f"worm turns back only for τ > 1/β: reversals {lo} at τ = 0.15 (β·ω = 0.75), {hi} at τ = 0.6 (3 cycles → 6)")
    # transition
    a.key("back"); a.key("back"); a.key("home"); a.key("right"); a.key("enter"); chk(a.selected == "O2", "re-opened O2")
    while a.target != "O3": a.key("target")
    for mode in B.TRANSITION_MODES:
        while a.tmode != mode: a.key("mapping")
        a.key("enter"); chk(a.mode == "transition" and a.alpha == 0.0, f"Enter started TRANSITION O2→O3 [{mode}]")
        z0 = a.frame["panels"][0]["z"]; chk(np.allclose(z0, B.transition("O2", "O3", 0, mode, tau_a=a.tau)[0]), f"[{mode}] α = 0 frame equals the transition formula")
        mid = []
        n = 0
        while a.tick(): n += 1; mid.append(a.frame["panels"][0]["z"])
        z1 = a.frame["panels"][0]["z"]; chk(a.alpha == 1.0 and n >= 5, f"[{mode}] α ran 0→1 in {n} ticks")
        chk(all(np.mean(np.abs(mid[i + 1][:min(len(mid[i]), len(mid[i + 1]))] - mid[i][:min(len(mid[i]), len(mid[i + 1]))])) > 1e-9 for i in range(len(mid) - 1)),
            f"[{mode}] the map changes at every α tick")
        if mode == "wheel":
            chk(np.allclose(z0, B.points("O2", a.tau, a.params), atol=1e-9), "wheel α = 0 is exactly the worm at the current τ")
            chk(np.allclose(z1, B.points("O3"), atol=1e-8), "wheel α = 1 is exactly the rosette")
        a.key("back"); chk(a.mode == "stage", "Esc returned to STAGE")
    chk(True, "all three transition mappings ran through the menu")
    if verbose:
        for s, m in log: print(f"[{s}] {m}")
    return log


def export(outdir):
    plt = _agg(); os.makedirs(outdir, exist_ok=True); from PIL import Image
    a = App(vw=5, vh=5); fig = plt.figure(figsize=(16, 12)); a.refresh(); render(a, fig); fig.savefig(os.path.join(outdir, "poster_all_patterns.png"), dpi=90, facecolor=BG)
    for pid in ("O2", "O3", "W1", "Z1", "R3"):
        fig = plt.figure(figsize=(14, 3.2)); fig.set_facecolor(BG); taus = (.1, .3, .6, 1.0)
        for i, t in enumerate(taus):
            ax = fig.add_subplot(1, 4, i + 1); _style(ax); z = B.points(pid, t); ax.plot(z.real, z.imag, lw=.4, color=ACC); ax.set_aspect("equal", "datalim"); ax.set_title(f"{pid} τ = {t}", color=INK, fontsize=8, loc="left")
        fig.savefig(os.path.join(outdir, f"strip_{pid}.png"), dpi=90, facecolor=BG)
    for pa, pb in (("O2", "O3"), ("O1", "O2"), ("R2", "R3"), ("G1", "Z1"), ("O2", "Z1")):
        for mode in B.TRANSITION_MODES:
            fig = plt.figure(figsize=(14, 3)); fig.set_facecolor(BG)
            for i, al in enumerate((0, .25, .5, .75, 1)):
                ax = fig.add_subplot(1, 5, i + 1); _style(ax); z, used = B.transition(pa, pb, al, mode); ax.plot(z.real, z.imag, lw=.4, color=ACC2); ax.set_aspect("equal", "datalim")
                ax.set_title(f"{pa}→{pb} α={al}  [{used}]", color=INK, fontsize=7, loc="left")
            fig.savefig(os.path.join(outdir, f"transition_{pa}_{pb}_{mode}.png"), dpi=80, facecolor=BG)
    def gif(name, make, n=48):
        frames = []
        for i in range(n):
            fig = plt.figure(figsize=(5, 5), dpi=80); fig.set_facecolor(BG); ax = fig.add_axes([.02, .02, .96, .94]); _style(ax); ax.set_xticks([]); ax.set_yticks([])
            z, title = make(i / (n - 1)); ax.plot(z.real, z.imag, lw=.5, color=ACC); ax.set_aspect("equal", "datalim"); ax.set_title(title, color=INK, fontsize=8, loc="left")
            fig.canvas.draw(); frames.append(Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:, :, :3])); plt.close(fig)
        frames[0].save(os.path.join(outdir, name), save_all=True, append=True and False or False, append_images=frames[1:], duration=60, loop=0)
    gif("motion_worm.gif", lambda u: (B.points("O2", .02 + .98 * u), f"worm  τ = {.02 + .98 * u:.2f}"))
    gif("motion_rosette.gif", lambda u: (B.points("O3", .02 + .98 * u), f"rosette  τ = {.02 + .98 * u:.2f}"))
    gif("motion_transition_wheel_O2_O3.gif", lambda u: (B.transition("O2", "O3", u, "wheel")[0], f"worm → rosette  α = {u:.2f}  [wheel]"))
    gif("motion_zeta_height_sweep.gif", lambda u: (B.points("Z1", u), f"ζ path log P_K, t = {3 + 27 * u:.1f}"))
    return sorted(os.listdir(outdir))


def run():
    import matplotlib.pyplot as plt
    app = App(); fig = plt.figure(figsize=(13, 8.4))
    def draw(): render(app, fig); fig.canvas.draw_idle()
    def on_key(ev):
        act = app.key(ev.key)
        if act == "save": fig.savefig(f"boundary_{app.mode}_{app.selected}.png", dpi=110, facecolor=BG)
        draw()
    def on_timer():
        if app.tick(): draw()
    fig.canvas.mpl_connect("key_press_event", on_key); t = fig.canvas.new_timer(interval=60); t.add_callback(on_timer); t.start(); draw(); plt.show()


if __name__ == "__main__":
    if "--selftest" in sys.argv: _agg(); selftest()
    elif "--export" in sys.argv: print(export(sys.argv[sys.argv.index("--export") + 1]))
    else: run()
