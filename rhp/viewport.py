# Vendored verbatim from ~/.claude/skills/nes-viewport/assets/viewport.py (author's own skill assets), 2026-11-03.
# Provenance: OURS (methodology skill); unchanged apart from this header.
# SPDX-License-Identifier: GPL-3.0-only
"""Framework-agnostic 2D viewport + scrollbar + the soundness test.

No dependencies. Port the maths verbatim into Android / curses / a canvas; keep
`assert_viewport_sound` as the acceptance gate. Units are "cells" (list rows,
tiles, screens); a cell can be 1 for a text list or e.g. 48 for a 48px row —
`row_px` carries that so the flush check can allow one partial cell.

    vp = Viewport(cw=1, ch=500, vw=1, vh=20)     # a 500-row list in a 20-row pane
    vp.page_down(); vp.move_to_bottom()
    bar = vp.vbar(track_len=18)                  # -> Thumb(pos=.., length=..)
    assert_viewport_sound(vp)
"""
from __future__ import annotations

from dataclasses import dataclass

MIN_THUMB = 1                      # cells; use px (e.g. 24) for a pixel bar


def clamp(v, lo, hi):
    return lo if v < lo else hi if v > hi else v


@dataclass
class Thumb:
    pos: int
    length: int
    visible: bool


class Viewport:
    def __init__(self, cw: int, ch: int, vw: int, vh: int,
                 row_px: int = 1, vh_px: int | None = None,
                 page_overlap: int = 1):
        self.cw, self.ch = int(cw), int(ch)
        self.vw, self.vh = int(vw), int(vh)
        self.row_px = int(row_px)
        self.vh_px = int(vh_px) if vh_px is not None else self.vh * self.row_px
        self.page_overlap = int(page_overlap)
        self.ox = 0
        self.oy = 0
        self.follow_tail = False

    # ── derived ──
    @property
    def max_ox(self) -> int:
        return max(0, self.cw - self.vw)

    @property
    def max_oy(self) -> int:
        return max(0, self.ch - self.vh)

    def first_visible_row(self) -> int:
        return self.oy

    def last_visible_row(self) -> int:
        return min(self.ch - 1, self.oy + self.vh - 1)

    def visible_rows(self) -> range:
        return range(self.first_visible_row(), self.last_visible_row() + 1)

    # ── mutation (always clamp; recompute nothing else — the bar is derived) ──
    def set_offset(self, x: int | None = None, y: int | None = None) -> None:
        if x is not None:
            self.ox = clamp(int(x), 0, self.max_ox)
        if y is not None:
            self.oy = clamp(int(y), 0, self.max_oy)
        if self.follow_tail:
            self.oy = self.max_oy

    def up(self, n=1):        self.set_offset(y=self.oy - n)
    def down(self, n=1):      self.set_offset(y=self.oy + n)
    def left(self, n=1):      self.set_offset(x=self.ox - n)
    def right(self, n=1):     self.set_offset(x=self.ox + n)

    def page_up(self):        self.set_offset(y=self.oy - max(1, self.vh - self.page_overlap))
    def page_down(self):      self.set_offset(y=self.oy + max(1, self.vh - self.page_overlap))
    def page_left(self):      self.set_offset(x=self.ox - max(1, self.vw - self.page_overlap))
    def page_right(self):     self.set_offset(x=self.ox + max(1, self.vw - self.page_overlap))

    def move_to_top(self):    self.set_offset(y=0)
    def move_to_bottom(self): self.set_offset(y=self.max_oy)
    def move_to_left(self):   self.set_offset(x=0)
    def move_to_right(self):  self.set_offset(x=self.max_ox)
    def jump(self, y: int):   self.set_offset(y=y)

    # data changed under us (rows appended / removed / relayout) ──────────────
    def content_resized(self, cw: int | None = None, ch: int | None = None) -> None:
        if cw is not None:
            self.cw = int(cw)
        if ch is not None:
            self.ch = int(ch)
        # re-apply the invariant; follow_tail re-pins to the new bottom
        self.set_offset(x=self.ox, y=self.max_oy if self.follow_tail else self.oy)

    def viewport_resized(self, vw: int | None = None, vh: int | None = None,
                         vh_px: int | None = None) -> None:
        if vw is not None:
            self.vw = int(vw)
        if vh is not None:
            self.vh = int(vh)
        if vh_px is not None:
            self.vh_px = int(vh_px)
        self.set_offset(x=self.ox, y=self.oy)

    # ── scrollbars (pure views of the offset) ──
    def vbar(self, track_len: int, min_thumb: int = MIN_THUMB) -> Thumb:
        return _bar(track_len, self.vh, self.ch, self.oy, self.max_oy, min_thumb)

    def hbar(self, track_len: int, min_thumb: int = MIN_THUMB) -> Thumb:
        return _bar(track_len, self.vw, self.cw, self.ox, self.max_ox, min_thumb)

    def vbar_to_offset(self, track_len: int, thumb_top: int,
                       min_thumb: int = MIN_THUMB) -> None:
        t = self.vbar(track_len, min_thumb)
        span = max(1, track_len - t.length)
        self.set_offset(y=round(self.max_oy * clamp(thumb_top, 0, span) / span))


def _bar(track_len: int, extent: int, total: int, off: int, max_off: int,
         min_thumb: int) -> Thumb:
    if total <= extent or max_off <= 0:
        return Thumb(pos=0, length=track_len, visible=False)
    length = clamp(round(track_len * extent / total), min_thumb, track_len)
    pos = round((track_len - length) * (off / max_off))
    return Thumb(pos=pos, length=length, visible=True)


# ── the acceptance test ────────────────────────────────────────────────────

def assert_viewport_sound(vp: Viewport, track_len: int = 20) -> None:
    """Raise AssertionError on the first defect. `track_len` is the scrollbar
    track in the same units as the thumb (rows for a text bar)."""

    # follow_tail pins the offset on purpose; suspend it for the geometry checks,
    # then restore and verify the pin contract separately.
    was_following = vp.follow_tail
    vp.follow_tail = False

    # 1. CLAMP
    assert vp.max_oy == max(0, vp.ch - vp.vh), (
        f"max_oy {vp.max_oy} != max(0, ch-vh) {max(0, vp.ch - vp.vh)}")
    vp.set_offset(y=-10**9)
    assert vp.oy == 0, f"under-clamp: oy={vp.oy}"
    vp.set_offset(y=10**9)
    assert vp.oy == vp.max_oy, f"over-clamp: oy={vp.oy} max_oy={vp.max_oy}"

    scrolls = vp.ch > vp.vh

    # 2. FLUSH
    vp.move_to_top()
    assert vp.first_visible_row() == 0, f"top not flush: first={vp.first_visible_row()}"
    vp.move_to_bottom()
    if scrolls:
        assert vp.last_visible_row() == vp.ch - 1, (
            f"bottom not flush: last={vp.last_visible_row()} ch-1={vp.ch - 1}")
        # last row's far pixel edge should sit on the viewport's far edge,
        # with at most one partial cell of slack when vh_px isn't a whole # of rows
        rows_shown = vp.last_visible_row() - vp.first_visible_row() + 1
        painted_px = rows_shown * vp.row_px
        slack = vp.vh_px % vp.row_px or 0
        assert vp.vh_px - slack <= painted_px <= vp.vh_px + vp.row_px, (
            f"last row clipped/overhangs: painted {painted_px}px into {vp.vh_px}px")

    # 3. THUMB
    if scrolls:
        vp.move_to_top()
        t0 = vp.vbar(track_len)
        assert t0.visible and t0.pos == 0, f"thumb not at top: {t0}"
        vp.move_to_bottom()
        t1 = vp.vbar(track_len)
        assert t1.pos == track_len - t1.length, (
            f"thumb not at bottom: pos={t1.pos} expected {track_len - t1.length}")
        assert t1.length == clamp(round(track_len * vp.vh / vp.ch), MIN_THUMB, track_len), (
            f"thumb length {t1.length} != proportional")
    else:
        assert not vp.vbar(track_len).visible, "bar should be hidden when ch <= vh"

    # 4. SWEEP
    seen: set[int] = set()
    vp.move_to_top()
    guard = 0
    while True:
        seen |= set(vp.visible_rows())
        if vp.oy >= vp.max_oy:
            break
        before = vp.oy
        vp.page_down()
        assert vp.oy <= vp.max_oy, f"page_down overshoot: oy={vp.oy} > max_oy={vp.max_oy}"
        assert vp.oy > before, "page_down made no progress"
        guard += 1
        assert guard < vp.ch + 2, "sweep did not terminate"
    seen |= set(vp.visible_rows())
    assert seen == set(range(vp.ch)), (
        f"sweep missed rows: {sorted(set(range(vp.ch)) - seen)[:8]}...")

    # 5. FOLLOW-TAIL contract (only if the pane was pinned coming in)
    vp.follow_tail = was_following
    if was_following:
        vp.move_to_bottom()
        assert vp.oy == vp.max_oy, "follow_tail: not pinned to the bottom"
        vp.content_resized(ch=vp.ch + 37)
        assert vp.oy == vp.max_oy, "follow_tail: growth did not re-pin"
        vp.content_resized(ch=vp.ch - 37)
        pinned = vp.oy
        vp.follow_tail = False
        vp.up(3)                                  # user scrolls away -> unpinned
        moved = vp.oy
        assert moved < pinned, "unpin: scrolling up had no effect"
        vp.content_resized(ch=vp.ch + 60)         # new content while unpinned
        assert vp.oy == moved, "unpin: appended content yanked the viewport"
        vp.content_resized(ch=vp.ch - 60)
        vp.follow_tail = True
        vp.move_to_bottom()


if __name__ == "__main__":
    for cw, ch, vw, vh in [(1, 500, 1, 20), (1, 20, 1, 20), (1, 21, 1, 20),
                           (40, 1, 10, 1), (40, 60, 12, 18), (1, 1, 1, 1),
                           (1, 3, 1, 20)]:
        v = Viewport(cw, ch, vw, vh)
        assert_viewport_sound(v)
        print(f"ok  content {cw}x{ch}  view {vw}x{vh}  max_o=({v.max_ox},{v.max_oy})")
    print("all self-tests passed")
