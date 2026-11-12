"""Shared helpers for the 复变函数 manim scenes (Manim Community, no LaTeX needed).

All formulas are drawn with plain `Text` + Unicode so the scenes render on a
machine without a TeX installation.
"""
from manim import *  # noqa: F401,F403

FONT = "Noto Sans CJK SC"

# palette (works on the default black background)
C_Z = YELLOW
C_W = TEAL
C_RE = RED_C
C_IM = BLUE_C
C_HL = ORANGE


def T(s: str, size: int = 28, color=WHITE, **kw) -> Text:
    """Chinese-capable text."""
    return Text(s, font=FONT, font_size=size, color=color, **kw)


def plane(x_range=(-4, 4, 1), y_range=(-3.5, 3.5, 1), unit=1.0, **kw) -> ComplexPlane:
    p = ComplexPlane(
        x_range=list(x_range),
        y_range=list(y_range),
        x_length=(x_range[1] - x_range[0]) * unit,
        y_length=(y_range[1] - y_range[0]) * unit,
        background_line_style={"stroke_color": BLUE_E, "stroke_width": 1, "stroke_opacity": 0.6},
        axis_config={"stroke_width": 2},
        **kw,
    )
    return p


def axis_labels(p: ComplexPlane, re="Re", im="Im", size=22):
    xr = p.x_range[1]
    yr = p.y_range[1]
    a = T(re, size, GREY_B).next_to(p.n2p(complex(xr, 0)), DOWN + LEFT * 0.3, buff=0.12)
    b = T(im, size, GREY_B).next_to(p.n2p(complex(0, yr)), RIGHT, buff=0.12)
    return VGroup(a, b)


def tick_numbers(p: ComplexPlane, xs=(), ys=(), size=18):
    g = VGroup()
    for x in xs:
        g.add(T(str(x), size, GREY_B).next_to(p.n2p(complex(x, 0)), DOWN, buff=0.08))
    for y in ys:
        g.add(T(f"{y}i", size, GREY_B).next_to(p.n2p(complex(0, y)), LEFT, buff=0.08))
    return g


def panel(*lines, width=None, size=26, buff=0.22, color=WHITE):
    """A vertical stack of text lines (strings or mobjects)."""
    mobs = [T(l, size, color) if isinstance(l, str) else l for l in lines]
    g = VGroup(*mobs).arrange(DOWN, aligned_edge=LEFT, buff=buff)
    return g


def boxed(m, color=GREY_D):
    r = SurroundingRectangle(m, color=color, buff=0.2, corner_radius=0.1)
    r.set_fill(BLACK, opacity=0.75)
    return VGroup(r, m)
