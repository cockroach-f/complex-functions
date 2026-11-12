#!/usr/bin/env python3
"""生成讲义中的静态示意图（matplotlib）。输出到 docs/assets/figs/。

用法：python tools/make_figures.py
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, FancyArrowPatch, Wedge, Polygon

plt.rcParams["font.sans-serif"] = ["Noto Sans CJK SC", "Noto Sans CJK JP", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
OUT = Path(__file__).resolve().parent.parent / "docs" / "assets" / "figs"
OUT.mkdir(parents=True, exist_ok=True)

FILL = "#4db6ac"
EDGE = "#00796b"


def axes(ax, lim=2.6, title=None):
    ax.axhline(0, color="#999", lw=0.8)
    ax.axvline(0, color="#999", lw=0.8)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    if title:
        ax.set_title(title, fontsize=12)


def fig_regions():
    fig, axs = plt.subplots(2, 3, figsize=(12, 8))
    # (a) 邻域 / 去心邻域
    ax = axs[0, 0]
    axes(ax, title="(a) δ 邻域 与 去心邻域")
    ax.add_patch(Circle((-1.1, 0.6), 0.9, fc=FILL, ec=EDGE, ls="--", alpha=0.5))
    ax.plot(-1.1, 0.6, "ko", ms=4)
    ax.text(-1.1, 0.75, "$z_0$", ha="center")
    ax.text(-1.1, -0.55, r"$|z-z_0|<\delta$", ha="center", fontsize=10)
    ax.add_patch(Circle((1.2, 0.6), 0.9, fc=FILL, ec=EDGE, ls="--", alpha=0.5))
    ax.plot(1.2, 0.6, "o", mfc="white", mec="k", ms=6)
    ax.text(1.2, 0.75, "$z_0$", ha="center")
    ax.text(1.2, -0.55, r"$0<|z-z_0|<\delta$", ha="center", fontsize=10)
    ax.text(0, -1.6, "虚线 = 不含边界；空心点 = 挖掉中心", ha="center", fontsize=9, color="#555")

    # (b) 内点、外点、边界点
    ax = axs[0, 1]
    axes(ax, title="(b) 内点 · 外点 · 边界点")
    t = np.linspace(0, 2 * np.pi, 200)
    blob = np.c_[1.3 * np.cos(t) + 0.2 * np.cos(3 * t), 1.1 * np.sin(t)]
    ax.add_patch(Polygon(blob, fc=FILL, ec=EDGE, alpha=0.5))
    for (x, y, r, name) in [(-0.2, 0.1, 0.35, "内点"), (1.95, 1.6, 0.3, "外点"), (1.5, 0.0, 0.3, "边界点")]:
        ax.add_patch(Circle((x, y), r, fc="none", ec="k", lw=1))
        ax.plot(x, y, "ko", ms=3)
        ax.text(x, y - r - 0.28, name, ha="center", fontsize=10)
    ax.text(-0.9, -1.6, "G", fontsize=13, color=EDGE)

    # (c) 圆环
    ax = axs[0, 2]
    axes(ax, title="(c) 1<|z|<2：二连通区域（中间有洞）")
    ax.add_patch(Wedge((0, 0), 2, 0, 360, width=1, fc=FILL, ec=EDGE, ls="--", alpha=0.5))
    th = np.linspace(0, 2 * np.pi, 100)
    ax.plot(1.5 * np.cos(th), 1.5 * np.sin(th), color="#d84315", lw=1.5)
    ax.text(0, -2.45, "红色圆周无法在区域内缩成一点", ha="center", fontsize=9, color="#d84315")

    # (d) 角形区域
    ax = axs[1, 0]
    axes(ax, title="(d) 0<arg z<π/3：角形区域（无界、单连通）")
    ax.add_patch(Wedge((0, 0), 3.5, 0, 60, fc=FILL, ec=EDGE, ls="--", alpha=0.5))
    ax.text(0.55, 0.18, "π/3", fontsize=10)

    # (e) 闭圆盘
    ax = axs[1, 1]
    axes(ax, title="(e) |z−(2+i)|≤1：闭区域（闭圆盘）")
    ax.set_xlim(-0.6, 3.6)
    ax.set_ylim(-1.1, 2.6)
    ax.add_patch(Circle((2, 1), 1, fc=FILL, ec=EDGE, lw=2, alpha=0.5))
    ax.plot(2, 1, "ko", ms=3)
    ax.text(2.05, 1.1, "2+i")
    ax.text(1.5, -0.6, "实线 = 包含边界", fontsize=9, color="#555")

    # (f) 不连通
    ax = axs[1, 2]
    axes(ax, title="(f) |z+1|<0.8 ∪ |z−1|<0.8：不是区域")
    ax.add_patch(Circle((-1.1, 0), 0.8, fc=FILL, ec=EDGE, ls="--", alpha=0.5))
    ax.add_patch(Circle((1.1, 0), 0.8, fc=FILL, ec=EDGE, ls="--", alpha=0.5))
    ax.text(0, -1.5, "开集，但不连通（两块之间连不起来）", ha="center", fontsize=9, color="#555")
    fig.tight_layout()
    fig.savefig(OUT / "ch1_regions.png", dpi=130, bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)


def fig_curves():
    fig, axs = plt.subplots(1, 4, figsize=(13, 3.9))
    t = np.linspace(0, 2 * np.pi, 400)
    data = [
        ("简单 · 闭", np.c_[np.cos(t) + 0.25 * np.cos(2 * t), 0.8 * np.sin(t)]),
        ("简单 · 不闭", np.c_[np.linspace(-1.2, 1.2, 400), 0.5 * np.sin(np.linspace(0, 3 * np.pi, 400))]),
        ("不简单 · 闭（∞字形）", np.c_[np.sin(t), np.sin(t) * np.cos(t)]),
        ("不简单 · 不闭", np.c_[np.linspace(-1.3, 1.3, 400) + 0.6 * np.sin(np.linspace(0, 4 * np.pi, 400)), 0.7 * np.cos(np.linspace(0, 4 * np.pi, 400)) - 0.1 * np.linspace(-1, 1, 400)]),
    ]
    for ax, (name, xy) in zip(axs, data):
        axes(ax, lim=1.5, title=name)
        ax.axhline(0, color="white", lw=2)
        ax.axvline(0, color="white", lw=2)
        ax.plot(xy[:, 0], xy[:, 1], color=EDGE, lw=2)
        k = len(xy) // 6
        ax.annotate("", xy=xy[k + 4], xytext=xy[k], arrowprops=dict(arrowstyle="-|>", color="#d84315", lw=1.5))
    fig.tight_layout()
    fig.savefig(OUT / "ch1_curves.png", dpi=130, bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)


def fig_orientation():
    fig, axs = plt.subplots(1, 2, figsize=(10, 4.2))
    ax = axs[0]
    axes(ax, lim=2.4, title="单连通：边界逆时针 = 区域在左手边")
    th = np.linspace(0, 2 * np.pi, 300)
    ax.fill(1.8 * np.cos(th), 1.5 * np.sin(th), color=FILL, alpha=0.4)
    ax.plot(1.8 * np.cos(th), 1.5 * np.sin(th), color=EDGE, lw=2)
    for a in (0.3, 2.4, 4.4):
        p0 = np.array([1.8 * np.cos(a), 1.5 * np.sin(a)])
        p1 = np.array([1.8 * np.cos(a + 0.25), 1.5 * np.sin(a + 0.25)])
        ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=18, color="#d84315"))
    ax.text(0, 0, "D", fontsize=16, ha="center", color=EDGE)
    ax = axs[1]
    axes(ax, lim=2.4, title="多连通：外边界逆时针，内边界顺时针")
    ax.fill(2 * np.cos(th), 1.7 * np.sin(th), color=FILL, alpha=0.4)
    ax.plot(2 * np.cos(th), 1.7 * np.sin(th), color=EDGE, lw=2)
    for (cx, cy, r) in [(-0.8, 0.3, 0.45), (0.8, -0.3, 0.5)]:
        ax.fill(cx + r * np.cos(th), cy + r * np.sin(th), color="white")
        ax.plot(cx + r * np.cos(th), cy + r * np.sin(th), color=EDGE, lw=2)
        for a in (0.5, 3.6):
            p0 = np.array([cx + r * np.cos(a), cy + r * np.sin(a)])
            p1 = np.array([cx + r * np.cos(a - 0.5), cy + r * np.sin(a - 0.5)])
            ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=14, color="#1565c0"))
    for a in (0.3, 2.4, 4.4):
        p0 = np.array([2 * np.cos(a), 1.7 * np.sin(a)])
        p1 = np.array([2 * np.cos(a + 0.25), 1.7 * np.sin(a + 0.25)])
        ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=18, color="#d84315"))
    ax.text(0, 1.1, "D（三连通）", fontsize=12, ha="center", color=EDGE)
    fig.tight_layout()
    fig.savefig(OUT / "ch1_orientation.png", dpi=130, bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)


def fig_curve_equations():
    fig, axs = plt.subplots(1, 5, figsize=(15, 3.9))
    th = np.linspace(0, 2 * np.pi, 400)
    items = [
        ("|z−i| = 2\n⇒ x²+(y−1)² = 4", lambda ax: ax.plot(2 * np.cos(th), 1 + 2 * np.sin(th), color=EDGE, lw=2), 3.2),
        ("|z+i| = |z−i|\n⇒ y = 0（实轴）", lambda ax: ax.plot([-3, 3], [0, 0], color=EDGE, lw=3), 3.2),
        ("|z−2i| = |z+2|\n⇒ y = −x", lambda ax: ax.plot([-3, 3], [3, -3], color=EDGE, lw=2), 3.2),
        ("|z+1|+|z−1| = 4\n⇒ x²/4 + y²/3 = 1", lambda ax: ax.plot(2 * np.cos(th), np.sqrt(3) * np.sin(th), color=EDGE, lw=2), 3.2),
        ("Re(z²) = 1\n⇒ x² − y² = 1", None, 3.2),
    ]
    for ax, (name, draw, lim) in zip(axs, items):
        axes(ax, lim=lim, title=name)
        if draw:
            draw(ax)
        else:
            s = np.linspace(-1.8, 1.8, 200)
            for sg in (1, -1):
                ax.plot(sg * np.cosh(s), np.sinh(s), color=EDGE, lw=2)
            ax.plot([-3, 3], [-3, 3], "--", color="#aaa", lw=0.8)
            ax.plot([-3, 3], [3, -3], "--", color="#aaa", lw=0.8)
    axs[0].plot(0, 1, "o", color="#d84315", ms=4)
    axs[1].plot([0, 0], [1, -1], "o", color="#d84315", ms=4)
    axs[2].plot([0, -2], [2, 0], "o", color="#d84315", ms=4)
    axs[3].plot([-1, 1], [0, 0], "o", color="#d84315", ms=4)
    fig.tight_layout()
    fig.savefig(OUT / "ch1_curve_equations.png", dpi=120, bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)


if __name__ == "__main__":
    fig_regions()
    fig_curves()
    fig_orientation()
    fig_curve_equations()
    print("figures ->", OUT)
