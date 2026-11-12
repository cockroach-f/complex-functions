"""第六章 共形映射 —— manim 动画源码（LaTeX 排版）."""
import numpy as np

from common import *  # noqa: F401,F403


def fitw(m, w):
    if m.width > w:
        m.scale_to_fit_width(w)
    return m


def path_on(p, zs, color=WHITE, sw=3, pad=0.0):
    """把复数采样点画到平面 p 上；超出平面范围的部分断开（多段子路径）。"""
    x0, x1 = p.x_range[0] - pad, p.x_range[1] + pad
    y0, y1 = p.y_range[0] - pad, p.y_range[1] + pad
    ok = np.isfinite(zs) & (zs.real >= x0) & (zs.real <= x1) & (zs.imag >= y0) & (zs.imag <= y1)
    vm = VMobject(stroke_color=color, stroke_width=sw)
    run = []
    runs = []
    for z, good in zip(zs, ok):
        if good:
            run.append(p.n2p(z))
        elif run:
            runs.append(run)
            run = []
    if run:
        runs.append(run)
    runs = [r for r in runs if len(r) > 1]
    if not runs:
        vm.set_points_as_corners([p.n2p(0j), p.n2p(0j)])
        vm.set_stroke(opacity=0)
        return vm
    for k, r in enumerate(runs):
        if k == 0:
            vm.set_points_as_corners(r)
        else:
            vm.start_new_path(r[0])
            vm.add_points_as_corners(r[1:])
    return vm


# --------------------------------------------------------------------------
# 6.1 保角性：f'(z0)≠0 时夹角保持；f'(z0)=0 时失效
# --------------------------------------------------------------------------
class C6_Conformality(Scene):
    def construct(self):
        title = title_bar(r"保角性：$f'(z_0)\neq0$ 时，过 $z_0$ 的两条曲线的夹角在映射下不变")
        self.play(FadeIn(title))
        zp = plane((-0.5, 2.5, 0.5), (-0.5, 2, 0.5), unit=1.6).shift(LEFT * 3.5 + DOWN * 0.5)
        wp = plane((-1.5, 3, 1), (-0.5, 4, 1), unit=0.85).shift(RIGHT * 3.4 + DOWN * 0.5)
        zl = CT(r"$z$ 平面", 26, GREY_A).next_to(zp, UP, buff=0.08)
        wl = CT(r"$w=z^2$", 26, GREY_A).next_to(wp, UP, buff=0.08)
        self.play(Create(zp), Create(wp), FadeIn(zl), FadeIn(wl))
        f = lambda z: z * z
        z0 = 1 + 0.6j
        t = np.linspace(-0.8, 0.8, 200)
        c1 = z0 + t                       # 水平直线
        c2 = z0 + 0.9 * np.exp(1j * (np.pi / 2 + 0.6)) * np.sin(t) + 0.35 * t ** 2  # 弯曲曲线
        g1 = path_on(zp, c1, C_Z, 4)
        g2 = path_on(zp, c2, C_W, 4)
        d0 = Dot(zp.n2p(z0))
        self.play(Create(g1), Create(g2), FadeIn(d0))

        def tangent_arrows(p, ff, scale):
            h = 1e-5
            a1 = (ff(z0 + h) - ff(z0)) / h
            v2 = (c2[101] - c2[99]) / (t[101] - t[99])
            a2 = (ff(z0 + h * v2) - ff(z0)) / h
            w0 = ff(z0)
            A1 = Arrow(p.n2p(w0), p.n2p(w0 + scale * a1 / abs(a1)), buff=0, color=C_Z, stroke_width=5)
            A2 = Arrow(p.n2p(w0), p.n2p(w0 + scale * a2 / abs(a2)), buff=0, color=C_W, stroke_width=5)
            ang = np.angle(a2 / a1)
            arc = Arc(radius=0.45, start_angle=np.angle(a1), angle=ang, arc_center=p.n2p(w0), color=C_HL, stroke_width=4)
            return VGroup(A1, A2, arc), ang

        tz, angz = tangent_arrows(zp, lambda z: z, 0.6)
        lz = MT(rf"\alpha\approx{np.degrees(angz):.1f}^\circ", 28, C_HL).next_to(zp.n2p(z0), UL, buff=0.35)
        self.play(FadeIn(tz), FadeIn(lz))
        h1 = path_on(wp, f(c1), C_Z, 4)
        h2 = path_on(wp, f(c2), C_W, 4)
        w0d = Dot(wp.n2p(f(z0)))
        self.play(TransformFromCopy(g1, h1), TransformFromCopy(g2, h2), TransformFromCopy(d0, w0d), run_time=2.2)
        tw, angw = tangent_arrows(wp, f, 1.0)
        lw = MT(rf"\alpha\approx{np.degrees(angw):.1f}^\circ", 28, C_HL).next_to(wp.n2p(f(z0)), UP, buff=0.75).shift(LEFT * 0.9)
        self.play(FadeIn(tw), FadeIn(lw))
        fp = 2 * z0
        cap = VGroup(
            CT(rf"$f'(z_0)=2z_0={fp.real:.0f}{fp.imag:+.1f}i$：切向量整体旋转 $\arg f'(z_0)\approx{np.degrees(np.angle(fp)):.0f}^\circ$，伸缩 $|f'(z_0)|\approx{abs(fp):.2f}$ 倍", 24),
            CT(r"两条切线转过同一个角 $\Rightarrow$ 夹角不变（大小与方向都不变）", 24, C_HL),
        ).arrange(DOWN, buff=0.12)
        fitw(cap, 13).to_edge(DOWN, buff=0.15)
        self.play(FadeIn(cap))
        self.wait(2)

        # f'(0)=0 的情形
        self.play(FadeOut(VGroup(g1, g2, d0, tz, lz, h1, h2, w0d, tw, lw, cap)))
        r = np.linspace(0, 1.4, 120)
        b = np.pi / 5
        k1 = path_on(zp, r + 0j, C_Z, 4)
        k2 = path_on(zp, r * np.exp(1j * b), C_W, 4)
        arcz = Arc(radius=0.5, start_angle=0, angle=b, arc_center=zp.n2p(0), color=C_HL, stroke_width=4)
        az = MT(r"36^\circ", 26, C_HL).next_to(arcz, RIGHT, buff=0.08).shift(UP * 0.1)
        self.play(Create(k1), Create(k2), Create(arcz), FadeIn(az))
        m1 = path_on(wp, (r + 0j) ** 2, C_Z, 4)
        m2 = path_on(wp, (r * np.exp(1j * b)) ** 2, C_W, 4)
        arcw = Arc(radius=0.5, start_angle=0, angle=2 * b, arc_center=wp.n2p(0), color=C_HL, stroke_width=4)
        aw = MT(r"72^\circ", 26, C_HL).next_to(arcw, RIGHT, buff=0.08).shift(UP * 0.15)
        self.play(TransformFromCopy(k1, m1), TransformFromCopy(k2, m2), run_time=2)
        self.play(Create(arcw), FadeIn(aw))
        cap2 = fitw(CT(r"在 $z_0=0$ 处 $f'(0)=0$：$w=z^2$ 把夹角\textbf{加倍}，保角性失效（临界点）", 26, RED_B), 13).to_edge(DOWN, buff=0.25)
        self.play(FadeIn(cap2))
        self.wait(3)


# --------------------------------------------------------------------------
# 6.3 反演 w=1/z：圆 <-> 直线，内部 <-> 外部
# --------------------------------------------------------------------------
class C6_Inversion(Scene):
    def construct(self):
        title = title_bar(r"$w=\dfrac1z$：过原点的圆 $\leftrightarrow$ 直线，单位圆内部 $\leftrightarrow$ 外部")
        self.play(FadeIn(title))
        zp = plane((-2.5, 2.5, 1), (-2.5, 2.5, 1), unit=1.0).shift(LEFT * 3.5 + DOWN * 0.55)
        wp = plane((-2.5, 2.5, 1), (-2.5, 2.5, 1), unit=1.0).shift(RIGHT * 3.5 + DOWN * 0.55)
        zl = CT(r"$z$ 平面", 26, GREY_A).next_to(zp, UP, buff=0.08)
        wl = CT(r"$w$ 平面", 26, GREY_A).next_to(wp, UP, buff=0.08)
        uc = VGroup(Circle(radius=1, color=GREY_B, stroke_width=2).move_to(zp.n2p(0)),
                    Circle(radius=1, color=GREY_B, stroke_width=2).move_to(wp.n2p(0)))
        self.play(Create(zp), Create(wp), FadeIn(zl), FadeIn(wl), Create(uc))
        inv = lambda z: 1 / z
        t = np.linspace(-60, 60, 4000)
        ys = np.concatenate([-np.logspace(2, -2, 1500), np.logspace(-2, 2, 1500)])
        cols = [C_Z, C_HL, RED_C, PURPLE_B]
        cs = [0.4, 0.7, 1.2, 2.0]
        src = VGroup(*[path_on(zp, c + 1j * ys, cols[k], 3) for k, c in enumerate(cs)])
        tgt = VGroup(*[path_on(wp, inv(c + 1j * ys), cols[k], 3) for k, c in enumerate(cs)])
        cap = fitw(CT(r"竖直线 $\operatorname{Re}z=c$ $\mapsto$ 过原点的圆 $\left|w-\frac1{2c}\right|=\frac1{2c}$", 26), 13).to_edge(DOWN, buff=0.2)
        self.play(Create(src), FadeIn(cap))
        self.play(TransformFromCopy(src, tgt), run_time=2.5)
        self.wait(1.2)
        self.play(FadeOut(src), FadeOut(tgt))

        # 圆 -> 直线
        th = np.linspace(0.002, 2 * np.pi - 0.002, 2000)
        circ = 1 + np.exp(1j * th)
        s1 = path_on(zp, circ, C_W, 4)
        t1 = path_on(wp, inv(circ), C_W, 4)
        cap2 = fitw(CT(r"圆 $|z-1|=1$（过原点）$\mapsto$ 直线 $\operatorname{Re}w=\frac12$；原点 $\mapsto\infty$", 26, C_W), 13).to_edge(DOWN, buff=0.2)
        self.play(Create(s1), Transform(cap, cap2))
        self.play(TransformFromCopy(s1, t1), run_time=2)
        self.wait(1)
        circ2 = 0.5 + 1.2 * np.exp(1j * th)
        s2 = path_on(zp, circ2, C_Z, 4)
        t2 = path_on(wp, inv(circ2), C_Z, 4)
        cap3 = fitw(CT(r"不过原点的圆 $\mapsto$ 圆（\textbf{保圆性}：把直线看成过 $\infty$ 的圆）", 26, C_Z), 13).to_edge(DOWN, buff=0.2)
        self.play(Create(s2), Transform(cap, cap3))
        self.play(TransformFromCopy(s2, t2), run_time=2)
        self.wait(1)
        self.play(FadeOut(VGroup(s1, t1, s2, t2)))

        # 内部 <-> 外部
        pts = [0.45 * np.exp(1j * a) for a in np.linspace(0.3, 2 * np.pi + 0.3, 7)[:-1]]
        dz = VGroup(*[Dot(zp.n2p(p), color=C_Z, radius=0.07) for p in pts])
        dw = VGroup(*[Dot(wp.n2p(1 / p), color=C_Z, radius=0.07) for p in pts])
        cap4 = fitw(CT(r"$|z|<1\mapsto|w|>1$；$\arg w=-\arg z$：反演 $+$ 关于实轴的对称", 26, C_HL), 13).to_edge(DOWN, buff=0.2)
        self.play(FadeIn(dz), Transform(cap, cap4))
        self.play(TransformFromCopy(dz, dw), run_time=2)
        self.wait(3)


# --------------------------------------------------------------------------
# 6.3 上半平面 -> 单位圆：w=(z-i)/(z+i)
# --------------------------------------------------------------------------
class C6_HalfPlaneToDisk(Scene):
    def construct(self):
        title = title_bar(r"$w=\dfrac{z-i}{z+i}$：上半平面 $\to$ 单位圆盘，实轴 $\to$ 单位圆周，$i\to0$")
        self.play(FadeIn(title))
        zp = plane((-3, 3, 1), (-0.5, 3.5, 1), unit=1.0).shift(LEFT * 3.5 + DOWN * 0.7)
        wp = plane((-1.3, 1.3, 0.5), (-1.3, 1.3, 0.5), unit=1.75).shift(RIGHT * 3.6 + DOWN * 0.55)
        zl = CT(r"$z$ 平面：$\operatorname{Im}z>0$", 26, GREY_A).next_to(zp, UP, buff=0.08)
        wl = CT(r"$w$ 平面：$|w|<1$", 26, GREY_A).next_to(wp, UP, buff=0.08)
        self.play(Create(zp), Create(wp), FadeIn(zl), FadeIn(wl))
        T = lambda z: (z - 1j) / (z + 1j)
        fill = Rectangle(width=6, height=3.5, stroke_width=0, fill_color=C_Z, fill_opacity=0.12)
        fill.move_to(zp.n2p(1.75j))
        self.play(FadeIn(fill))
        xs = np.linspace(-200, 200, 6000)
        xs = np.sinh(np.linspace(-6, 6, 3000))
        hs = [0.25, 0.5, 1.0, 2.0, 3.2]
        vs = [-2, -1, 0, 1, 2]
        yy = np.sinh(np.linspace(-6, 6.5, 3000))
        yy = yy[yy > 0]
        Hs = VGroup(*[path_on(zp, xs + 1j * h, C_W, 2.5) for h in hs])
        Vs = VGroup(*[path_on(zp, v + 1j * yy, C_HL, 2.5) for v in vs])
        Ax = path_on(zp, xs + 0j, RED_C, 5)
        self.play(Create(Hs), Create(Vs), Create(Ax), run_time=1.5)
        Hw = VGroup(*[path_on(wp, T(xs + 1j * h), C_W, 2.5) for h in hs])
        Vw = VGroup(*[path_on(wp, T(v + 1j * yy), C_HL, 2.5) for v in vs])
        Aw = path_on(wp, T(xs + 0j), RED_C, 5)
        disk = Circle(radius=wp.get_x_unit_size(), stroke_width=0, fill_color=C_Z, fill_opacity=0.12).move_to(wp.n2p(0))
        self.play(TransformFromCopy(Hs, Hw), TransformFromCopy(Vs, Vw), TransformFromCopy(Ax, Aw), FadeIn(disk), run_time=3)
        pts = [(1j, r"i", r"0"), (0j, r"0", r"-1"), (1 + 0j, r"1", r"-i")]
        g = VGroup()
        for z, a, b in pts:
            d1 = Dot(zp.n2p(z), color=C_Z, radius=0.08)
            d2 = Dot(wp.n2p(T(z)), color=C_Z, radius=0.08)
            g.add(d1, d2, MT(a, 24, C_Z).next_to(d1, UR, buff=0.04), MT(b, 24, C_Z).next_to(d2, UR, buff=0.04))
        self.play(FadeIn(g))
        cap = VGroup(
            CT(r"水平线 $\operatorname{Im}z=c$ $\mapsto$ 与单位圆内切于 $w=1$ 的圆；竖直线 $\mapsto$ 过 $w=1$ 的圆弧", 24),
            CT(r"$i$ 与 $-i$ 关于实轴对称 $\Rightarrow$ 像 $0$ 与 $\infty$ 关于单位圆对称（\textbf{保对称性}）", 24, C_HL),
        ).arrange(DOWN, buff=0.12)
        fitw(cap, 13).to_edge(DOWN, buff=0.15)
        self.play(FadeIn(cap))
        self.wait(3.5)


# --------------------------------------------------------------------------
# 6.4 组合映射：上半单位圆 -> 第一象限 -> 上半平面 -> 单位圆
# --------------------------------------------------------------------------
class C6_Pipeline(Scene):
    def construct(self):
        title = title_bar(r"构造共形映射：上半单位圆盘 $\to$ 第一象限 $\to$ 上半平面 $\to$ 单位圆盘", 32)
        self.play(FadeIn(title))
        P = [
            plane((-1.2, 1.2, 1), (-1.2, 1.2, 1), unit=1.1),
            plane((-0.4, 2.6, 1), (-0.4, 2.6, 1), unit=0.88),
            plane((-2.6, 2.6, 1), (-0.4, 2.6, 1), unit=0.5),
            plane((-1.2, 1.2, 1), (-1.2, 1.2, 1), unit=1.1),
        ]
        xs_c = [-5.4, -1.85, 1.65, 5.35]
        for p, x in zip(P, xs_c):
            p.move_to(RIGHT * x + DOWN * 0.3)
        names = [r"$z$", r"$\zeta=\frac{1+z}{1-z}$", r"$\eta=\zeta^2$", r"$w=\frac{\eta-i}{\eta+i}$"]
        labs = VGroup(*[CT(n, 24, GREY_A).next_to(p, UP, buff=0.1) for n, p in zip(names, P)])
        arrows = VGroup(*[Arrow(P[k].get_right(), P[k + 1].get_left(), buff=0.08, color=GREY_B, stroke_width=3) for k in range(3)])
        maps = [lambda z: z, lambda z: (1 + z) / (1 - z), lambda z: ((1 + z) / (1 - z)) ** 2,
                lambda z: (((1 + z) / (1 - z)) ** 2 - 1j) / (((1 + z) / (1 - z)) ** 2 + 1j)]
        self.play(*[Create(p) for p in P], FadeIn(labs), FadeIn(arrows))

        th = np.linspace(1e-3, np.pi - 1e-3, 1500)
        xr = np.concatenate([np.linspace(-0.999, 0.98, 1500), 1 - np.logspace(-1.7, -7, 1500)])
        rr = np.concatenate([np.linspace(0.0005, 0.98, 1200), 1 - np.logspace(-1.7, -7, 800)])

        def fam(k):
            f = maps[k]
            p = P[k]
            g = VGroup()
            g.add(path_on(p, f(xr + 0j), RED_C, 4))              # 直径 [-1,1]
            g.add(path_on(p, f(np.exp(1j * th)), BLUE_C, 4))      # 上半圆周
            for r in (0.3, 0.6, 0.85):
                g.add(path_on(p, f(r * np.exp(1j * th)), C_Z, 2))
            for a in (np.pi / 6, np.pi / 3, np.pi / 2, 2 * np.pi / 3, 5 * np.pi / 6):
                g.add(path_on(p, f(rr * np.exp(1j * a)), C_W, 2))
            return g

        fills = [
            Sector(radius=P[0].get_x_unit_size(), angle=PI, fill_opacity=0.15, fill_color=C_Z, stroke_width=0).move_arc_center_to(P[0].n2p(0)),
            Polygon(P[1].n2p(0), P[1].n2p(2.6), P[1].n2p(2.6 + 2.6j), P[1].n2p(2.6j), fill_opacity=0.15, fill_color=C_Z, stroke_width=0),
            Polygon(P[2].n2p(-2.6), P[2].n2p(2.6), P[2].n2p(2.6 + 2.6j), P[2].n2p(-2.6 + 2.6j), fill_opacity=0.15, fill_color=C_Z, stroke_width=0),
            Circle(radius=P[3].get_x_unit_size(), fill_opacity=0.15, fill_color=C_Z, stroke_width=0).move_to(P[3].n2p(0)),
        ]
        caps = [
            r"上半单位圆盘：边界为直径（红）与半圆周（蓝），两条边界交于 $\pm1$，夹角 $90^\circ$",
            r"分式线性映射把交点 $-1\mapsto0$、$1\mapsto\infty$：两条圆弧变成从原点出发的两条射线（第一象限）",
            r"幂函数 $\zeta^2$ 把张角 $\frac{\pi}{2}$ 放大为 $\pi$：第一象限 $\to$ 上半平面",
            r"标准映射：上半平面 $\to$ 单位圆盘。复合得 $w=\dfrac{(1+z)^2-i(1-z)^2}{(1+z)^2+i(1-z)^2}$",
        ]
        prev = fam(0)
        self.play(FadeIn(fills[0]), Create(prev), run_time=1.5)
        cap = fitw(CT(caps[0], 24), 13).to_edge(DOWN, buff=0.25)
        self.play(FadeIn(cap))
        self.wait(1.2)
        for k in range(1, 4):
            nxt = fam(k)
            ncap = fitw(CT(caps[k], 24, C_HL if k == 3 else WHITE), 13).to_edge(DOWN, buff=0.25)
            self.play(TransformFromCopy(prev, nxt), FadeIn(fills[k]), Transform(cap, ncap), run_time=2.5)
            self.wait(1.3)
            prev = nxt
        self.wait(3)
