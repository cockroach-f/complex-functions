"""第四章 级数 —— manim 动画源码（LaTeX 排版）."""
import numpy as np

from common import *  # noqa: F401,F403


def fitw(m, w):
    if m.width > w:
        m.scale_to_fit_width(w)
    return m


# --------------------------------------------------------------------------
# 4.2 几何级数的部分和：|z|<1 收敛到 1/(1-z)，|z|>1 发散
# --------------------------------------------------------------------------
class C4_GeometricPartialSums(Scene):
    def construct(self):
        title = title_bar(r"几何级数 $\displaystyle\sum_{n=0}^\infty z^n$：部分和 $S_N=1+z+\cdots+z^N$ 走向何方？")
        self.play(FadeIn(title))
        zp = plane((-1.5, 1.5, 0.5), (-1.5, 1.5, 0.5), unit=1.3).shift(LEFT * 4.6 + DOWN * 0.6)
        sp = plane((-1, 3, 1), (-1, 3, 1), unit=1.3).shift(RIGHT * 2.0 + DOWN * 0.6)
        zl = CT(r"$z$ 的位置", 24, GREY_A).next_to(zp, UP, buff=0.1)
        sl = CT(r"部分和 $S_N$ 的轨迹", 24, GREY_A).next_to(sp, UP, buff=0.1)
        unit = Circle(radius=zp.get_x_unit_size(), color=GREY_B, stroke_width=2).move_to(zp.n2p(0))
        self.play(Create(zp), Create(sp), Create(unit), FadeIn(zl), FadeIn(sl))

        def show(z, N, cap, col):
            zd = Dot(zp.n2p(z), color=col)
            zt = MT("z", 28, col).next_to(zd, UR, buff=0.05)
            cap_m = fitw(CT(cap, 26, col), 13).to_edge(DOWN, buff=0.2)
            self.play(FadeIn(zd), FadeIn(zt), FadeIn(cap_m))
            S = np.concatenate([[0], np.cumsum(z ** np.arange(N + 1))])
            pts = [sp.n2p(s) for s in S]
            arrows = VGroup()
            for k in range(N + 1):
                if abs(S[k + 1].real) > 3.3 or abs(S[k + 1].imag) > 3.3 or S[k + 1].real < -1.3 or S[k + 1].imag < -1.3:
                    break
                arrows.add(Arrow(pts[k], pts[k + 1], buff=0, stroke_width=3, max_tip_length_to_length_ratio=0.25,
                                 color=interpolate_color(col, WHITE, min(k / 20, 0.7))))
            extra = VGroup()
            if abs(z) < 1:
                lim = 1 / (1 - z)
                ld = Dot(sp.n2p(lim), color=RED, radius=0.09)
                ll = MT(r"\frac{1}{1-z}", 28, RED).next_to(ld, RIGHT, buff=0.1)
                extra.add(ld, ll)
                self.play(FadeIn(extra))
            self.play(LaggedStart(*[GrowArrow(a) for a in arrows], lag_ratio=0.85), run_time=min(0.25 * len(arrows), 5))
            self.wait(1.2)
            return VGroup(zd, zt, cap_m, arrows, extra)

        g = show(0.5 + 0.5j, 25, r"$|z|=0.71<1$：$S_N=\dfrac{1-z^{N+1}}{1-z}\to\dfrac1{1-z}=1+i$", C_Z)
        self.play(FadeOut(g))
        g = show(0.92 * np.exp(0.9j), 60, r"$|z|=0.92<1$：绕着极限点盘旋，越转越近（收敛较慢）", C_W)
        self.play(FadeOut(g))
        g = show(1.08 * np.exp(0.6j), 40, r"$|z|=1.08>1$：$|z^N|\to\infty$，通项不趋于 $0$，发散", C_HL)
        self.play(FadeOut(g))
        box = boxed(VGroup(
            MT(r"\sum_{n=0}^{\infty}z^n=\frac{1}{1-z}\quad(|z|<1)", 40),
            CT(r"$|z|\ge1$ 时发散；收敛域是一个圆盘——幂级数的“收敛圆”", 26, GREY_B),
        ).arrange(DOWN, buff=0.2)).move_to(DOWN * 0.4)
        self.play(FadeIn(box))
        self.wait(3)


# --------------------------------------------------------------------------
# 4.3 收敛半径 = 到最近奇点的距离：1/(1+x^2) 的谜题
# --------------------------------------------------------------------------
class C4_RadiusOfConvergence(Scene):
    def construct(self):
        title = title_bar(r"为什么 $\dfrac{1}{1+x^2}$ 的泰勒级数只在 $|x|<1$ 收敛？")
        self.play(FadeIn(title))
        ax = Axes(x_range=[-2.2, 2.2, 1], y_range=[-1, 2, 1], x_length=6.2, y_length=4.2,
                  axis_config={"stroke_width": 2, "include_numbers": False}).shift(LEFT * 3.4 + DOWN * 0.7)
        nums = VGroup(*[MT(str(v), 22, GREY_B).next_to(ax.c2p(v, 0), DOWN, buff=0.1) for v in (-2, -1, 1, 2)])
        self.play(Create(ax), FadeIn(nums))
        g = ax.plot(lambda x: 1 / (1 + x * x), x_range=[-2.2, 2.2], color=C_Z, stroke_width=4)
        gl = MT(r"\frac{1}{1+x^2}", 30, C_Z).next_to(ax.c2p(1.5, 0.5), UP, buff=0.1)
        self.play(Create(g), FadeIn(gl))
        cap = fitw(CT(r"$\dfrac1{1+x^2}=1-x^2+x^4-x^6+\cdots$：在实轴上处处光滑，为何 $|x|\ge1$ 就发散？", 26), 13).to_edge(DOWN, buff=0.2)
        self.play(FadeIn(cap))

        def taylor(N):
            return lambda x: sum((-1) ** k * x ** (2 * k) for k in range(N + 1))

        cur = None
        for N in (1, 2, 3, 5, 8, 12, 18):
            gr = ax.plot(taylor(N), x_range=[-1.5, 1.5, 0.01], color=C_W, stroke_width=3, use_smoothing=False)
            gr = gr.copy()
            # 截掉超出视窗的部分
            pts = [p for p in gr.points if -1.0 <= ax.p2c(p)[1] <= 2.0]
            lab = MT(rf"N={2*N}", 26, C_W).next_to(ax, UP, buff=0.05).shift(RIGHT * 1.5)
            curve = VMobject(color=C_W, stroke_width=3)
            if len(pts) > 4:
                xs = np.linspace(-1.5, 1.5, 600)
                ys = taylor(N)(xs)
                seg = [ax.c2p(x, y) for x, y in zip(xs, ys) if -1.0 <= y <= 2.0]
                curve.set_points_as_corners(seg)
            grp = VGroup(curve, lab)
            if cur is None:
                self.play(Create(curve), FadeIn(lab), run_time=0.8)
            else:
                self.play(Transform(cur, grp), run_time=0.6)
            cur = cur or grp
            self.wait(0.2)
        bad = VGroup(DashedLine(ax.c2p(-1, -1), ax.c2p(-1, 2), color=RED), DashedLine(ax.c2p(1, -1), ax.c2p(1, 2), color=RED))
        self.play(Create(bad))

        zp = plane((-2, 2, 1), (-2, 2, 1), unit=1.15).shift(RIGHT * 3.8 + DOWN * 0.5)
        zl = CT(r"在复平面上看", 24, GREY_A).next_to(zp, UP, buff=0.1)
        self.play(Create(zp), FadeIn(zl))
        poles = VGroup(Dot(zp.n2p(1j), color=RED, radius=0.09), Dot(zp.n2p(-1j), color=RED, radius=0.09))
        pl = VGroup(MT("i", 26, RED).next_to(poles[0], RIGHT, buff=0.08), MT("-i", 26, RED).next_to(poles[1], RIGHT, buff=0.08))
        cap2 = fitw(CT(r"$\dfrac1{1+z^2}$ 在 $z=\pm i$ 处有奇点（分母为零）", 26, RED), 13).to_edge(DOWN, buff=0.2)
        self.play(FadeIn(poles), FadeIn(pl), Transform(cap, cap2))
        disc = Circle(radius=zp.get_x_unit_size(), color=C_Z, fill_opacity=0.2, stroke_width=3).move_to(zp.n2p(0))
        r = Line(zp.n2p(0), zp.n2p(1j), color=C_HL, stroke_width=4)
        rl = MT("R=1", 26, C_HL).next_to(r, LEFT, buff=0.08)
        self.play(GrowFromCenter(disc), Create(r), FadeIn(rl))
        seg = Line(zp.n2p(-1), zp.n2p(1), color=C_Z, stroke_width=7)
        self.play(Create(seg))
        cap3 = fitw(CT(r"收敛圆最大只能“撑到”最近的奇点：$R=$ 展开中心到最近奇点的距离", 26, C_Z), 13).to_edge(DOWN, buff=0.2)
        self.play(Transform(cap, cap3))
        self.wait(1.5)
        # 换一个展开中心
        c0 = 1.0
        disc2 = Circle(radius=np.sqrt(2) * zp.get_x_unit_size(), color=C_W, fill_opacity=0.15, stroke_width=3).move_to(zp.n2p(c0))
        r2 = Line(zp.n2p(c0), zp.n2p(1j), color=C_W, stroke_width=4)
        cd = Dot(zp.n2p(c0), color=C_W)
        cap4 = fitw(CT(r"在 $z_0=1$ 展开：$R=|1-i|=\sqrt2$（奇点不在实轴上）", 26, C_W), 13).to_edge(DOWN, buff=0.2)
        self.play(FadeIn(cd), Create(r2), GrowFromCenter(disc2), Transform(cap, cap4))
        self.wait(3)


# --------------------------------------------------------------------------
# 4.4 洛朗级数：同一个函数，不同圆环，不同展开式
# --------------------------------------------------------------------------
class C4_LaurentAnnuli(Scene):
    def construct(self):
        title = title_bar(r"$f(z)=\dfrac{1}{(z-1)(z-2)}=\dfrac{1}{1-z}-\dfrac{1}{2-z}$：奇点把平面切成圆环")
        self.play(FadeIn(title))
        zp = plane((-3, 3, 1), (-3, 3, 1), unit=0.95).shift(LEFT * 3.9 + DOWN * 0.55)
        self.play(Create(zp))
        u = zp.get_x_unit_size()
        sing = VGroup(Dot(zp.n2p(1), color=RED, radius=0.08), Dot(zp.n2p(2), color=RED, radius=0.08))
        sl = VGroup(MT("1", 24, RED).next_to(sing[0], DOWN, buff=0.08), MT("2", 24, RED).next_to(sing[1], DOWN, buff=0.08))
        self.play(FadeIn(sing), FadeIn(sl))
        circ = VGroup(DashedLine(ORIGIN, ORIGIN))  # placeholder
        c1 = DashedVMobject(Circle(radius=u, color=GREY_B).move_to(zp.n2p(0)), num_dashes=40)
        c2 = DashedVMobject(Circle(radius=2 * u, color=GREY_B).move_to(zp.n2p(0)), num_dashes=60)
        self.play(Create(c1), Create(c2))

        X = 1.3  # 右侧文字左边界
        W = 5.6

        def region(mob, lines, col):
            txt = VGroup(*[fitw(l, W) for l in lines]).arrange(DOWN, aligned_edge=LEFT, buff=0.28)
            txt.set_x(X, LEFT).set_y(0.0)
            self.play(FadeIn(mob))
            self.play(LaggedStart(*[FadeIn(t) for t in txt], lag_ratio=0.5), run_time=2.2)
            self.wait(2.0)
            self.play(FadeOut(mob), FadeOut(txt))

        r1 = Circle(radius=u, color=C_Z, fill_opacity=0.35, stroke_width=0).move_to(zp.n2p(0))
        region(r1, [
            CT(r"\textbf{(1)} $|z|<1$：$|z|<1,\ |z/2|<1$", 28, C_Z),
            MT(r"\frac{1}{1-z}=\sum_{n=0}^{\infty}z^n,\quad \frac{1}{2-z}=\frac12\sum_{n=0}^{\infty}\frac{z^n}{2^n}", 30),
            MT(r"f(z)=\sum_{n=0}^{\infty}\Big(1-\frac{1}{2^{n+1}}\Big)z^n", 32, C_Z),
            CT(r"只有正幂：就是泰勒级数", 24, GREY_B),
        ], C_Z)
        r2 = Annulus(inner_radius=u, outer_radius=2 * u, color=C_W, fill_opacity=0.35, stroke_width=0).move_to(zp.n2p(0))
        region(r2, [
            CT(r"\textbf{(2)} $1<|z|<2$：$|1/z|<1,\ |z/2|<1$", 28, C_W),
            MT(r"\frac{1}{1-z}=-\frac1z\cdot\frac{1}{1-\frac1z}=-\sum_{n=1}^{\infty}\frac{1}{z^{n}}", 30),
            MT(r"f(z)=-\sum_{n=1}^{\infty}\frac{1}{z^{n}}-\sum_{n=0}^{\infty}\frac{z^n}{2^{n+1}}", 32, C_W),
            CT(r"正、负幂都有：洛朗级数", 24, GREY_B),
        ], C_W)
        r3 = Annulus(inner_radius=2 * u, outer_radius=3.0 * u, color=C_HL, fill_opacity=0.3, stroke_width=0).move_to(zp.n2p(0))
        region(r3, [
            CT(r"\textbf{(3)} $|z|>2$：$|1/z|<1,\ |2/z|<1$", 28, C_HL),
            MT(r"\frac{1}{2-z}=-\frac1z\cdot\frac{1}{1-\frac2z}=-\sum_{n=0}^{\infty}\frac{2^n}{z^{n+1}}", 30),
            MT(r"f(z)=\sum_{n=2}^{\infty}\frac{2^{n-1}-1}{z^{n}}", 32, C_HL),
            CT(r"只有负幂（且从 $z^{-2}$ 开始）", 24, GREY_B),
        ], C_HL)
        # 以 z=1 为中心
        c3 = DashedVMobject(Circle(radius=u, color=GREY_B).move_to(zp.n2p(1)), num_dashes=40)
        self.play(FadeOut(c1), FadeOut(c2), Create(c3))
        r4 = Annulus(inner_radius=0.04, outer_radius=u, color=PURPLE_B, fill_opacity=0.4, stroke_width=0).move_to(zp.n2p(1))
        region(r4, [
            CT(r"\textbf{(4)} $0<|z-1|<1$（以奇点 $z=1$ 为中心的去心邻域）", 28, PURPLE_B),
            MT(r"f(z)=-\frac{1}{z-1}\cdot\frac{1}{1-(z-1)}=-\sum_{n=-1}^{\infty}(z-1)^{n}", 32, PURPLE_B),
            CT(r"$(z-1)^{-1}$ 的系数 $c_{-1}=-1$：第 5 章的“留数”", 26, C_HL),
        ], PURPLE_B)
        box = boxed(CT(r"同一函数在不同圆环内的洛朗展开\textbf{不同}；在同一圆环内展开式\textbf{唯一}", 28)).move_to(DOWN * 0.5)
        self.play(FadeIn(box))
        self.wait(3)
