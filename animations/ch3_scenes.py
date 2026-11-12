"""第三章 复变函数的积分 —— manim 动画源码（LaTeX 排版）."""
import numpy as np

from common import *  # noqa: F401,F403


# --------------------------------------------------------------------------
# 3.1 复积分 = 黎曼和的极限：小向量 f(ζ_k)Δz_k 首尾相接
# --------------------------------------------------------------------------
class C3_RiemannSum(Scene):
    def construct(self):
        title = title_bar(r"复积分 $\displaystyle\int_C f(z)\,\mathrm dz=\lim\sum_k f(\zeta_k)\,\Delta z_k$")
        self.play(FadeIn(title))
        zp = plane((-0.5, 1.5, 0.5), (-0.5, 1.5, 0.5), unit=2.2).shift(LEFT * 3.6 + DOWN * 0.6)
        sp = plane((-1, 0.5, 0.25), (-1, 0.5, 0.25), unit=3.0).shift(RIGHT * 3.4 + DOWN * 0.6)
        zl = CT(r"$z$ 平面：曲线 $C$（单位圆弧，从 $1$ 到 $i$）", 24, GREY_A).next_to(zp, UP, buff=0.1)
        sl = CT(r"累加 $\sum f(\zeta_k)\Delta z_k$", 24, GREY_A).next_to(sp, UP, buff=0.1)
        self.play(Create(zp), Create(sp), FadeIn(zl), FadeIn(sl))
        f = lambda z: z * z
        exact = (1j ** 3 - 1) / 3  # ∫ z^2 dz = (i^3 - 1)/3
        C = ParametricFunction(lambda t: zp.n2p(np.exp(1j * t)), t_range=[0, PI / 2], color=C_Z, stroke_width=4)
        flab = CT(r"$f(z)=z^2$", 30, C_Z).next_to(zp, DOWN, buff=0.15)
        self.play(Create(C), FadeIn(flab))
        ex = Dot(sp.n2p(exact), color=RED, radius=0.08)
        exl = MT(r"\tfrac{-1-i}{3}", 28, RED).next_to(ex, DL, buff=0.08)
        self.play(FadeIn(ex), FadeIn(exl))

        prev = VGroup()
        for n in (3, 6, 12, 40):
            ts = np.linspace(0, PI / 2, n + 1)
            zs = np.exp(1j * ts)
            mids = np.exp(1j * (ts[:-1] + ts[1:]) / 2)  # ζ_k 取弧中点
            chords = VGroup(*[Line(zp.n2p(zs[k]), zp.n2p(zs[k + 1]), color=WHITE, stroke_width=2) for k in range(n)])
            dots = VGroup(*[Dot(zp.n2p(m), radius=0.04, color=C_W) for m in mids])
            terms = [f(mids[k]) * (zs[k + 1] - zs[k]) for k in range(n)]
            partial = np.concatenate([[0], np.cumsum(terms)])
            arrows = VGroup(*[
                Arrow(sp.n2p(partial[k]), sp.n2p(partial[k + 1]), buff=0, stroke_width=3 if n < 20 else 2,
                      max_tip_length_to_length_ratio=0.3, color=interpolate_color(C_W, C_Z, k / max(n - 1, 1)))
                for k in range(n)
            ])
            nl = CT(rf"$n={n}$：和 $\approx {partial[-1].real:.4f}{partial[-1].imag:+.4f}i$", 26).to_edge(DOWN, buff=0.25)
            self.play(FadeOut(prev), run_time=0.3)
            self.play(Create(chords), FadeIn(dots), run_time=0.6)
            self.play(LaggedStart(*[GrowArrow(a) for a in arrows], lag_ratio=0.6 if n < 20 else 0.1), FadeIn(nl), run_time=1.8 if n < 20 else 1.5)
            self.wait(0.6)
            prev = VGroup(chords, dots, arrows, nl)
        res = CT(r"极限 $=\displaystyle\int_C z^2\,\mathrm dz=\frac{z^3}{3}\Big|_1^{i}=\frac{-1-i}{3}$（与路径无关：$z^2$ 解析）", 28, C_Z).to_edge(DOWN, buff=0.25)
        self.play(FadeOut(prev[3]), FadeOut(flab), FadeIn(res))
        self.wait(2.5)


# --------------------------------------------------------------------------
# 3.1 最重要的积分：∮ dz/(z-z0)^n
# --------------------------------------------------------------------------
class C3_OneOverZ(Scene):
    def construct(self):
        title = title_bar(r"$\displaystyle\oint_{|z|=1}\frac{\mathrm dz}{z^n}$：$n=1$ 时每一小段都贡献 $i\,\mathrm d\theta$")
        self.play(FadeIn(title))
        zp = plane((-1.6, 1.6, 0.5), (-1.6, 1.6, 0.5), unit=1.6).shift(LEFT * 3.6 + DOWN * 0.5)
        sp = plane((-1, 2.5, 0.5), (-1.2, 7, 1), unit=0.68).shift(RIGHT * 3.8 + DOWN * 0.5)
        circ = Circle(radius=zp.get_x_unit_size(), color=C_Z).move_to(zp.n2p(0))
        self.play(Create(zp), Create(sp), Create(circ))
        tp = MT(r"2\pi i", 26, RED).next_to(sp.n2p(2 * PI * 1j), RIGHT, buff=0.1)
        tick = Line(sp.n2p(-0.15 + 2 * PI * 1j), sp.n2p(0.15 + 2 * PI * 1j), color=RED)
        self.play(FadeIn(tp), Create(tick))

        def run(n, label, N=24):
            th = np.linspace(0, 2 * PI, N + 1)
            zs = np.exp(1j * th)
            mids = np.exp(1j * (th[:-1] + th[1:]) / 2)
            terms = [(zs[k + 1] - zs[k]) / mids[k] ** n for k in range(N)]
            part = np.concatenate([[0], np.cumsum(terms)])
            lab = CT(label, 28, C_Z).to_edge(DOWN, buff=0.25)
            self.play(FadeIn(lab))
            grp = VGroup()
            for k in range(N):
                dz = Arrow(zp.n2p(zs[k]), zp.n2p(zs[k + 1]), buff=0, color=C_W, stroke_width=3, max_tip_length_to_length_ratio=0.4)
                fdz = Arrow(sp.n2p(part[k]), sp.n2p(part[k + 1]), buff=0, color=C_Z, stroke_width=3, max_tip_length_to_length_ratio=0.35)
                grp.add(dz, fdz)
                self.play(GrowArrow(dz), GrowArrow(fdz), run_time=0.12)
            res = MT(rf"\oint = {part[-1].real:.2f}{part[-1].imag:+.2f}i", 30, WHITE).move_to(RIGHT * 0.7 + UP * 1.2)
            self.play(FadeIn(res))
            self.wait(1.2)
            return VGroup(lab, grp, res)

        g1 = run(1, r"$n=1$：$z=e^{i\theta}$，$\dfrac{\mathrm dz}{z}=\dfrac{ie^{i\theta}\mathrm d\theta}{e^{i\theta}}=i\,\mathrm d\theta$，全部指向 $i$，累加到 $2\pi i$")
        self.play(FadeOut(g1))
        g2 = run(2, r"$n=2$：$\dfrac{\mathrm dz}{z^2}=ie^{-i\theta}\mathrm d\theta$，方向转了一整圈，首尾相接回到 $0$")
        self.play(FadeOut(g2))
        box = boxed(VGroup(
            MT(r"\oint_{|z-z_0|=r}\frac{\mathrm dz}{(z-z_0)^{n}}=\begin{cases}2\pi i,& n=1\\[2pt] 0,& n\neq1\end{cases}", 40),
            CT(r"（与圆心 $z_0$、半径 $r$ 无关）", 26, GREY_B),
        ).arrange(DOWN, buff=0.2))
        box.move_to(DOWN * 0.3)
        self.play(FadeIn(box))
        self.wait(3)


# --------------------------------------------------------------------------
# 3.2 闭路变形原理与复合闭路定理
# --------------------------------------------------------------------------
def fitw(m, w=5.5):
    if m.width > w:
        m.scale_to_fit_width(w)
    return m


class C3_Deformation(Scene):
    def construct(self):
        title = title_bar(r"闭路变形：只要不越过奇点，回路怎么变，积分都不变")
        self.play(FadeIn(title))
        zp = plane((-3, 3, 1), (-2.2, 2.2, 1), unit=1.15).shift(LEFT * 2.4 + DOWN * 0.55)
        self.play(Create(zp))
        a = 0.6 + 0.3j
        sing = Dot(zp.n2p(a), color=RED, radius=0.08)
        sl = MT(r"z_0", 28, RED).next_to(sing, UR, buff=0.05)
        self.play(FadeIn(sing), FadeIn(sl))
        f = lambda z: 1 / (z - a)

        s = ValueTracker(0.0)
        c = ValueTracker(0.0)  # 平移量，用于“越过奇点”

        def shape(t, s_, c_, shift=-1.3 - 0.4j):
            c_ = c_ * shift
            # 大的不规则回路 -> 绕 a 的小圆
            big = (1.9 + 0.35 * np.cos(3 * t) + 0.2 * np.sin(2 * t)) * np.exp(1j * t)
            small = a + 0.35 * np.exp(1j * t)
            return (1 - s_) * big + s_ * small + c_

        def integral(s_, c_, N=600):
            t = np.linspace(0, 2 * PI, N + 1)
            z = shape(t, s_, c_)
            mid = (z[1:] + z[:-1]) / 2
            return np.sum(f(mid) * np.diff(z))

        curve = always_redraw(lambda: ParametricFunction(lambda t: zp.n2p(shape(t, s.get_value(), c.get_value())), t_range=[0, 2 * PI], color=C_Z, stroke_width=4))
        val = always_redraw(lambda: CT(
            r"$\displaystyle\oint_C\frac{\mathrm dz}{z-z_0}\approx %+.3f%+.3fi$" % (integral(s.get_value(), c.get_value()).real, integral(s.get_value(), c.get_value()).imag), 28
        ).set_x(1.35, LEFT).set_y(1.6))
        self.play(Create(curve), FadeIn(val))
        self.wait(0.5)
        cap = fitw(CT(r"连续收缩回路（不碰奇点）：值保持 $2\pi i\approx6.283i$", 26, C_Z)).move_to(UP * 0.6, aligned_edge=LEFT).set_x(1.35, LEFT)
        self.play(FadeIn(cap))
        self.play(s.animate.set_value(1.0), run_time=3.5, rate_func=smooth)
        self.wait(0.6)
        cap2 = fitw(CT(r"把回路移开、不再包围 $z_0$：值变成 $0$（柯西定理）", 26, C_HL)).next_to(cap, DOWN, buff=0.3, aligned_edge=LEFT)
        self.play(FadeIn(cap2))
        self.play(c.animate.set_value(1.0), run_time=3, rate_func=smooth)
        self.wait(1.2)
        self.play(FadeOut(VGroup(cap, cap2, val, curve, sing, sl)))

        # 复合闭路
        title2 = title_bar(r"复合闭路定理：$\displaystyle\oint_C f\,\mathrm dz=\sum_{k}\oint_{C_k}f\,\mathrm dz$（全部逆时针）")
        self.play(FadeTransform(title, title2))
        p1, p2 = -0.9 + 0.2j, 1.0 - 0.3j
        d1, d2 = Dot(zp.n2p(p1), color=RED), Dot(zp.n2p(p2), color=RED)
        C = ParametricFunction(lambda t: zp.n2p((2.3 + 0.2 * np.cos(2 * t)) * np.cos(t) + 1j * 1.75 * np.sin(t)), t_range=[0, 2 * PI], color=C_Z, stroke_width=4)
        C1 = Circle(radius=0.45 * zp.get_x_unit_size(), color=C_W).move_to(zp.n2p(p1))
        C2 = Circle(radius=0.45 * zp.get_x_unit_size(), color=C_W).move_to(zp.n2p(p2))
        labs = VGroup(MT("C", 30, C_Z).next_to(zp.n2p(2.5 + 0.0j), RIGHT, buff=0.05),
                      MT("C_1", 28, C_W).next_to(C1, UP, buff=0.05), MT("C_2", 28, C_W).next_to(C2, UP, buff=0.05))
        self.play(Create(C), FadeIn(d1), FadeIn(d2), FadeIn(labs[0]))
        self.play(Create(C1), Create(C2), FadeIn(labs[1:]))
        reg = fitw(CT(r"在 $C$ 与 $C_1,C_2$ 之间的多连通域内 $f$ 解析", 26)).set_x(1.35, LEFT).set_y(1.4)
        self.play(FadeIn(reg))
        cut1 = Line(zp.n2p(p1 - 0.45), zp.n2p(-2.35 + 0.05j), color=GREY_B)
        cut2 = Line(zp.n2p(p2 + 0.45), zp.n2p(2.45 - 0.2j), color=GREY_B)
        cutl = fitw(CT(r"用割缝连起来 $\Rightarrow$ 单连通；沿割缝来回各走一次，相互抵消", 24, GREY_B)).next_to(reg, DOWN, buff=0.3, aligned_edge=LEFT)
        self.play(Create(cut1), Create(cut2), FadeIn(cutl))
        rule = fitw(CT(r"边界正向：外边界逆时针，内边界顺时针（区域在左手边）", 24, C_HL)).next_to(cutl, DOWN, buff=0.3, aligned_edge=LEFT)
        arrs = VGroup(
            Arrow(C.point_from_proportion(0.12), C.point_from_proportion(0.16), buff=0, color=C_Z),
            Arrow(C1.point_from_proportion(0.30), C1.point_from_proportion(0.22), buff=0, color=C_W),
            Arrow(C2.point_from_proportion(0.30), C2.point_from_proportion(0.22), buff=0, color=C_W),
        )
        self.play(FadeIn(rule), FadeIn(arrs))
        res = fitw(CT(r"$\Rightarrow\ \displaystyle\oint_C f\,\mathrm dz=\oint_{C_1}f\,\mathrm dz+\oint_{C_2}f\,\mathrm dz$", 30, C_Z)).next_to(rule, DOWN, buff=0.35, aligned_edge=LEFT)
        self.play(FadeIn(res))
        self.wait(3)
