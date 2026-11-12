"""第九章 拉普拉斯变换 —— manim 动画源码（LaTeX 排版）."""
import numpy as np

from common import *  # noqa: F401,F403


def fitw(m, w):
    if m.width > w:
        m.scale_to_fit_width(w)
    return m


def ax(xr, yr, xl, yl, **kw):
    return Axes(
        x_range=list(xr), y_range=list(yr), x_length=xl, y_length=yl,
        axis_config={"stroke_width": 2, "include_tip": True, "tip_length": 0.18, "include_ticks": False},
        y_axis_config={"include_tip": False},
        **kw,
    )


def splane(xr=(-3, 3, 1), yr=(-3, 3, 1), unit=0.8):
    p = plane(xr, yr, unit=unit)
    return p


def s_labels(p, size=24):
    a = MT(r"\sigma", size, GREY_B).next_to(p.n2p(complex(p.x_range[1], 0)), DOWN + LEFT * 0.3, buff=0.12)
    b = MT(r"j\omega", size, GREY_B).next_to(p.n2p(complex(0, p.y_range[1])), RIGHT, buff=0.12)
    return VGroup(a, b)


def cross(p, z, color=RED, s=0.11):
    c = p.n2p(z)
    return VGroup(Line(c + (UL) * s, c + (DR) * s, color=color, stroke_width=5),
                  Line(c + (UR) * s, c + (DL) * s, color=color, stroke_width=5))


# --------------------------------------------------------------------------
# 9.1 拉普拉斯变换 = 乘衰减因子 e^{-σt} 后的傅里叶变换；收敛域是右半平面
# --------------------------------------------------------------------------
class C9_DampingROC(Scene):
    def construct(self):
        title = title_bar(r"$\mathcal L[f]=\mathcal F\!\left[f(t)u(t)e^{-\sigma t}\right]$：乘上衰减因子，让增长的信号也可积")
        self.play(FadeIn(title))
        A = ax((0, 6, 1), (-0.1, 3.2, 1), 5.8, 4.2).move_to(LEFT * 3.4 + DOWN * 0.6)
        self.play(Create(A))
        self.add(MT("t", 24, GREY_B).next_to(A.x_axis.get_end(), DOWN, buff=0.1))
        a0 = 0.5
        orig = A.plot(lambda t: np.exp(a0 * t), x_range=[0, np.log(3.2) / a0], color=GREY_B, stroke_width=2.5)
        ol = MT(r"f(t)=e^{0.5t}", 28, GREY_B).move_to(A.c2p(1.4, 2.9))
        self.play(Create(orig), FadeIn(ol))
        note0 = CT(r"$\int_0^{\infty}|f(t)|\,\mathrm dt=\infty$：没有古典傅里叶变换", 24, GREY_B).next_to(A, UP, buff=0.12)
        self.play(FadeIn(note0))
        sig = ValueTracker(0.0)

        def damped():
            k = a0 - sig.get_value()
            xmax = 6 if k <= 0 else min(6, np.log(3.2) / k)
            col = C_W if sig.get_value() > a0 else RED_C
            return A.plot(lambda t: np.exp(k * t), x_range=[0, xmax, 0.02], color=col, stroke_width=4)

        D = always_redraw(damped)
        dl = MT(r"f(t)e^{-\sigma t}", 28, C_W).move_to(A.c2p(4.6, 2.3))
        P = splane((-2, 2.5, 1), (-2.2, 2.2, 1), unit=1.0).move_to(RIGHT * 3.5 + DOWN * 0.4)
        self.play(Create(P), FadeIn(s_labels(P)))
        pole = cross(P, a0)
        pl = MT(r"s=0.5", 22, RED_C).next_to(P.n2p(a0), DOWN + RIGHT * 0.3, buff=0.1)
        self.play(FadeIn(pole), FadeIn(pl))
        line = always_redraw(lambda: DashedLine(P.n2p(complex(sig.get_value(), -2.2)), P.n2p(complex(sig.get_value(), 2.2)),
                                                color=C_W if sig.get_value() > a0 else RED_C, stroke_width=3))
        sl = always_redraw(lambda: VGroup(MT(r"\sigma=", 28), DecimalNumber(sig.get_value(), num_decimal_places=2, font_size=28))
                           .arrange(RIGHT, buff=0.08).next_to(P, UP, buff=0.15))
        self.add(D, line, sl)
        self.play(FadeIn(dl))
        self.play(sig.animate.set_value(0.35), run_time=2)
        bad = CT(r"$\sigma\le0.5$：$f(t)e^{-\sigma t}$ 仍不衰减，积分发散", 24, RED_C).to_edge(DOWN, buff=0.35)
        self.play(FadeIn(bad))
        self.wait(0.8)
        self.play(sig.animate.set_value(1.5), run_time=3)
        good = CT(r"$\sigma>0.5$：指数衰减，积分收敛 $\Rightarrow F(s)=\displaystyle\int_0^{\infty}f(t)e^{-st}\,\mathrm dt$ 存在", 24, C_W).to_edge(DOWN, buff=0.35)
        self.play(FadeTransform(bad, good))
        roc = Polygon(P.n2p(complex(a0, -2.2)), P.n2p(complex(2.5, -2.2)), P.n2p(complex(2.5, 2.2)), P.n2p(complex(a0, 2.2)),
                      stroke_width=0, fill_color=C_W, fill_opacity=0.22)
        rl = CT(r"收敛域 $\operatorname{Re}s>0.5$", 24, C_W).move_to(P.n2p(complex(1.5, 1.6)))
        self.play(FadeIn(roc), FadeIn(rl))
        self.play(sig.animate.set_value(0.9), run_time=1.5)
        fin = VGroup(
            MT(r"F(s)=\frac{1}{s-0.5},\quad \operatorname{Re}s>0.5", 30, C_Z),
            CT(r"$F(s)$ 在收敛半平面内解析；唯一的奇点（极点）在边界左侧", 22, C_Z),
        ).arrange(DOWN, buff=0.15)
        fin.next_to(P, DOWN, buff=0.12)
        fitw(fin, 6.0)
        self.play(FadeOut(good), FadeIn(fin))
        self.wait(3)


# --------------------------------------------------------------------------
# 9.4 极点位置 ↔ 时间响应：左半平面衰减、虚轴等幅、右半平面发散
# --------------------------------------------------------------------------
class C9_PoleResponse(Scene):
    def construct(self):
        title = title_bar(r"极点 $s=\sigma_0\pm j\omega_0$ 决定 $\mathcal L^{-1}$ 的形状：$e^{\sigma_0 t}\cos\omega_0 t$")
        self.play(FadeIn(title))
        P = splane((-2, 2, 1), (-3.4, 3.4, 1), unit=0.85).move_to(LEFT * 4.3 + DOWN * 0.45)
        lhp = Polygon(P.n2p(-2 - 3.4j), P.n2p(0 - 3.4j), P.n2p(0 + 3.4j), P.n2p(-2 + 3.4j), stroke_width=0, fill_color=GREEN, fill_opacity=0.16)
        rhp = Polygon(P.n2p(0 - 3.4j), P.n2p(2 - 3.4j), P.n2p(2 + 3.4j), P.n2p(0 + 3.4j), stroke_width=0, fill_color=RED, fill_opacity=0.16)
        self.play(FadeIn(lhp), FadeIn(rhp), Create(P), FadeIn(s_labels(P)))
        self.add(CT("稳定", 22, GREEN).move_to(P.n2p(-1.2 + 3.0j)), CT("不稳定", 22, RED_C).move_to(P.n2p(1.1 + 3.0j)))
        A = ax((0, 10, 1), (-2.3, 2.3, 1), 7.6, 4.6).move_to(RIGHT * 2.55 + DOWN * 0.45)
        self.play(Create(A))
        self.add(MT("t", 24, GREY_B).next_to(A.x_axis.get_end(), DOWN, buff=0.1))
        sx, sw = ValueTracker(-0.5), ValueTracker(0.0)
        poles = always_redraw(lambda: VGroup(cross(P, complex(sx.get_value(), sw.get_value()), C_Z),
                                             cross(P, complex(sx.get_value(), -sw.get_value()), C_Z)))

        def resp():
            a, w = sx.get_value(), sw.get_value()
            ts = np.linspace(0, 10, 500)
            ys = np.exp(a * ts) * np.cos(w * ts)
            pts = []
            for t, y in zip(ts, ys):
                if abs(y) > 2.3:
                    break
                pts.append(A.c2p(t, y))
            m = VMobject(color=C_Z, stroke_width=3.5)
            if len(pts) > 1:
                m.set_points_smoothly(pts) if len(pts) < 4 else m.set_points_as_corners(pts)
            return m

        R = always_redraw(resp)
        lab = always_redraw(lambda: MT(rf"\sigma_0={sx.get_value():.2f},\ \omega_0={sw.get_value():.2f}", 26, C_Z).next_to(A, UP, buff=0.1))
        self.add(poles, R, lab)
        cap = None
        seq = [
            ((-0.5, 0.0), r"负实轴上的极点：单调衰减 $e^{\sigma_0 t}$"),
            ((-0.3, 2.5), r"左半平面共轭极点：衰减振荡（欠阻尼）"),
            ((0.0, 2.5), r"虚轴上的极点：等幅振荡（临界）"),
            ((0.12, 2.5), r"右半平面：振幅指数增长——不稳定"),
            ((0.25, 0.0), r"正实轴：单调发散"),
            ((-0.8, 1.2), r"离虚轴越远衰减越快；离实轴越远振荡越快"),
        ]
        for k, ((a, w), text) in enumerate(seq):
            c = CT(text, 26, C_HL).to_edge(DOWN, buff=0.3)
            if k == 0:
                self.play(FadeIn(c), run_time=0.6)
            else:
                self.play(sx.animate.set_value(a), sw.animate.set_value(w), FadeTransform(cap, c), run_time=2.2)
            cap = c
            self.wait(1.6)
        fin = CT(r"系统稳定 $\iff$ $H(s)$ 的极点全在左半平面 $\operatorname{Re}s<0$", 28, C_Z).to_edge(DOWN, buff=0.3)
        self.play(FadeTransform(cap, fin))
        self.wait(3)


# --------------------------------------------------------------------------
# 9.3 反演积分：Bromwich 直线 + 左半大圆弧（约当引理）⇒ 留数之和
# --------------------------------------------------------------------------
class C9_Bromwich(Scene):
    def construct(self):
        title = title_bar(r"反演积分 $f(t)=\dfrac{1}{2\pi j}\displaystyle\int_{\beta-j\infty}^{\beta+j\infty}F(s)e^{st}\,\mathrm ds=\sum_k\operatorname{Res}\!\left[F(s)e^{st},s_k\right]$", 30)
        self.play(FadeIn(title))
        P = splane((-4.2, 2.2, 1), (-3.3, 3.3, 1), unit=0.95).move_to(LEFT * 3.3 + DOWN * 0.5)
        self.play(Create(P), FadeIn(s_labels(P)))
        poles = VGroup(cross(P, -1), cross(P, -2))
        pls = VGroup(MT("-1", 22, RED_C).next_to(P.n2p(-1), DOWN, buff=0.15), MT("-2", 22, RED_C).next_to(P.n2p(-2), DOWN, buff=0.15))
        F = MT(r"F(s)=\frac{1}{(s+1)(s+2)}", 32, C_Z).to_edge(RIGHT, buff=0.6).shift(UP * 2.2)
        self.play(FadeIn(poles), FadeIn(pls), FadeIn(F))
        beta = 1.0
        R = ValueTracker(1.2)
        Lline = always_redraw(lambda: Line(P.n2p(complex(beta, -R.get_value())), P.n2p(complex(beta, R.get_value())), color=C_W, stroke_width=4).add_tip(tip_length=0.18))
        bl = MT(r"\operatorname{Re}s=\beta", 24, C_W).next_to(P.n2p(complex(beta, 2.6)), RIGHT, buff=0.08)
        t1 = CT(r"积分路径：收敛域内、所有奇点右侧的竖直线", 24).next_to(F, DOWN, buff=0.4)
        t1.set_x(0.3 + t1.width / 2)
        self.play(Create(Lline), FadeIn(bl), FadeIn(t1))
        self.wait(0.8)
        arc = always_redraw(lambda: ParametricFunction(lambda th: P.n2p(beta + R.get_value() * np.exp(1j * th)), t_range=[PI / 2, 3 * PI / 2], color=C_HL, stroke_width=4))
        al = always_redraw(lambda: MT(r"C_R", 26, C_HL).next_to(P.n2p(beta - R.get_value() * 0.72 + 0.72j * R.get_value()), UL, buff=0.05))
        t2 = CT(r"$t>0$：向\textbf{左}补半圆 $C_R$，$|e^{st}|=e^{t\operatorname{Re}s}$ 在左边很小", 24, C_HL)
        fitw(t2, 6.4).next_to(t1, DOWN, buff=0.3).align_to(t1, LEFT)
        self.play(Create(arc), FadeIn(al), FadeIn(t2))
        self.play(R.animate.set_value(3.2), run_time=3)
        t3 = CT(r"约当引理：$R\to\infty$ 时 $\displaystyle\int_{C_R}F(s)e^{st}\,\mathrm ds\to0$", 24, C_HL)
        fitw(t3, 6.4).next_to(t2, DOWN, buff=0.3).align_to(t1, LEFT)
        self.play(FadeIn(t3), *[Indicate(p, color=C_Z) for p in poles])
        r1 = MT(r"\operatorname{Res}[F e^{st},-1]=\frac{e^{-t}}{-1+2}=e^{-t}", 28)
        r2 = MT(r"\operatorname{Res}[F e^{st},-2]=\frac{e^{-2t}}{-2+1}=-e^{-2t}", 28)
        r3 = MT(r"f(t)=e^{-t}-e^{-2t}\quad(t>0)", 34, C_Z)
        g = VGroup(r1, r2, r3).arrange(DOWN, aligned_edge=LEFT, buff=0.25)
        fitw(g, 6.0).next_to(t3, DOWN, buff=0.4).align_to(t1, LEFT)
        self.play(FadeIn(r1))
        self.play(FadeIn(r2))
        box = SurroundingRectangle(r3, color=C_Z, buff=0.1)
        self.play(FadeIn(r3), Create(box))
        self.wait(1.5)
        # t<0：向右闭合
        arcR = ParametricFunction(lambda th: P.n2p(beta + 1.15 * np.exp(1j * th)), t_range=[-PI / 2, PI / 2], color=GREY_B, stroke_width=3)
        arcR = DashedVMobject(arcR, num_dashes=24)
        t4 = CT(r"$t<0$：改向右补圆弧，内部无奇点 $\Rightarrow f(t)=0$（因果性）", 24, GREY_B)
        fitw(t4, 6.4).to_edge(DOWN, buff=0.3).align_to(t1, LEFT)
        self.play(Create(arcR), FadeIn(t4))
        self.wait(3)


# --------------------------------------------------------------------------
# 9.1/9.4 |H(s)| 曲面：极点处隆起；沿虚轴切一刀就是频率响应 |H(jω)|
# --------------------------------------------------------------------------
def max_den(d):
    return d if abs(d) > 1e-6 else 1e-6


class C9_SurfaceFreqResponse(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=62 * DEGREES, theta=-35 * DEGREES, zoom=0.9)
        title = title_bar(r"$|H(s)|$，$H(s)=\dfrac{1}{(s+0.5)^2+4}$：极点 $s=-0.5\pm2j$ 处隆起", 30)
        self.add_fixed_in_frame_mobjects(title)
        axes = ThreeDAxes(x_range=[-2, 1.5, 1], y_range=[-3.5, 3.5, 1], z_range=[0, 1.6, 0.5], x_length=3.5 * 1.3, y_length=7 * 0.9, z_length=2.6,
                          axis_config={"include_ticks": False, "stroke_width": 2})
        H = lambda s: 1 / max_den((s + 0.5) ** 2 + 4)
        zc = 1.5

        def surf(u, v):
            val = min(abs(H(complex(u, v))), zc)
            return axes.c2p(u, v, val)

        S = Surface(surf, u_range=[-2, 1.5], v_range=[-3.5, 3.5], resolution=(42, 70), fill_opacity=0.75, stroke_width=0.2,
                    stroke_color=GREY_B, checkerboard_colors=[TEAL_D, TEAL_E])
        S.set_fill_by_value(axes=axes, colorscale=[(BLUE_E, 0), (TEAL, 0.4), (YELLOW, 1.0), (ORANGE, 1.5)], axis=2)
        lx = MT(r"\sigma", 30).move_to(axes.c2p(1.75, 0, 0))
        ly = MT(r"j\omega", 30).move_to(axes.c2p(0, 3.9, 0))
        self.add_fixed_orientation_mobjects(lx, ly)
        self.play(FadeIn(title), Create(axes), FadeIn(lx), FadeIn(ly))
        self.play(Create(S), run_time=3)
        cap = CT(r"曲面高度 $=|H(s)|$；在极点处趋于 $\infty$（图中截断）", 26, C_W).to_edge(DOWN, buff=0.5)
        self.add_fixed_in_frame_mobjects(cap)
        self.play(FadeIn(cap))
        self.begin_ambient_camera_rotation(rate=0.12)
        self.wait(3)
        curve = ParametricFunction(lambda v: axes.c2p(0, v, abs(H(1j * v))), t_range=[-3.5, 3.5], color=RED, stroke_width=6)
        jw = Line3D(axes.c2p(0, -3.5, 0), axes.c2p(0, 3.5, 0), color=RED, thickness=0.02)
        cap2 = CT(r"沿虚轴 $s=j\omega$ 切开：$|H(j\omega)|$ 就是频率响应（傅里叶变换）", 26, RED_B).to_edge(DOWN, buff=0.5)
        self.add_fixed_in_frame_mobjects(cap2)
        self.play(FadeOut(cap), FadeIn(cap2), Create(jw), Create(curve), run_time=2.5)
        self.wait(3)
        cap3 = CT(r"极点越靠近虚轴，$|H(j\omega)|$ 在 $\omega\approx\omega_0$ 处的谐振峰越尖", 26, C_Z).to_edge(DOWN, buff=0.5)
        self.add_fixed_in_frame_mobjects(cap3)
        self.play(FadeOut(cap2), FadeIn(cap3))
        self.wait(4)
        self.stop_ambient_camera_rotation()


# --------------------------------------------------------------------------
# 9.2 卷积：翻转、平移、相乘、积分；L[f1*f2]=F1F2
# --------------------------------------------------------------------------
class C9_Convolution(Scene):
    def construct(self):
        title = title_bar(r"卷积 $(f_1*f_2)(t)=\displaystyle\int_0^{t}f_1(\tau)f_2(t-\tau)\,\mathrm d\tau$：翻转、平移、相乘、求面积", 30)
        self.play(FadeIn(title))
        W = 1.5  # 矩形脉冲宽度
        f1 = lambda x: np.exp(-x) if x >= 0 else 0.0
        f2 = lambda x: 1.0 if 0 <= x <= W else 0.0

        def y(t):
            if t <= 0:
                return 0.0
            if t <= W:
                return 1 - np.exp(-t)
            return np.exp(-(t - W)) - np.exp(-t)

        A = ax((-1.6, 6, 1), (0, 1.3, 1), 10.5, 2.2).move_to(UP * 1.25 + LEFT * 0.4)
        B = ax((0, 6, 1), (0, 1.0, 1), 8.7, 2.0).move_to(DOWN * 2.15 + RIGHT * 0.5)
        self.play(Create(A), Create(B))
        self.add(MT(r"\tau", 26, GREY_B).next_to(A.x_axis.get_end(), DOWN, buff=0.1),
                 MT("t", 26, GREY_B).next_to(B.x_axis.get_end(), DOWN, buff=0.1))
        g1 = A.plot(lambda x: np.exp(-x), x_range=[0, 6, 0.02], color=C_W, stroke_width=4)
        l1 = MT(r"f_1(\tau)=e^{-\tau}", 28, C_W).move_to(A.c2p(4.3, 0.95))
        self.play(Create(g1), FadeIn(l1))
        t = ValueTracker(-0.8)

        def win():
            tt = t.get_value()
            a, b = tt - W, tt
            pts = [A.c2p(max(a, -1.6), 0), A.c2p(max(a, -1.6), 1), A.c2p(b, 1), A.c2p(b, 0)]
            return VMobject(color=C_Z, stroke_width=4).set_points_as_corners(pts)

        G2 = always_redraw(win)
        l2 = always_redraw(lambda: MT(r"f_2(t-\tau)", 26, C_Z).next_to(A.c2p(t.get_value() - W / 2, 1.0), UP, buff=0.08))
        cap = CT(r"$f_2$ 是宽 $1.5$ 的矩形脉冲；把它\textbf{翻转}成 $f_2(-\tau)$，再\textbf{平移} $t$", 24, C_Z).to_edge(DOWN, buff=0.05).shift(UP * 3.55)
        cap.next_to(A, DOWN, buff=0.15)
        self.play(Create(G2), FadeIn(l2), FadeIn(cap))

        def area():
            tt = t.get_value()
            a, b = max(0.0, tt - W), tt
            if b <= 0.0 or b - a < 1e-3:
                return VMobject()
            return A.get_area(g1, x_range=[a, b], color=C_HL, opacity=0.55)

        AR = always_redraw(area)
        trace = always_redraw(lambda: B.plot(y, x_range=[0, max(t.get_value(), 0.001), 0.02], color=C_HL, stroke_width=4))
        dot = always_redraw(lambda: Dot(B.c2p(max(t.get_value(), 0), y(t.get_value())), color=C_HL, radius=0.06))
        lb = MT(r"y(t)=(f_1*f_2)(t)=\text{阴影面积}", 26, C_HL, tex_template=CTEX).move_to(B.c2p(4.3, 0.85))
        self.add(AR, trace, dot)
        self.play(FadeIn(lb))
        self.play(t.animate.set_value(W), run_time=3.5, rate_func=linear)
        c2 = CT(r"$0<t\le1.5$：面积 $=\int_0^t e^{-\tau}\,\mathrm d\tau=1-e^{-t}$（充电）", 24, C_Z).move_to(cap)
        self.play(FadeTransform(cap, c2))
        self.play(t.animate.set_value(5.6), l1.animate.set_opacity(0), run_time=5, rate_func=linear)
        c3 = CT(r"$t>1.5$：面积 $=\int_{t-1.5}^{t}e^{-\tau}\,\mathrm d\tau=e^{-(t-1.5)}-e^{-t}$（放电）", 24, C_Z).move_to(cap)
        self.play(FadeTransform(c2, c3))
        self.wait(1)
        fin = MT(r"\mathcal L[f_1*f_2]=F_1(s)F_2(s)=\frac{1}{s+1}\cdot\frac{1-e^{-1.5s}}{s}", 28, C_Z)
        fin = boxed(fin).move_to(B.c2p(4.4, 0.82))
        self.play(FadeOut(lb), FadeIn(fin))
        self.wait(3)
