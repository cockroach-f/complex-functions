"""第二章 解析函数 —— manim 动画源码（LaTeX 排版）."""
import numpy as np

from common import *  # noqa: F401,F403


# --------------------------------------------------------------------------
# 2.1 C-R 方程的几何意义：可导 = 局部只做“旋转 + 伸缩”
# --------------------------------------------------------------------------
class C2_CauchyRiemannGeometry(Scene):
    def construct(self):
        title = title_bar(r"复可导 $\Longleftrightarrow$ 局部看，$f$ 只做旋转 + 伸缩")
        self.play(FadeIn(title))
        zp = plane((-0.5, 2.5, 0.5), (-1, 2, 0.5), unit=1.55).shift(LEFT * 3.6 + DOWN * 0.6)
        wp = plane((-1.5, 3.5, 1), (-1.5, 3.5, 1), unit=0.9).shift(RIGHT * 3.4 + DOWN * 0.6)
        zl = CT(r"$z$ 平面", 28, GREY_A).next_to(zp, UP, buff=0.1)
        wl = CT(r"$w$ 平面", 28, GREY_A).next_to(wp, UP, buff=0.1)
        self.play(Create(zp), Create(wp), FadeIn(zl), FadeIn(wl), run_time=1.2)

        z0 = 1 + 0.5j
        h = 0.55
        # 小网格方块
        def grid_square(c, s, n=4):
            g = VGroup()
            for k in range(n + 1):
                t = -s / 2 + s * k / n
                g.add(Line(zp.n2p(c + complex(-s / 2, t)), zp.n2p(c + complex(s / 2, t)), stroke_width=2, color=GREY_B))
                g.add(Line(zp.n2p(c + complex(t, -s / 2)), zp.n2p(c + complex(t, s / 2)), stroke_width=2, color=GREY_B))
            for m in g:
                m.insert_n_curves(20)
            return g

        sq = grid_square(z0, h)
        ex = Arrow(zp.n2p(z0), zp.n2p(z0 + 0.5), buff=0, color=C_RE, stroke_width=5, max_tip_length_to_length_ratio=0.25)
        ey = Arrow(zp.n2p(z0), zp.n2p(z0 + 0.5j), buff=0, color=C_IM, stroke_width=5, max_tip_length_to_length_ratio=0.25)
        lx = MT(r"\Delta z=\Delta x", 24, C_RE).next_to(ex, DOWN, buff=0.08)
        ly = MT(r"\Delta z=i\Delta y", 24, C_IM).next_to(ey, LEFT, buff=0.08)
        z0d = Dot(zp.n2p(z0), color=C_Z)
        self.play(Create(sq), GrowArrow(ex), GrowArrow(ey), FadeIn(lx), FadeIn(ly), FadeIn(z0d))

        cases = [
            (r"f(z)=z^2", lambda z: z * z, r"像仍是小正方形：红箭头逆时针转 $90^\circ$ 恰为蓝箭头", r"满足 C-R：$u_x=v_y,\ u_y=-v_x$", C_Z),
            (r"f(z)=\bar z", lambda z: np.conj(z), r"像被翻转：蓝箭头在红箭头的顺时针方向", r"$u_x=1\neq v_y=-1$，不满足 C-R", RED_C),
            (r"f(z)=x+2iy", lambda z: z.real + 2j * z.imag, r"像被压扁成长方形：两方向伸缩率不同", r"$u_x=1\neq v_y=2$，不满足 C-R", RED_C),
        ]
        prev = VGroup()
        for tex, f, cap1, cap2, col in cases:
            def F(p, f=f):
                return wp.n2p(f(zp.p2n(p)))
            img = sq.copy()
            fl = MT(tex, 36, col).move_to(DOWN * 3.35 + LEFT * 3.6)
            w0 = f(z0)
            # 方向导数：沿 Δx 与 iΔy 的差商
            dx = (f(z0 + 1e-6) - w0) / 1e-6
            dy = (f(z0 + 1e-6j) - w0) / 1e-6
            ax = Arrow(wp.n2p(w0), wp.n2p(w0 + 0.5 * dx), buff=0, color=C_RE, stroke_width=5, max_tip_length_to_length_ratio=0.2)
            ay = Arrow(wp.n2p(w0), wp.n2p(w0 + 0.5 * dy), buff=0, color=C_IM, stroke_width=5, max_tip_length_to_length_ratio=0.2)
            caps = VGroup(CT(cap1, 24), CT(cap2, 24, col)).arrange(DOWN, aligned_edge=LEFT, buff=0.1).move_to(DOWN * 3.3 + RIGHT * 2.9)
            self.play(FadeOut(prev), FadeIn(fl), run_time=0.5)
            self.play(Transform(img, img.copy().apply_function(F)), TransformFromCopy(ex, ax), TransformFromCopy(ey, ay), run_time=2)
            self.play(FadeIn(caps))
            self.wait(1.6)
            prev = VGroup(img, ax, ay, caps, fl)
        self.play(FadeOut(prev))
        summary = VGroup(
            CT(r"把 $f$ 在 $z_0$ 附近线性化：$\Delta w\approx\begin{pmatrix}u_x&u_y\\v_x&v_y\end{pmatrix}\begin{pmatrix}\Delta x\\\Delta y\end{pmatrix}$", 28),
            CT(r"复可导要求它等于乘以一个复数 $f'(z_0)=a+ib$：$\begin{pmatrix}a&-b\\b&a\end{pmatrix}$", 28),
            CT(r"$\Longrightarrow\ u_x=v_y,\quad u_y=-v_x$ \quad（柯西-黎曼方程）", 30, C_Z),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.18).to_edge(DOWN, buff=0.2)
        bg = BackgroundRectangle(summary, fill_opacity=0.9, buff=0.15)
        self.play(FadeIn(bg), FadeIn(summary))
        self.wait(3)


# --------------------------------------------------------------------------
# 2.2 调和共轭：u = c 与 v = c 两族曲线正交
# --------------------------------------------------------------------------
class C2_HarmonicOrthogonal(Scene):
    def construct(self):
        title = title_bar(r"解析函数 $f=u+iv$：等值线 $u=c_1$ 与 $v=c_2$ 处处正交")
        self.play(FadeIn(title))
        ax = plane((-2.5, 2.5, 1), (-2.5, 2.5, 1), unit=1.25).shift(LEFT * 3.2 + DOWN * 0.45)
        self.play(Create(ax))

        def fam(func, levels, color):
            g = VGroup()
            for c in levels:
                m = ImplicitFunction(lambda x, y, c=c: func(x, y) - c, x_range=[-2.5, 2.5], y_range=[-2.5, 2.5], color=color, stroke_width=2.5)
                m.scale(ax.get_x_unit_size(), about_point=ORIGIN).shift(ax.n2p(0))
                g.add(m)
            return g

        lab = CT(r"$f(z)=z^2$：$u=x^2-y^2,\ v=2xy$", 30, C_Z).to_edge(RIGHT, buff=0.4).shift(UP * 1.6)
        self.play(FadeIn(lab))
        U = fam(lambda x, y: x * x - y * y, [-3, -2, -1, -0.4, 0.4, 1, 2, 3], C_RE)
        V = fam(lambda x, y: 2 * x * y, [-3, -2, -1, -0.4, 0.4, 1, 2, 3], C_IM)
        lu = CT(r"红：$u=c_1$（等势线）", 28, C_RE).next_to(lab, DOWN, buff=0.35, aligned_edge=LEFT)
        lv = CT(r"蓝：$v=c_2$（电力线 / 流线）", 28, C_IM).next_to(lu, DOWN, buff=0.15, aligned_edge=LEFT)
        self.play(Create(U), FadeIn(lu), run_time=2)
        self.play(Create(V), FadeIn(lv), run_time=2)
        why = VGroup(
            CT(r"梯度 $\nabla u=(u_x,u_y)$，$\nabla v=(v_x,v_y)$", 26),
            CT(r"由 C-R：$\nabla u\cdot\nabla v=u_xv_x+u_yv_y$", 26),
            MT(r"=u_xv_x+(-v_x)(u_x)=0", 30, C_Z),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.15).next_to(lv, DOWN, buff=0.45, aligned_edge=LEFT)
        self.play(FadeIn(why))
        self.wait(2)
        self.play(FadeOut(VGroup(U, V, lab, why)))
        lab2 = CT(r"$f(z)=\dfrac1z$：$u=\dfrac{x}{x^2+y^2},\ v=\dfrac{-y}{x^2+y^2}$", 30, C_Z)
        lab2.scale_to_fit_width(min(lab2.width, 6.0)).move_to(lab, aligned_edge=LEFT).shift(UP * 0.1)
        U2 = fam(lambda x, y: x / (x * x + y * y + 1e-9), [-2, -1, -0.5, 0.5, 1, 2], C_RE)
        V2 = fam(lambda x, y: -y / (x * x + y * y + 1e-9), [-2, -1, -0.5, 0.5, 1, 2], C_IM)
        cap = CT(r"两族相切于原点的圆，同样正交（原点是奇点）", 26, GREY_B).next_to(lv, DOWN, buff=0.45, aligned_edge=LEFT)
        self.play(FadeIn(lab2))
        self.play(Create(U2), Create(V2), run_time=2.5)
        self.play(FadeIn(cap))
        self.wait(2.5)


# --------------------------------------------------------------------------
# 2.3 指数函数 w = e^z 的映射
# --------------------------------------------------------------------------
class C2_ExpMap(Scene):
    def construct(self):
        title = title_bar(r"$w=e^z=e^x e^{iy}$：实部管半径，虚部管角度")
        self.play(FadeIn(title))
        zp = plane((-2, 1.5, 1), (-0.5, 7, 1), unit=0.85).shift(LEFT * 4.3 + DOWN * 0.35)
        wp = plane((-4.5, 4.5, 1), (-4.5, 4.5, 1), unit=0.68).shift(RIGHT * 2.6 + DOWN * 0.35)
        zl = CT(r"$z$ 平面", 26, GREY_A).next_to(zp, RIGHT, buff=0.15).align_to(zp, UP)
        wl = CT(r"$w$ 平面", 26, GREY_A).next_to(wp, UP, buff=0.05)
        self.play(Create(zp), Create(wp), FadeIn(zl), FadeIn(wl))
        tp = MT(r"2\pi", 24, GREY_B).next_to(zp.n2p(2 * PI * 1j), LEFT, buff=0.05)
        pi_ = MT(r"\pi", 24, GREY_B).next_to(zp.n2p(PI * 1j), LEFT, buff=0.05)
        self.play(FadeIn(tp), FadeIn(pi_))

        def F(p):
            return wp.n2p(np.exp(zp.p2n(p)))

        # 竖线 x = c -> 圆周 |w| = e^c
        vl = VGroup(*[Line(zp.n2p(complex(c, 0)), zp.n2p(complex(c, 2 * PI)), color=C_RE, stroke_width=3) for c in (-1, -0.5, 0, 0.5, 1, 1.4)])
        for m in vl:
            m.insert_n_curves(80)
        cap = CT(r"竖线 $x=c\ \longmapsto$ 圆周 $|w|=e^c$", 28, C_RE).to_edge(DOWN, buff=0.25)
        self.play(Create(vl), FadeIn(cap))
        self.play(Transform(vl.copy(), vl.copy().apply_function(F)), run_time=2.5)
        self.wait(0.6)
        # 横线 y = c -> 射线 arg w = c
        hl = VGroup(*[Line(zp.n2p(complex(-2, c)), zp.n2p(complex(1.45, c)), color=C_IM, stroke_width=3) for c in np.linspace(0.3, 2 * PI - 0.3, 9)])
        for m in hl:
            m.insert_n_curves(40)
        cap2 = CT(r"横线 $y=c\ \longmapsto$ 射线 $\arg w=c$", 28, C_IM).to_edge(DOWN, buff=0.25)
        self.play(Create(hl), Transform(cap, cap2))
        self.play(Transform(hl.copy(), hl.copy().apply_function(F)), run_time=2.5)
        self.wait(0.8)
        self.play(*[FadeOut(m) for m in self.mobjects if m not in (zp, wp, zl, wl, title, tp, pi_)])

        # 带形 0 < Im z < π -> 上半平面
        band = Polygon(zp.n2p(-2 + 0.01j), zp.n2p(1.45 + 0.01j), zp.n2p(1.45 + PI * 1j), zp.n2p(-2 + PI * 1j), color=C_Z, fill_opacity=0.35, stroke_width=0)
        band.insert_n_curves(150)
        cap3 = CT(r"带形 $0<\operatorname{Im}z<\pi\ \longmapsto$ 上半平面 $\operatorname{Im}w>0$（第 6 章常用）", 28, C_Z).to_edge(DOWN, buff=0.25)
        self.play(FadeIn(band), FadeIn(cap3))
        self.play(Transform(band.copy(), band.copy().apply_function(F)), run_time=2.5)
        self.wait(0.8)
        per = CT(r"周期 $2\pi i$：带形 $\pi<\operatorname{Im}z<2\pi$ 盖住下半平面；再往上就重复覆盖", 28, C_W).to_edge(DOWN, buff=0.25)
        band2 = Polygon(zp.n2p(-2 + PI * 1j), zp.n2p(1.45 + PI * 1j), zp.n2p(1.45 + 2 * PI * 1j - 0.01j), zp.n2p(-2 + 2 * PI * 1j - 0.01j), color=C_W, fill_opacity=0.35, stroke_width=0)
        band2.insert_n_curves(150)
        self.play(FadeIn(band2), Transform(cap3, per))
        self.play(Transform(band2.copy(), band2.copy().apply_function(F)), run_time=2.5)
        self.wait(2)


# --------------------------------------------------------------------------
# 2.3 对数函数的多值性：Ln z 的虚部 = Arg z 是一个“螺旋楼梯”
# --------------------------------------------------------------------------
class C2_LogHelix(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=70 * DEGREES, theta=-60 * DEGREES, zoom=0.85)
        title = title_bar(r"$\operatorname{Ln}z=\ln|z|+i\operatorname{Arg}z$：绕原点一圈，虚部升高 $2\pi$")
        self.add_fixed_in_frame_mobjects(title)
        grid = NumberPlane(x_range=[-3, 3, 1], y_range=[-3, 3, 1], background_line_style={"stroke_color": BLUE_E, "stroke_width": 1, "stroke_opacity": 0.5})
        k = 0.35  # 高度缩放：z 方向 = k * Arg z
        helic = Surface(
            lambda r, t: np.array([r * np.cos(t), r * np.sin(t), k * t]),
            u_range=[0.25, 2.6], v_range=[-3 * PI, 3 * PI], resolution=(6, 72),
            fill_opacity=0.35, stroke_width=0.3, stroke_color=GREY_B, checkerboard_colors=[TEAL_D, TEAL_E],
        )
        self.play(FadeIn(title), Create(grid))
        self.play(Create(helic), run_time=2.5)
        cap = CT(r"曲面高度 $=\operatorname{Arg}z$（多值：每一层是一个分支）", 28, C_W).to_edge(DOWN, buff=0.7)
        self.add_fixed_in_frame_mobjects(cap)
        self.play(FadeIn(cap))
        th = ValueTracker(-PI + 0.05)
        R = 1.8
        pdot = always_redraw(lambda: Dot3D(np.array([R * np.cos(th.get_value()), R * np.sin(th.get_value()), 0]), color=C_Z, radius=0.07))
        hdot = always_redraw(lambda: Dot3D(np.array([R * np.cos(th.get_value()), R * np.sin(th.get_value()), k * th.get_value()]), color=RED, radius=0.08))
        stick = always_redraw(lambda: DashedLine(np.array([R * np.cos(th.get_value()), R * np.sin(th.get_value()), 0]), np.array([R * np.cos(th.get_value()), R * np.sin(th.get_value()), k * th.get_value()]), color=GREY_B))
        trail = always_redraw(lambda: ParametricFunction(lambda t: np.array([R * np.cos(t), R * np.sin(t), k * t]), t_range=[-PI + 0.05, max(th.get_value(), -PI + 0.06)], color=RED, stroke_width=4))
        self.add(trail, stick, pdot, hdot)
        self.begin_ambient_camera_rotation(rate=0.1)
        self.play(th.animate.set_value(3 * PI - 0.05), run_time=6, rate_func=linear)
        cap2 = CT(r"$z$ 回到原处，$\operatorname{Ln}z$ 的值却升了 $2\pi i$：多值函数", 28, C_HL).to_edge(DOWN, buff=0.7)
        self.add_fixed_in_frame_mobjects(cap2)
        self.play(FadeOut(cap), FadeIn(cap2))
        self.wait(1)
        # 支割线：负实轴；主值分支 -π < arg z ≤ π
        cut = Line3D(np.array([-2.8, 0, 0]), ORIGIN, color=RED, thickness=0.03)
        main = Surface(
            lambda r, t: np.array([r * np.cos(t), r * np.sin(t), k * t]),
            u_range=[0.25, 2.6], v_range=[-PI, PI], resolution=(6, 36),
            fill_opacity=0.7, stroke_width=0.3, checkerboard_colors=[YELLOW_D, YELLOW_E],
        )
        cap3 = CT(r"沿负实轴剪开（支割线），取 $-\pi<\arg z\le\pi$ 这一层：主值 $\ln z$", 28, C_Z).to_edge(DOWN, buff=0.7)
        self.add_fixed_in_frame_mobjects(cap3)
        self.play(FadeOut(cap2), FadeIn(cap3), Create(cut), FadeIn(main), run_time=2)
        self.wait(3)
        self.stop_ambient_camera_rotation()
