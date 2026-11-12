"""第一章 复数与复变函数 —— manim 动画源码.

渲染（720p30）：
    manim -qm ch1_scenes.py C1_VectorAndMultiply
或者直接运行仓库根目录的 tools/render_videos.sh
"""
import numpy as np

from common import *  # noqa: F401,F403


# --------------------------------------------------------------------------
# 1.1 复数 = 点 = 箭头；三角/指数形式；乘法 = 旋转 + 伸缩
# --------------------------------------------------------------------------
class C1_VectorAndMultiply(Scene):
    def construct(self):
        p = plane((-3, 3, 1), (-3, 3, 1), unit=1.2).shift(LEFT * 2.9 + DOWN * 0.45)
        labels = axis_labels(p)
        ticks = tick_numbers(p, xs=(1, 2), ys=(1, 2))
        title = T("复数 z = x + iy：平面上的一个点，也是一根箭头", 30).to_edge(UP, buff=0.25)
        self.play(Create(p), FadeIn(labels), FadeIn(ticks), Write(title), run_time=1.6)

        z = np.sqrt(3) + 1j
        zp = p.n2p(z)
        dot = Dot(zp, color=C_Z)
        zlab = T("z = √3 + i", 26, C_Z).next_to(dot, UR, buff=0.1)
        self.play(FadeIn(dot, scale=1.5), Write(zlab))

        # 实部、虚部投影
        vx = DashedLine(zp, p.n2p(z.real), color=C_RE)
        vy = DashedLine(zp, p.n2p(1j * z.imag), color=C_IM)
        rx = T("x = Re z = √3", 22, C_RE).next_to(p.n2p(z.real), DOWN, buff=0.35)
        ry = T("y = Im z = 1", 22, C_IM).next_to(p.n2p(1j * z.imag), LEFT, buff=0.45)
        self.play(Create(vx), Create(vy), FadeIn(rx), FadeIn(ry))

        arrow = Arrow(p.n2p(0), zp, buff=0, color=C_Z, stroke_width=5)
        ang = Angle(p.x_axis, Line(p.n2p(0), zp), radius=0.7, color=C_HL)
        th = T("θ = arg z = π/6", 22, C_HL).next_to(ang, RIGHT, buff=0.12).shift(UP * 0.05 + RIGHT * 0.25)
        rl = T("r = |z| = 2", 22, C_Z).rotate(np.pi / 6).move_to(p.n2p(z * 0.55) + UL * 0.32)
        self.play(GrowArrow(arrow), Create(ang), FadeIn(th), FadeIn(rl))

        right = panel(
            T("同一个复数，三种写法", 26, GREY_A),
            T("代数式  z = x + iy", 24),
            T("三角式  z = r(cosθ + i sinθ)", 24),
            T("指数式  z = r·e^(iθ)", 24, C_Z),
            T("x = r cosθ,  y = r sinθ", 22, GREY_B),
            T("r = √(x² + y²)", 22, GREY_B),
        ).to_edge(RIGHT, buff=0.35).shift(UP * 0.8)
        for m in right:
            self.play(FadeIn(m, shift=LEFT * 0.2), run_time=0.55)
        self.wait(1.2)

        # ---------------- 乘法 = 旋转 + 伸缩 ----------------
        self.play(FadeOut(VGroup(vx, vy, rx, ry, right, rl, th, ang)), run_time=0.8)
        title2 = T("乘以 w = 1 + i：长度 × √2，角度 + π/4", 30).to_edge(UP, buff=0.25)
        self.play(FadeTransform(title, title2)); title = title2

        wv = 1 + 1j
        warrow = Arrow(p.n2p(0), p.n2p(wv), buff=0, color=C_W, stroke_width=4)
        wlab = T("w = 1 + i = √2·e^(iπ/4)", 22, C_W).next_to(p.n2p(wv), UL, buff=0.1).shift(LEFT * 0.2)
        self.play(GrowArrow(warrow), FadeIn(wlab))

        rule = panel(
            T("z·w = r₁r₂ · e^(i(θ₁+θ₂))", 26, WHITE),
            T("模相乘：|zw| = |z|·|w|", 24, C_Z),
            T("辐角相加：Arg(zw) = Arg z + Arg w", 24, C_HL),
            T("(√3 + i)(1 + i) = (√3−1) + (√3+1)i", 22, GREY_B),
        ).to_edge(RIGHT, buff=0.3).shift(UP * 0.6)
        self.play(FadeIn(rule, shift=LEFT * 0.2))

        t = ValueTracker(0.0)
        r0, a0 = abs(z), np.angle(z)

        def cur():
            s = t.get_value()
            return r0 * (np.sqrt(2) ** s) * np.exp(1j * (a0 + s * np.pi / 4))

        moving = always_redraw(lambda: Arrow(p.n2p(0), p.n2p(cur()), buff=0, color=C_Z, stroke_width=5))
        trail = always_redraw(
            lambda: ParametricFunction(
                lambda u: p.n2p(r0 * (np.sqrt(2) ** u) * np.exp(1j * (a0 + u * np.pi / 4))),
                t_range=[0, max(t.get_value(), 1e-3)],
                color=GREY_B,
                stroke_width=2,
            )
        )
        ghost = arrow.copy().set_opacity(0.35)
        self.add(ghost)
        self.remove(arrow)
        self.add(trail, moving)
        self.play(FadeOut(zlab), FadeOut(dot))
        self.play(t.animate.set_value(1.0), run_time=2.5, rate_func=smooth)
        zw = z * wv
        zwlab = T("zw", 26, C_Z).next_to(p.n2p(zw), UP, buff=0.1)
        self.play(FadeIn(zwlab))
        self.wait(1.2)

        # ---------------- 乘以 i = 逆时针转 90° ----------------
        self.play(FadeOut(VGroup(rule, warrow, wlab, zwlab, ghost)), run_time=0.7)
        self.remove(moving, trail)
        title3 = T("特别地：乘以 i 就是逆时针旋转 90°", 30).to_edge(UP, buff=0.25)
        self.play(FadeTransform(title, title3)); title = title3
        circle = Circle(radius=p.get_x_unit_size() * 2, color=GREY_C, stroke_width=1.5).move_to(p.n2p(0))
        self.play(Create(circle))
        phi = ValueTracker(a0)
        rot = always_redraw(lambda: Arrow(p.n2p(0), p.n2p(2 * np.exp(1j * phi.get_value())), buff=0, color=C_Z, stroke_width=5))
        self.add(rot)
        names = ["z", "iz", "i²z = −z", "i³z = −iz", "i⁴z = z"]
        lab = T(names[0], 24, C_Z).next_to(p.n2p(2 * np.exp(1j * a0)), UR, buff=0.05)
        self.play(FadeIn(lab))
        for k in range(1, 5):
            self.play(phi.animate.set_value(a0 + k * np.pi / 2), run_time=1.0)
            new = T(names[k], 24, C_Z).next_to(p.n2p(2.25 * np.exp(1j * (a0 + k * np.pi / 2))), ORIGIN)
            self.play(Transform(lab, new), run_time=0.4)
        note = boxed(panel(T("i² = −1 的几何意义：", 24, GREY_A), T("转两次 90° = 掉头 180°", 26, C_HL))).to_edge(RIGHT, buff=0.3)
        self.play(FadeIn(note))
        self.wait(2)


# --------------------------------------------------------------------------
# 1.2 乘幂与方根：棣莫弗公式、单位根、方根的多值性
# --------------------------------------------------------------------------
class C1_PowersAndRoots(Scene):
    def construct(self):
        p = plane((-2.5, 2.5, 1), (-2.5, 2.5, 1), unit=1.45).shift(LEFT * 3.0 + DOWN * 0.25)
        title = T("乘幂 zⁿ = rⁿ e^(inθ)：模乘 n 次方，角度乘 n", 30).to_edge(UP, buff=0.25)
        ticks = tick_numbers(p, xs=(1, 2), ys=(1, 2))
        unit = Circle(radius=p.get_x_unit_size(), color=GREY_C, stroke_width=1.5).move_to(p.n2p(0))
        self.play(Create(p), FadeIn(ticks), Create(unit), Write(title), run_time=1.5)

        z = 1.12 * np.exp(1j * np.pi / 7)
        pts = [z ** k for k in range(0, 8)]
        dots = VGroup(*[Dot(p.n2p(w), radius=0.06, color=interpolate_color(C_Z, RED, k / 7)) for k, w in enumerate(pts)])
        names = ["1", "z", "z²", "z³", "z⁴", "z⁵", "z⁶", "z⁷"]
        info = panel(
            T("z = 1.12·e^(iπ/7)", 24, C_Z),
            T("每乘一次 z：", 24),
            T("  · 箭头转 π/7", 22, C_HL),
            T("  · 长度 ×1.12", 22, C_HL),
            T("⇒ 点列沿对数螺线向外", 22, GREY_B),
        ).to_edge(RIGHT, buff=0.4).shift(UP * 0.6)
        self.play(FadeIn(info))
        self.add(dots[0])
        prev = dots[0]
        for k in range(1, 8):
            seg = Line(p.n2p(pts[k - 1]), p.n2p(pts[k]), color=GREY_B, stroke_width=2)
            lab = T(names[k], 20, dots[k].get_color()).next_to(dots[k], (dots[k].get_center() - p.n2p(0)) / max(np.linalg.norm(dots[k].get_center() - p.n2p(0)), 1e-6), buff=0.08)
            self.play(Create(seg), FadeIn(dots[k]), FadeIn(lab), run_time=0.45)
        self.wait(0.8)
        self.play(*[FadeOut(m) for m in self.mobjects if m not in (p, ticks, unit, title)], run_time=0.7)

        # ---------------- 单位根 ----------------
        title2 = T("zⁿ = 1 的 n 个根：均匀分布在单位圆上", 30).to_edge(UP, buff=0.25)
        self.play(FadeTransform(title, title2)); title = title2
        formula = panel(
            T("ⁿ√1 = e^(i·2kπ/n)", 28, C_Z),
            T("k = 0, 1, …, n−1", 24, GREY_B),
            T("相邻两根夹角 2π/n", 24, C_HL),
        ).to_edge(RIGHT, buff=0.5).shift(UP * 0.8)
        self.play(FadeIn(formula))

        def roots_group(n, r=1.0, phase=0.0, color=C_Z):
            ws = [r * np.exp(1j * (phase + 2 * np.pi * k / n)) for k in range(n)]
            poly = Polygon(*[p.n2p(w) for w in ws], color=color, stroke_width=2.5).set_fill(color, 0.12)
            ds = VGroup(*[Dot(p.n2p(w), color=color, radius=0.07) for w in ws])
            return VGroup(poly, ds)

        nlab = T("n = 3", 30, C_HL).next_to(formula, DOWN, buff=0.6, aligned_edge=LEFT)
        g = roots_group(3)
        self.play(Create(g), FadeIn(nlab))
        for n in (4, 5, 6, 8):
            self.play(Transform(g, roots_group(n)), Transform(nlab, T(f"n = {n}", 30, C_HL).move_to(nlab, aligned_edge=LEFT)), run_time=0.9)
            self.wait(0.35)
        self.wait(0.6)
        self.play(FadeOut(g), FadeOut(nlab), FadeOut(formula))

        # ---------------- 方根的多值性：z 绕一圈，三个立方根轮换 ----------------
        title3 = T("多值性：z 绕原点一圈，∛z 的三个值依次轮换", 30).to_edge(UP, buff=0.25)
        self.play(FadeTransform(title, title3)); title = title3
        big = Circle(radius=p.get_x_unit_size() * 2, color=GREY_D, stroke_width=1.5).move_to(p.n2p(0))
        phi = ValueTracker(0.0)
        R = 2.0  # |z| = 2  →  |∛z| = 2^(1/3)
        rr = R ** (1 / 3)
        zdot = always_redraw(lambda: Dot(p.n2p(R * np.exp(1j * phi.get_value())), color=C_W, radius=0.09))
        zarc = always_redraw(lambda: Arc(radius=p.get_x_unit_size() * R, start_angle=0, angle=max(phi.get_value(), 1e-3), arc_center=p.n2p(0), color=C_W, stroke_width=3))
        cols = [C_Z, RED_C, GREEN_C]
        rdots = [
            always_redraw(lambda k=k: Dot(p.n2p(rr * np.exp(1j * (phi.get_value() + 2 * np.pi * k) / 3)), color=cols[k], radius=0.09))
            for k in range(3)
        ]
        rtrails = [
            always_redraw(lambda k=k: Arc(radius=p.get_x_unit_size() * rr, start_angle=2 * np.pi * k / 3, angle=max(phi.get_value() / 3, 1e-3), arc_center=p.n2p(0), color=cols[k], stroke_width=3))
            for k in range(3)
        ]
        leg = panel(
            T("z = 2e^(iφ)   (青色)", 24, C_W),
            T("w_k = ∛2 · e^(i(φ+2kπ)/3)", 24, C_Z),
            T("φ 走 2π，每个 w_k 只走 2π/3", 22, GREY_B),
            T("→ w₀ 走到 w₁ 的位置", 22, GREY_B),
            T("这就是第二章“支割线”的伏笔", 22, C_HL),
        ).to_edge(RIGHT, buff=0.35).shift(UP * 0.5)
        self.play(Create(big), FadeIn(leg))
        self.add(zarc, zdot, *rtrails, *rdots)
        self.play(phi.animate.set_value(2 * np.pi), run_time=5, rate_func=linear)
        self.wait(0.6)
        self.play(phi.animate.set_value(4 * np.pi), run_time=4, rate_func=linear)
        self.wait(1.5)


# --------------------------------------------------------------------------
# 1.5 复变函数 = 映射：w = z² 把网格变成抛物线族、把四分之一圆盘变成上半圆盘
# --------------------------------------------------------------------------
class C1_MapSquare(Scene):
    def construct(self):
        zp = plane((-1.5, 1.5, 0.5), (-1.5, 1.5, 0.5), unit=1.45).shift(LEFT * 3.5 + DOWN * 0.45)
        wp = plane((-3, 3, 1), (-3, 3, 1), unit=0.75).shift(RIGHT * 3.5 + DOWN * 0.45)
        zt = T("z 平面", 26, GREY_A).next_to(zp, UP, buff=0.12)
        wt = T("w 平面", 26, GREY_A).next_to(wp, UP, buff=0.12)
        arrow = Arrow(LEFT * 0.9, RIGHT * 0.9, color=WHITE).shift(DOWN * 0.45)
        fl = T("w = z²", 30, C_Z).next_to(arrow, UP, buff=0.1)
        title = T("复变函数 = 把一个平面“搬”到另一个平面的映射", 30).to_edge(UP, buff=0.2)
        self.play(Write(title), Create(zp), Create(wp), FadeIn(zt), FadeIn(wt), GrowArrow(arrow), Write(fl), run_time=1.8)

        def f(point):
            z = zp.p2n(point)
            return wp.n2p(z * z)

        # 竖线 x = c （红）与 横线 y = c （蓝）
        cs = [0.25, 0.5, 0.75, 1.0]
        vlines = VGroup(*[Line(zp.n2p(complex(c, -1.2)), zp.n2p(complex(c, 1.2)), color=C_RE, stroke_width=3) for c in cs])
        hlines = VGroup(*[Line(zp.n2p(complex(-1.2, c)), zp.n2p(complex(1.2, c)), color=C_IM, stroke_width=3) for c in cs])
        for m in (*vlines, *hlines):
            m.insert_n_curves(60)
        self.play(Create(vlines), Create(hlines), run_time=1.2)
        uv = T("u = x² − y²,   v = 2xy", 24, GREY_A).to_edge(DOWN, buff=0.25)
        self.play(FadeIn(uv))
        img_v = vlines.copy()
        img_h = hlines.copy()
        self.play(ApplyPointwiseFunction(f, img_v), ApplyPointwiseFunction(f, img_h), run_time=3)
        cap = T("直线 x=c → 向左开口的抛物线；y=c → 向右开口的抛物线", 22, GREY_B).next_to(uv, UP, buff=0.12)
        self.play(FadeIn(cap))
        self.wait(1.8)
        self.play(FadeOut(VGroup(vlines, hlines, img_v, img_h, cap)))

        # 极坐标网格：射线角度翻倍，圆周半径平方
        cap2 = T("用极坐标看最清楚：w = r²·e^(i2θ) —— 半径平方，角度翻倍", 22, GREY_B).next_to(uv, UP, buff=0.12)
        self.play(Transform(uv, T("z = r e^(iθ)  →  w = r² e^(i2θ)", 24, GREY_A).move_to(uv)), FadeIn(cap2))
        rays = VGroup(*[Line(zp.n2p(0), zp.n2p(1.0 * np.exp(1j * a)), color=C_RE, stroke_width=3) for a in np.linspace(0, np.pi / 2, 5)])
        arcs = VGroup(*[ParametricFunction(lambda t, r=r: zp.n2p(r * np.exp(1j * t)), t_range=[0, np.pi / 2], color=C_IM, stroke_width=3) for r in (0.4, 0.7, 1.0)])
        sector = Sector(radius=zp.get_x_unit_size() * 1.0, start_angle=0, angle=PI / 2, arc_center=zp.n2p(0), color=C_Z, fill_opacity=0.35, stroke_width=0)
        sector.insert_n_curves(80)
        for m in (*rays,):
            m.insert_n_curves(30)
        dlab = T("D: 0<|z|<1, 0<arg z<π/2", 22, C_Z).next_to(zp.n2p(0.5 + 1.05j), UP, buff=0.05)
        self.play(FadeIn(sector), Create(rays), Create(arcs), FadeIn(dlab), run_time=1.4)
        img = VGroup(sector.copy(), rays.copy(), arcs.copy())
        self.play(ApplyPointwiseFunction(f, img), run_time=3)
        glab = T("G: |w|<1, Im w>0（上半圆盘）", 22, C_Z).next_to(wp.n2p(0 + 1.05j), UP, buff=0.05)
        self.play(FadeIn(glab))
        # 一个动点：角度翻倍
        th = ValueTracker(0.05)
        zd = always_redraw(lambda: Dot(zp.n2p(0.85 * np.exp(1j * th.get_value())), color=WHITE))
        wd = always_redraw(lambda: Dot(wp.n2p(0.85 ** 2 * np.exp(2j * th.get_value())), color=WHITE))
        self.add(zd, wd)
        self.play(th.animate.set_value(np.pi / 2 - 0.05), run_time=2.5, rate_func=there_and_back_with_pause)
        self.wait(1.5)


# --------------------------------------------------------------------------
# 1.5 极限：z → 0 的方式是任意的。f(z) = z / z̄ 沿不同射线趋于不同的值
# --------------------------------------------------------------------------
class C1_LimitPaths(Scene):
    def construct(self):
        zp = plane((-1.5, 1.5, 0.5), (-1.5, 1.5, 0.5), unit=1.45).shift(LEFT * 3.5 + DOWN * 0.4)
        wp = plane((-1.5, 1.5, 0.5), (-1.5, 1.5, 0.5), unit=1.45).shift(RIGHT * 3.5 + DOWN * 0.4)
        zt = T("z 平面", 26, GREY_A).next_to(zp, UP, buff=0.12)
        wt = T("w = f(z) 平面", 26, GREY_A).next_to(wp, UP, buff=0.12)
        title = T("极限要求“从任何方向逼近都得到同一个值”", 30).to_edge(UP, buff=0.2)
        uc = Circle(radius=zp.get_x_unit_size(), color=GREY_C, stroke_width=1.5).move_to(wp.n2p(0))
        self.play(Write(title), Create(zp), Create(wp), FadeIn(zt), FadeIn(wt), Create(uc), run_time=1.6)
        fl = T("f(z) = z / z̄ ，沿射线 z = r·e^(iθ) 有 f = e^(i2θ)（与 r 无关）", 24, C_Z).next_to(title, DOWN, buff=0.15)
        self.play(FadeIn(fl))

        angles = [0, np.pi / 6, np.pi / 4, np.pi / 2, 3 * np.pi / 4]
        cols = [RED_C, ORANGE, YELLOW, GREEN_C, BLUE_C]
        s = ValueTracker(1.2)
        mobs = []
        for a, c in zip(angles, cols):
            ray = DashedLine(zp.n2p(1.3 * np.exp(1j * a)), zp.n2p(0), color=c, stroke_width=2)
            zd = always_redraw(lambda a=a, c=c: Dot(zp.n2p(s.get_value() * np.exp(1j * a)), color=c))
            wd = always_redraw(lambda a=a, c=c: Dot(wp.n2p(np.exp(2j * a)), color=c, radius=0.1))
            mobs += [ray, zd, wd]
        self.play(*[Create(m) for m in mobs[0::3]], run_time=1)
        self.add(*mobs[1::3], *mobs[2::3])
        self.play(s.animate.set_value(0.05), run_time=3.5, rate_func=smooth)
        txt = VGroup(
            T("沿 x 轴（θ=0）：f → 1", 22, RED_C),
            T("沿 y 轴（θ=π/2）：f → −1", 22, GREEN_C),
            T("沿 θ=π/4：f → i", 22, YELLOW),
            T("值依赖方向 ⇒ lim f(z) 不存在", 24, WHITE),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.12).to_edge(DOWN, buff=0.15).shift(RIGHT * 0.2)
        bg = BackgroundRectangle(txt, fill_opacity=0.85, buff=0.15)
        self.play(FadeIn(bg), FadeIn(txt))
        self.wait(2)

        # 对比：|z|² 连续
        self.play(FadeOut(VGroup(bg, txt)), *[FadeOut(m) for m in mobs], FadeOut(fl))
        title2 = T("对比：f(z) = |z|² —— 不管从哪个方向，都趋于 0", 30).to_edge(UP, buff=0.2)
        self.play(FadeTransform(title, title2)); title = title2
        s.set_value(1.2)
        mobs2 = []
        for a, c in zip(angles, cols):
            ray = DashedLine(zp.n2p(1.3 * np.exp(1j * a)), zp.n2p(0), color=c, stroke_width=2)
            zd = always_redraw(lambda a=a, c=c: Dot(zp.n2p(s.get_value() * np.exp(1j * a)), color=c))
            wd = always_redraw(lambda a=a, c=c: Dot(wp.n2p(s.get_value() ** 2 + 0j) + UP * 0.1 * (angles.index(a) - 2), color=c, radius=0.08))
            mobs2 += [ray, zd, wd]
        self.play(*[Create(m) for m in mobs2[0::3]], run_time=0.8)
        self.add(*mobs2[1::3], *mobs2[2::3])
        self.play(s.animate.set_value(0.0), run_time=3, rate_func=smooth)
        end = T("|z|² 在 0 点连续：lim f(z) = 0 = f(0)", 24).to_edge(DOWN, buff=0.3)
        self.play(FadeIn(end))
        self.wait(2)


# --------------------------------------------------------------------------
# 1.4 复球面：球极投影，无穷远点 = 北极 N
# --------------------------------------------------------------------------
def _stereo(x, y, R=1.0):
    """平面点 (x, y, 0) 与北极 N=(0,0,2R) 连线交球面（与平面切于原点）的点."""
    t = 4 * R * R / (x * x + y * y + 4 * R * R)
    return np.array([t * x, t * y, 2 * R * (1 - t)])


class C1_RiemannSphere(ThreeDScene):
    def construct(self):
        R = 1.0
        self.set_camera_orientation(phi=68 * DEGREES, theta=-50 * DEGREES, zoom=0.95)
        grid = NumberPlane(
            x_range=[-6, 6, 1], y_range=[-6, 6, 1],
            background_line_style={"stroke_color": BLUE_E, "stroke_width": 1, "stroke_opacity": 0.5},
        )
        sphere = Surface(
            lambda u, v: np.array([R * np.cos(u) * np.sin(v), R * np.sin(u) * np.sin(v), R + R * np.cos(v)]),
            u_range=[0, TAU], v_range=[0, PI], resolution=(24, 12),
            fill_opacity=0.18, stroke_width=0.4, stroke_color=GREY_B, checkerboard_colors=[BLUE_D, BLUE_E],
        )
        N = Dot3D(np.array([0, 0, 2 * R]), color=RED, radius=0.07)
        S = Dot3D(ORIGIN, color=GREEN, radius=0.06)
        title = T("复球面：用“过北极 N 的直线”把平面点放到球面上", 28).to_edge(UP, buff=0.2)
        self.add_fixed_in_frame_mobjects(title)
        self.play(Create(grid), FadeIn(title), run_time=1.2)
        self.play(Create(sphere), FadeIn(N), FadeIn(S), run_time=1.5)
        nl = T("N 北极", 22, RED).to_corner(UL, buff=0.4).shift(DOWN * 0.7)
        sl = T("S 南极 = 原点 0", 22, GREEN).next_to(nl, DOWN, aligned_edge=LEFT)
        self.add_fixed_in_frame_mobjects(nl, sl)
        self.play(FadeIn(nl), FadeIn(sl))

        # 一个平面点沿螺线向外跑
        s = ValueTracker(0.0)

        def P():
            u = s.get_value()
            r = 0.6 * np.exp(u)
            return np.array([r * np.cos(2.2 * u), r * np.sin(2.2 * u), 0.0])

        pdot = always_redraw(lambda: Dot3D(P(), color=C_Z, radius=0.07))
        qdot = always_redraw(lambda: Dot3D(_stereo(P()[0], P()[1], R), color=C_Z, radius=0.06))
        ray = always_redraw(lambda: Line3D(np.array([0, 0, 2 * R]), P(), color=C_Z, thickness=0.01))
        self.add(ray, pdot, qdot)
        cap1 = T("平面上的点 z  ↔  球面上除 N 以外的点（一一对应）", 22, C_Z).to_edge(DOWN, buff=0.6)
        self.add_fixed_in_frame_mobjects(cap1)
        self.play(FadeIn(cap1))
        self.play(s.animate.set_value(1.6), run_time=4, rate_func=smooth)
        self.begin_ambient_camera_rotation(rate=0.12)
        cap2 = T("|z| 越大，像点越靠近 N；把 N 本身对应到“无穷远点 ∞”", 22, C_HL).next_to(cap1, DOWN, buff=0.12)
        self.add_fixed_in_frame_mobjects(cap2)
        self.play(FadeIn(cap2))
        self.play(s.animate.set_value(2.25), run_time=3, rate_func=smooth)
        self.wait(0.5)

        # ∞ 的邻域：|z| > M  ↔  N 附近的小球冠
        self.play(FadeOut(ray), FadeOut(pdot), FadeOut(qdot), FadeOut(cap1), FadeOut(cap2))
        M = 3.0
        circ = ParametricFunction(lambda t: np.array([M * np.cos(t), M * np.sin(t), 0]), t_range=[0, TAU], color=C_HL, stroke_width=3)
        lat = ParametricFunction(lambda t: _stereo(M * np.cos(t), M * np.sin(t), R), t_range=[0, TAU], color=C_HL, stroke_width=4)
        cap3 = T("圆 |z| = M 外部（∞ 的邻域）  ↔  N 附近的一顶小“帽子”", 22, C_HL).to_edge(DOWN, buff=0.6)
        self.add_fixed_in_frame_mobjects(cap3)
        self.play(Create(circ), FadeIn(cap3))
        self.play(TransformFromCopy(circ, lat), run_time=2)
        cap4 = T("扩充复平面 ℂ ∪ {∞} 就是一个完整的球面", 24, WHITE).next_to(cap3, DOWN, buff=0.12)
        self.add_fixed_in_frame_mobjects(cap4)
        self.play(FadeIn(cap4))
        self.wait(3)
        self.stop_ambient_camera_rotation()
