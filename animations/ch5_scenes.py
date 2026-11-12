"""第五章 留数及其应用 —— manim 动画源码（LaTeX 排版）."""
import colorsys

import numpy as np

from common import *  # noqa: F401,F403


def fitw(m, w):
    if m.width > w:
        m.scale_to_fit_width(w)
    return m


def domain_color(f, half=2.0, n=360, center=0j):
    """域着色：色相 = arg f，亮度随 |f| 增大而变亮，并叠加 |f| 的对数等高环。"""
    xs = np.linspace(-half, half, n)
    X, Y = np.meshgrid(xs, xs[::-1])
    Z = X + 1j * Y + center
    with np.errstate(all="ignore"):
        W = f(Z)
    W = np.where(np.isfinite(W), W, 1e8)
    H = (np.angle(W) / (2 * np.pi)) % 1.0
    M = np.abs(W)
    L = 2 / np.pi * np.arctan(np.sqrt(M))          # 0 -> 黑, ∞ -> 白
    ring = (np.log2(M + 1e-30) % 1.0)
    V = 0.55 + 0.45 * ring
    rgb = np.zeros(Z.shape + (3,))
    hsv = np.vectorize(colorsys.hls_to_rgb)
    r, g, b = hsv(H, 0.08 + 0.84 * L * (0.85 + 0.15 * V), 0.95)
    rgb[..., 0], rgb[..., 1], rgb[..., 2] = r, g, b
    return (np.clip(rgb, 0, 1) * 255).astype(np.uint8)


def dc_image(f, size=3.6, **kw):
    img = ImageMobject(domain_color(f, **kw))
    img.set_resampling_algorithm(RESAMPLING_ALGORITHMS["bilinear"])
    img.height = size
    return img


# --------------------------------------------------------------------------
# 5.1 域着色看三类孤立奇点
# --------------------------------------------------------------------------
class C5_DomainColoring(Scene):
    def construct(self):
        title = title_bar(r"域着色：颜色 $=\arg f(z)$，亮度 $=|f(z)|$（黑 $\to0$，白 $\to\infty$）")
        self.play(FadeIn(title))
        S = 3.6
        cases = [
            (lambda z: np.sin(z) / z, r"\dfrac{\sin z}{z}", r"可去奇点：颜色平滑，$\lim f=1$"),
            (lambda z: 1 / z ** 2, r"\dfrac{1}{z^{2}}", r"二阶极点：变白（$\to\infty$），颜色绕 $2$ 圈"),
            (lambda z: np.exp(1 / z), r"e^{1/z}", r"本性奇点：任意小邻域内取遍几乎所有值"),
        ]
        imgs = Group()
        labs = VGroup()
        caps = VGroup()
        for k, (f, tex, cap) in enumerate(cases):
            im = dc_image(f, size=S, half=1.5)
            im.move_to(RIGHT * (k - 1) * 4.3 + UP * 0.1)
            imgs.add(im)
            labs.add(MT(tex, 34).next_to(im, UP, buff=0.15))
            caps.add(fitw(CT(cap, 22), 4.1).next_to(im, DOWN, buff=0.2))
        for im, lb, cp in zip(imgs, labs, caps):
            dot = Dot(im.get_center(), radius=0.04, color=WHITE)
            self.play(FadeIn(im), FadeIn(lb), FadeIn(dot), run_time=0.9)
            self.play(FadeIn(cp), run_time=0.6)
            self.wait(1.0)
        self.wait(1.0)

        # 放大本性奇点
        self.play(FadeOut(imgs[:2]), FadeOut(labs[:2]), FadeOut(caps[:2]), *[FadeOut(m) for m in self.mobjects if isinstance(m, Dot)])
        big = Group(imgs[2], labs[2], caps[2])
        self.play(big.animate.move_to(LEFT * 2.6 + DOWN * 0.2).scale(1.35))
        cur = imgs[2]
        zl = None
        for half in (0.5, 0.15, 0.05):
            nim = dc_image(lambda z: np.exp(1 / z), size=cur.height, half=half).move_to(cur)
            nl = CT(rf"放大：窗口 $|x|,|y|\le{half}$", 24, C_HL).next_to(cur, RIGHT, buff=0.4).shift(UP * 1.2)
            self.play(FadeIn(nim), FadeOut(cur), FadeTransform(zl, nl) if zl else FadeIn(nl), run_time=1.2)
            cur, zl = nim, nl
            self.wait(0.6)
        note = VGroup(
            CT(r"越放大越“乱”：没有极限，也不趋于 $\infty$", 26),
            CT(r"洛朗级数 $e^{1/z}=\sum_{n\ge0}\frac{1}{n!\,z^{n}}$ 有无穷多负幂项", 26, C_Z),
            CT(r"（皮卡定理：除至多一个值外，取遍所有复数无穷多次）", 22, GREY_B),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.25)
        fitw(note, 6.0).next_to(zl, DOWN, buff=0.5, aligned_edge=LEFT)
        self.play(LaggedStart(*[FadeIn(n) for n in note], lag_ratio=0.5), run_time=2)
        self.wait(2.5)


# --------------------------------------------------------------------------
# 5.2 留数定理 = 复合闭路 + 每个小圆只留下 2πi c_{-1}
# --------------------------------------------------------------------------
class C5_ResidueTheorem(Scene):
    def construct(self):
        title = title_bar(r"留数定理：$\displaystyle\oint_C f\,\mathrm dz=2\pi i\sum_k\operatorname{Res}[f,z_k]$")
        self.play(FadeIn(title))
        zp = plane((-3, 3, 1), (-2.2, 2.2, 1), unit=1.1).shift(LEFT * 2.9 + DOWN * 0.55)
        self.play(Create(zp))
        pts = [-1.3 + 0.6j, 0.4 - 0.7j, 1.5 + 0.8j]
        res_val = [r"\operatorname{Res}[f,z_1]", r"\operatorname{Res}[f,z_2]", r"\operatorname{Res}[f,z_3]"]
        dots = VGroup(*[Dot(zp.n2p(p), color=RED) for p in pts])
        dl = VGroup(*[MT(f"z_{k+1}", 26, RED).next_to(d, UR, buff=0.04) for k, d in enumerate(dots)])
        self.play(FadeIn(dots), FadeIn(dl))
        s = ValueTracker(0.0)
        r = 0.32

        def big(t):
            return (2.5 + 0.25 * np.cos(3 * t)) * np.cos(t) + 1j * (1.9 + 0.1 * np.sin(2 * t)) * np.sin(t)

        C = ParametricFunction(lambda t: zp.n2p(big(t)), t_range=[0, 2 * PI], color=C_Z, stroke_width=4)
        Cl = MT("C", 30, C_Z).next_to(zp.n2p(big(0.3)), RIGHT, buff=0.08)
        self.play(Create(C), FadeIn(Cl))
        X = 0.75

        def put(m, y):
            fitw(m, 5.6)
            m.set_x(X + m.width / 2 + 0.2).set_y(y)
            return m

        t1 = put(CT(r"$f$ 在 $C$ 内除 $z_1,z_2,z_3$ 外解析", 26), 1.5)
        self.play(FadeIn(t1))
        circles = VGroup(*[Circle(radius=r * zp.get_x_unit_size(), color=C_W).move_to(zp.n2p(p)) for p in pts])
        self.play(*[TransformFromCopy(C, c) for c in circles], run_time=2)
        t2 = put(CT(r"复合闭路定理：$\displaystyle\oint_C=\sum_{k}\oint_{|z-z_k|=r}$", 26, C_W), 0.7)
        self.play(FadeIn(t2))
        # 小圆放大展示洛朗级数逐项积分
        t3 = put(CT(r"在 $z_k$ 附近：$f=\cdots+\dfrac{c_{-2}}{(z-z_k)^2}+\dfrac{c_{-1}}{z-z_k}+c_0+c_1(z-z_k)+\cdots$", 24), -0.15)
        t4 = put(CT(r"逐项积分（例题 3.2）：只有 $\dfrac{c_{-1}}{z-z_k}$ 一项给出 $2\pi i\,c_{-1}$", 24, C_HL), -0.95)
        self.play(FadeIn(t3))
        self.play(*[Indicate(c, color=C_HL) for c in circles], FadeIn(t4))
        labs = VGroup(*[MT(rf"2\pi i\,{res_val[k]}", 20, C_HL).next_to(circles[k], DOWN, buff=0.06) for k in range(3)])
        self.play(FadeIn(labs))
        res = put(CT(r"$\displaystyle\oint_Cf\,\mathrm dz=2\pi i\big(\operatorname{Res}[f,z_1]+\operatorname{Res}[f,z_2]+\operatorname{Res}[f,z_3]\big)$", 26, C_Z), -1.9)
        box = SurroundingRectangle(res, color=C_Z, buff=0.12)
        self.play(FadeIn(res), Create(box))
        self.wait(1.5)
        warn = put(CT(r"只算 $C$ \textbf{内部}的奇点；外部的奇点不参与", 24, RED_B), -2.75)
        far = Dot(zp.n2p(2.8 + 1.9j), color=GREY_B)
        fl = MT("z_4", 24, GREY_B).next_to(far, LEFT, buff=0.05)
        self.play(FadeIn(far), FadeIn(fl), FadeIn(warn))
        self.wait(3)


# --------------------------------------------------------------------------
# 5.3 上半圆周补线法：∫ dx/(1+x^2) = π
# --------------------------------------------------------------------------
class C5_SemicircleContour(Scene):
    def construct(self):
        title = title_bar(r"补上半圆周：$\displaystyle\int_{-\infty}^{+\infty}\frac{\mathrm dx}{1+x^2}=2\pi i\operatorname{Res}\Big[\frac{1}{1+z^2},i\Big]=\pi$")
        self.play(FadeIn(title))
        zp = plane((-5, 5, 1), (-1, 5, 1), unit=0.72).shift(LEFT * 2.6 + DOWN * 0.9)
        self.play(Create(zp))
        pole = Dot(zp.n2p(1j), color=RED)
        pole2 = Dot(zp.n2p(-1j), color=GREY_B).set_opacity(0.6)
        pl = MT("i", 26, RED).next_to(pole, RIGHT, buff=0.06)
        self.play(FadeIn(pole), FadeIn(pl), FadeIn(pole2))
        R = ValueTracker(1.6)
        f = lambda z: 1 / (1 + z * z)

        def arc_int(Rv, N=800):
            t = np.linspace(0, PI, N + 1)
            z = Rv * np.exp(1j * t)
            return np.sum(f((z[1:] + z[:-1]) / 2) * np.diff(z))

        def line_int(Rv):
            return 2 * np.arctan(Rv)

        seg = always_redraw(lambda: Line(zp.n2p(-R.get_value()), zp.n2p(R.get_value()), color=C_W, stroke_width=5))
        arc = always_redraw(lambda: ParametricFunction(lambda t: zp.n2p(R.get_value() * np.exp(1j * t)), t_range=[0, PI], color=C_Z, stroke_width=5))
        self.play(Create(seg), Create(arc))
        X = 2.2

        def panel_():
            Rv = R.get_value()
            ai = arc_int(Rv)
            g = VGroup(
                MT(rf"R={Rv:.1f}", 30),
                MT(rf"\int_{{-R}}^{{R}}\frac{{\mathrm dx}}{{1+x^2}}={line_int(Rv):.4f}", 28, C_W),
                MT(rf"\int_{{C_R}}\frac{{\mathrm dz}}{{1+z^2}}={ai.real:+.4f}{ai.imag:+.4f}i", 28, C_Z),
                CT(rf"两者之和 $={line_int(Rv)+ai.real:.4f}\ (=\pi)$", 28, C_HL),
            ).arrange(DOWN, aligned_edge=LEFT, buff=0.25)
            g.set_x(X + g.width / 2).set_y(1.0)
            return g

        pan = always_redraw(panel_)
        self.play(FadeIn(pan))
        cap = fitw(CT(r"闭路 $=[-R,R]+C_R$，内含奇点 $i$：$\oint=2\pi i\cdot\dfrac{1}{2i}=\pi$，与 $R$ 无关", 24), 13).to_edge(DOWN, buff=0.2)
        self.play(FadeIn(cap))
        self.play(R.animate.set_value(4.8), run_time=6, rate_func=linear)
        self.wait(0.6)
        est = fitw(CT(r"大圆弧估计：$\left|\displaystyle\int_{C_R}\right|\le\dfrac{\pi R}{R^2-1}\to0$，所以实积分 $=\pi$", 26, C_Z), 13).to_edge(DOWN, buff=0.2)
        self.play(Transform(cap, est))
        self.wait(3)


# --------------------------------------------------------------------------
# 5.3 z=e^{iθ}：把 [0,2π] 卷成单位圆
# --------------------------------------------------------------------------
class C5_UnitCircleSub(Scene):
    def construct(self):
        title = title_bar(r"$\displaystyle\int_0^{2\pi}R(\cos\theta,\sin\theta)\,\mathrm d\theta$：令 $z=e^{i\theta}$，区间卷成单位圆")
        self.play(FadeIn(title))
        zp = plane((-1.6, 1.6, 0.5), (-1.6, 1.6, 0.5), unit=1.45).shift(RIGHT * 3.6 + DOWN * 0.6)
        self.play(Create(zp))
        u = zp.get_x_unit_size()
        L0 = LEFT * 6.2 + DOWN * 2.6
        s = ValueTracker(0.0)

        def pos(th, sv):
            # 线段上的点 -> 圆上的点
            a = L0 + RIGHT * th / (2 * PI) * 5.5
            b = zp.n2p(np.exp(1j * th))
            return (1 - sv) * a + sv * b

        band = always_redraw(lambda: VMobject(stroke_color=C_Z, stroke_width=5).set_points_smoothly(
            [pos(t, s.get_value()) for t in np.linspace(0, 2 * PI, 120)]))
        ticks = always_redraw(lambda: VGroup(*[Dot(pos(t, s.get_value()), radius=0.05, color=C_W) for t in np.linspace(0, 2 * PI, 9)[:-1]]))
        l0 = MT(r"0", 26).next_to(L0, DOWN, buff=0.1)
        l2 = MT(r"2\pi", 26).next_to(L0 + RIGHT * 5.5, DOWN, buff=0.1)
        tl = MT(r"\theta", 28).next_to(L0 + RIGHT * 2.75, DOWN, buff=0.1)
        self.play(Create(band), FadeIn(ticks), FadeIn(l0), FadeIn(l2), FadeIn(tl))
        self.play(s.animate.set_value(1), FadeOut(l0), FadeOut(l2), FadeOut(tl), run_time=3)
        uc = MT(r"|z|=1", 28, C_Z).next_to(zp.n2p(0.75 + 0.75j), UR, buff=0.05)
        self.play(FadeIn(uc))
        rules = VGroup(
            MT(r"z=e^{i\theta},\qquad \mathrm d\theta=\frac{\mathrm dz}{iz}", 32),
            MT(r"\cos\theta=\frac{z+z^{-1}}{2},\qquad \sin\theta=\frac{z-z^{-1}}{2i}", 32),
            MT(r"\int_0^{2\pi}R(\cos\theta,\sin\theta)\,\mathrm d\theta=\oint_{|z|=1}R\Big(\tfrac{z+z^{-1}}{2},\tfrac{z-z^{-1}}{2i}\Big)\frac{\mathrm dz}{iz}", 28, C_Z),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.35)
        fitw(rules, 6.6).move_to(LEFT * 3.2 + UP * 0.9)
        self.play(LaggedStart(*[FadeIn(r) for r in rules], lag_ratio=0.6), run_time=2.5)
        self.wait(1)
        ex = VGroup(
            CT(r"例：$\displaystyle\int_0^{2\pi}\frac{\mathrm d\theta}{5+3\sin\theta}=\oint_{|z|=1}\frac{2\,\mathrm dz}{3z^2+10iz-3}$", 26),
            CT(r"极点 $-\tfrac i3$（圆内）与 $-3i$（圆外）", 26),
            CT(r"$=2\pi i\cdot\dfrac{2}{6z+10i}\Big|_{z=-i/3}=2\pi i\cdot\dfrac{1}{4i}=\dfrac{\pi}{2}$", 26, C_HL),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.25)
        fitw(ex, 6.6).next_to(rules, DOWN, buff=0.45, aligned_edge=LEFT)
        p_in = Dot(zp.n2p(-1j / 3), color=RED)
        p_in_l = MT(r"-\tfrac i3", 24, RED).next_to(p_in, RIGHT, buff=0.06)
        arrow_out = Arrow(zp.n2p(-1.2j), zp.n2p(-1.6j), buff=0, color=GREY_B)
        out_l = CT(r"$-3i$ 在下方圆外", 20, GREY_B).next_to(arrow_out, LEFT, buff=0.05)
        self.play(FadeIn(ex[0]))
        self.play(FadeIn(ex[1]), FadeIn(p_in), FadeIn(p_in_l), GrowArrow(arrow_out), FadeIn(out_l))
        self.play(FadeIn(ex[2]))
        self.wait(3)
