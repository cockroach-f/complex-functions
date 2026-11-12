"""第八章 傅里叶变换 —— manim 动画源码（LaTeX 排版）."""
import numpy as np

from common import *  # noqa: F401,F403


def fitw(m, w):
    if m.width > w:
        m.scale_to_fit_width(w)
    return m


def ax(xr, yr, xl, yl, **kw):
    a = Axes(
        x_range=list(xr), y_range=list(yr), x_length=xl, y_length=yl,
        axis_config={"stroke_width": 2, "include_tip": True, "tip_length": 0.18, "include_ticks": False},
        y_axis_config={"include_tip": False},
        **kw,
    )
    return a


# --------------------------------------------------------------------------
# 8.1 周期 T → ∞：离散谱线变密，傅里叶级数 → 傅里叶积分
# --------------------------------------------------------------------------
class C8_SeriesToIntegral(Scene):
    def construct(self):
        title = title_bar(r"周期 $T\to+\infty$：离散频谱 $\to$ 连续频谱，级数求和 $\to$ 积分")
        self.play(FadeIn(title))
        ta = ax((-10, 10, 1), (-0.2, 1.4, 1), 8.0, 1.6).move_to(LEFT * 2.3 + UP * 1.55)
        tl = MT("t", 24, GREY_B).next_to(ta.x_axis.get_end(), DOWN, buff=0.1)
        fa = ax((-8.5, 8.5, 1), (-0.7, 2.3, 1), 8.0, 3.0).move_to(LEFT * 2.3 + DOWN * 1.75)
        fl = MT(r"\omega", 26, GREY_B).next_to(fa.x_axis.get_end(), DOWN, buff=0.1)
        self.play(Create(ta), Create(fa), FadeIn(tl), FadeIn(fl))
        env = fa.plot(lambda w: 2 * np.sinc(w / np.pi), x_range=[-8.4, 8.4, 0.02], color=C_HL, stroke_width=2)
        env = DashedVMobject(env, num_dashes=90)
        envl = MT(r"F(\omega)=\frac{2\sin\omega}{\omega}", 28, C_HL).move_to(fa.c2p(5.6, 1.6))

        def pulses(T):
            g = VGroup()
            ks = range(-int(12 / T) - 1, int(12 / T) + 2)
            for k in ks:
                c = k * T
                if c - 1 < -10 or c + 1 > 10:
                    continue
                pts = [ta.c2p(c - 1, 0), ta.c2p(c - 1, 1), ta.c2p(c + 1, 1), ta.c2p(c + 1, 0)]
                g.add(VMobject(color=C_Z if k == 0 else C_Z, stroke_width=3, stroke_opacity=1 if k == 0 else 0.55).set_points_as_corners(pts))
            return g

        def stems(T):
            w0 = 2 * PI / T
            g = VGroup()
            n = int(8.4 / w0)
            for k in range(-n, n + 1):
                w = k * w0
                h = 2.0 if k == 0 else 2 * np.sin(w) / w
                g.add(Line(fa.c2p(w, 0), fa.c2p(w, h), color=C_W, stroke_width=2.5 if T < 20 else 1.8))
                g.add(Dot(fa.c2p(w, h), radius=0.035 if T < 20 else 0.025, color=C_W))
            return g

        side = VGroup(
            CT(r"$f_T(t)$：宽 $2$ 的方波，周期 $T$", 24, C_Z),
            CT(r"谱线位置 $\omega_n=n\omega_0$，$\omega_0=\dfrac{2\pi}{T}$", 24, C_W),
            CT(r"谱线高度 $T c_n=\displaystyle\int_{-T/2}^{T/2}f_T(t)e^{-j\omega_n t}\,\mathrm dt$", 24, C_W),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.3)
        fitw(side, 4.6).to_edge(RIGHT, buff=0.3).shift(UP * 1.4)
        self.play(FadeIn(side))
        cur_p = cur_s = lab = None
        for T in (4, 8, 16, 32):
            p, s = pulses(T), stems(T)
            l = CT(rf"$T={T}$，$\omega_0=\dfrac{{2\pi}}{{{T}}}\approx{2*np.pi/T:.2f}$", 28, WHITE).next_to(side, DOWN, buff=0.5, aligned_edge=LEFT)
            if cur_p is None:
                self.play(Create(p), FadeIn(l))
                self.play(LaggedStart(*[GrowFromPoint(m, fa.c2p(0, 0)) for m in s], lag_ratio=0.02), run_time=1.5)
            else:
                self.play(ReplacementTransform(cur_p, p), ReplacementTransform(cur_s, s), ReplacementTransform(lab, l), run_time=1.6)
            cur_p, cur_s, lab = p, s, l
            if T == 8:
                self.play(Create(env), FadeIn(envl))
            self.wait(1.0)
        note = VGroup(
            CT(r"谱线间隔 $\Delta\omega=\omega_0\to0$，顶点落在同一条包络 $F(\omega)$ 上", 24, C_HL),
            MT(r"f_T(t)=\frac{1}{2\pi}\sum_{n}F(\omega_n)e^{j\omega_n t}\Delta\omega", 28),
            MT(r"\xrightarrow{\ T\to\infty\ }\ f(t)=\frac{1}{2\pi}\int_{-\infty}^{+\infty}F(\omega)e^{j\omega t}\,\mathrm d\omega", 28, C_Z),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.28)
        fitw(note, 4.4).next_to(lab, DOWN, buff=0.4, aligned_edge=LEFT)
        self.play(LaggedStart(*[FadeIn(n) for n in note], lag_ratio=0.6), run_time=2.5)
        self.wait(3)


# --------------------------------------------------------------------------
# 8.2 δ 函数：面积为 1 的窄脉冲的极限；频谱趋于常数 1；筛选性质
# --------------------------------------------------------------------------
class C8_DeltaLimit(Scene):
    def construct(self):
        title = title_bar(r"$\delta(t)=\lim\limits_{\varepsilon\to0}\delta_\varepsilon(t)$：宽 $\varepsilon$、高 $\frac1\varepsilon$、面积恒为 $1$")
        self.play(FadeIn(title))
        ta = ax((-2, 2, 1), (0, 8.5, 1), 5.4, 5.0).move_to(LEFT * 3.5 + DOWN * 0.55)
        fa = ax((-30, 30, 10), (-0.3, 1.3, 1), 5.6, 2.4).move_to(RIGHT * 3.3 + UP * 0.85)
        self.play(Create(ta), Create(fa))
        self.add(MT("t", 24, GREY_B).next_to(ta.x_axis.get_end(), DOWN, buff=0.1),
                 MT(r"\omega", 24, GREY_B).next_to(fa.x_axis.get_end(), DOWN, buff=0.1))
        one = DashedLine(fa.c2p(-30, 1), fa.c2p(30, 1), color=GREY_B, stroke_width=1.5)
        onel = MT("1", 22, GREY_B).next_to(fa.c2p(-30, 1), LEFT, buff=0.1)
        self.add(one, onel)
        eps = ValueTracker(2.0)

        def rect():
            e = eps.get_value()
            r = Polygon(ta.c2p(-e / 2, 0), ta.c2p(-e / 2, 1 / e), ta.c2p(e / 2, 1 / e), ta.c2p(e / 2, 0),
                        color=C_Z, fill_color=C_Z, fill_opacity=0.45, stroke_width=3)
            return r

        R = always_redraw(rect)
        spec = always_redraw(lambda: fa.plot(lambda w: np.sinc(w * eps.get_value() / 2 / np.pi), x_range=[-30, 30, 0.1], color=C_W, stroke_width=3))
        el = always_redraw(lambda: VGroup(MT(r"\varepsilon=", 28), DecimalNumber(eps.get_value(), num_decimal_places=3, font_size=28)).arrange(RIGHT, buff=0.1).move_to(ta.c2p(1.25, 7.6)))
        area = CT(r"面积 $=\varepsilon\cdot\frac1\varepsilon=1$", 26, C_Z).move_to(ta.c2p(1.2, 6.5))
        sl = MT(r"\mathcal F[\delta_\varepsilon](\omega)=\frac{\sin(\omega\varepsilon/2)}{\omega\varepsilon/2}", 28, C_W).next_to(fa, UP, buff=0.15)
        self.play(Create(R), Create(spec), FadeIn(el), FadeIn(area), FadeIn(sl))
        self.play(eps.animate.set_value(0.125), run_time=6, rate_func=smooth)
        res = VGroup(
            CT(r"脉冲越窄，频谱越平：$\mathcal F[\delta(t)]=1$", 26, C_HL),
            CT(r"（均匀频谱 / 白色频谱：所有频率等幅）", 22, GREY_B),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.15).next_to(fa, DOWN, buff=0.25)
        self.play(FadeIn(res))
        self.wait(1.5)
        # 筛选性质
        self.play(eps.animate.set_value(1.6), run_time=1.5)
        f = lambda t: 2 + np.cos(3 * t)
        fg = ta.plot(f, x_range=[-2, 2, 0.02], color=C_RE, stroke_width=3)
        flab = MT(r"f(t)=2+\cos 3t", 24, C_RE).move_to(ta.c2p(-1.05, 3.6))
        self.play(Create(fg), FadeIn(flab))
        val = always_redraw(lambda: VGroup(
            MT(r"\int_{-\infty}^{+\infty}\delta_\varepsilon(t)f(t)\,\mathrm dt=", 28),
            DecimalNumber(2 + np.sin(1.5 * eps.get_value()) / (1.5 * eps.get_value()), num_decimal_places=4, font_size=30, color=C_HL),
        ).arrange(RIGHT, buff=0.12).next_to(res, DOWN, buff=0.3))
        self.play(FadeIn(val))
        self.play(eps.animate.set_value(0.13), run_time=5, rate_func=smooth)
        sift = CT(r"筛选性质：$\displaystyle\int_{-\infty}^{+\infty}\delta(t-t_0)f(t)\,\mathrm dt=f(t_0)$，这里 $f(0)=3$", 26, C_Z)
        fitw(sift, 6.4).next_to(val, DOWN, buff=0.3)
        self.play(FadeIn(sift))
        self.wait(3)


# --------------------------------------------------------------------------
# 8.3 相似性质：时域压缩 ⇔ 频域展宽
# --------------------------------------------------------------------------
class C8_Scaling(Scene):
    def construct(self):
        title = title_bar(r"相似性质 $\mathcal F[f(at)]=\dfrac{1}{|a|}F\!\left(\dfrac{\omega}{a}\right)$：时域越窄，频谱越宽")
        self.play(FadeIn(title))
        ta = ax((-4, 4, 1), (0, 1.25, 1), 5.6, 3.2).move_to(LEFT * 3.4 + DOWN * 0.3)
        fa = ax((-12, 12, 4), (0, 3.8, 1), 5.6, 3.2).move_to(RIGHT * 3.4 + DOWN * 0.3)
        self.play(Create(ta), Create(fa))
        self.add(MT("t", 24, GREY_B).next_to(ta.x_axis.get_end(), DOWN, buff=0.1),
                 MT(r"\omega", 24, GREY_B).next_to(fa.x_axis.get_end(), DOWN, buff=0.1))
        a = ValueTracker(1.0)
        g = always_redraw(lambda: ta.plot(lambda t: np.exp(-(a.get_value() * t) ** 2), x_range=[-4, 4, 0.02], color=C_Z, stroke_width=3.5))
        G = always_redraw(lambda: fa.plot(lambda w: np.sqrt(np.pi) / a.get_value() * np.exp(-w ** 2 / (4 * a.get_value() ** 2)), x_range=[-12, 12, 0.04], color=C_W, stroke_width=3.5))
        gl = MT(r"f(at)=e^{-a^2t^2}", 30, C_Z).next_to(ta, UP, buff=0.25)
        Gl = MT(r"\frac{1}{|a|}F\!\left(\frac{\omega}{a}\right)=\frac{\sqrt\pi}{a}e^{-\omega^2/(4a^2)}", 30, C_W).next_to(fa, UP, buff=0.15)
        al = always_redraw(lambda: VGroup(MT("a=", 32, C_HL), DecimalNumber(a.get_value(), num_decimal_places=2, font_size=32, color=C_HL)).arrange(RIGHT, buff=0.1).to_edge(DOWN, buff=1.3))
        # 半高宽标记
        def wbar(axes, half, y, col):
            return Line(axes.c2p(-half, y), axes.c2p(half, y), color=col, stroke_width=5)
        tb = always_redraw(lambda: wbar(ta, np.sqrt(np.log(2)) / a.get_value(), 0.5, C_HL))
        fb = always_redraw(lambda: wbar(fa, min(2 * a.get_value() * np.sqrt(np.log(2)), 12), np.sqrt(np.pi) / a.get_value() / 2, C_HL))
        self.play(Create(g), Create(G), FadeIn(gl), FadeIn(Gl), FadeIn(al))
        self.play(FadeIn(tb), FadeIn(fb))
        self.wait(0.5)
        self.play(a.animate.set_value(2.5), run_time=3)
        c1 = CT(r"$a>1$：信号压缩（播放加快），频带展宽、幅度降低", 26, C_HL).to_edge(DOWN, buff=0.5)
        self.play(FadeIn(c1))
        self.wait(1)
        self.play(a.animate.set_value(0.5), run_time=4)
        c2 = CT(r"$0<a<1$：信号拉长（播放变慢），频带变窄、幅度升高", 26, C_HL).to_edge(DOWN, buff=0.5)
        self.play(FadeTransform(c1, c2))
        self.wait(1)
        self.play(a.animate.set_value(1.0), run_time=2)
        c3 = CT(r"时宽 $\times$ 带宽 $=$ 常数：不能同时任意窄（测不准关系）", 26, C_Z).to_edge(DOWN, buff=0.5)
        self.play(FadeTransform(c2, c3))
        self.wait(3)


# --------------------------------------------------------------------------
# 8.3 卷积：翻转、平移、相乘、积分
# --------------------------------------------------------------------------
class C8_Convolution(Scene):
    def construct(self):
        title = title_bar(r"卷积 $f*g(t)=\displaystyle\int_{-\infty}^{+\infty}f(\tau)g(t-\tau)\,\mathrm d\tau$：翻转 $\to$ 平移 $\to$ 相乘 $\to$ 求面积")
        self.play(FadeIn(title))
        A = ax((-2.5, 5, 1), (0, 1.25, 1), 7.2, 2.4).move_to(LEFT * 2.95 + UP * 1.0)
        B = ax((-2.5, 5, 1), (0, 0.32, 0.1), 7.2, 2.0).move_to(LEFT * 2.95 + DOWN * 2.2)
        self.play(Create(A), Create(B))
        self.add(MT(r"\tau", 26, GREY_B).next_to(A.x_axis.get_end(), DOWN, buff=0.1),
                 MT("t", 26, GREY_B).next_to(B.x_axis.get_end(), DOWN, buff=0.1))
        f = lambda x: np.where(x >= 0, np.exp(-x), 0.0)
        g = lambda x: np.where(x >= 0, np.exp(-2 * x), 0.0)
        h = lambda t: np.where(t >= 0, np.exp(-t) - np.exp(-2 * t), 0.0)
        fg = VGroup(A.plot(lambda x: 0, x_range=[-2.5, 0], color=C_Z, stroke_width=3.5),
                    Line(A.c2p(0, 0), A.c2p(0, 1), color=C_Z, stroke_width=3.5),
                    A.plot(lambda x: np.exp(-x), x_range=[0, 4.9], color=C_Z, stroke_width=3.5))
        fl = MT(r"f(\tau)=e^{-\tau}u(\tau)", 28, C_Z).move_to(RIGHT * 2.0 + UP * 2.5, aligned_edge=LEFT)
        self.play(Create(fg), FadeIn(fl))
        steps = VGroup(
            CT(r"\textcircled{1} 翻转：$g(\tau)\to g(-\tau)$", 28),
            CT(r"\textcircled{2} 平移：$g(t-\tau)$", 28),
            CT(r"\textcircled{3} 相乘：$f(\tau)g(t-\tau)$", 28),
            CT(r"\textcircled{4} 积分：阴影面积 $=f*g(t)$", 28),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.28)
        fitw(steps, 4.8).move_to(RIGHT * 2.0 + UP * 0.3, aligned_edge=LEFT)
        g0 = A.plot(lambda x: float(g(x)), x_range=[0.001, 4.9], color=C_W, stroke_width=3.5)
        g0l = MT(r"g(\tau)=e^{-2\tau}u(\tau)", 28, C_W).next_to(fl, DOWN, buff=0.2, aligned_edge=LEFT)
        self.play(Create(g0), FadeIn(g0l))
        self.wait(0.5)
        self.play(FadeIn(steps[0]))
        gflip = A.plot(lambda x: float(g(-x)), x_range=[-2.5, -0.001], color=C_W, stroke_width=3.5)
        self.play(Transform(g0, gflip), run_time=1.5)
        t = ValueTracker(0.0)
        self.remove(g0)

        def gshift():
            tt = t.get_value()
            lo = max(-2.5, tt - 6)
            return A.plot(lambda x: float(g(tt - x)), x_range=[lo, min(tt, 4.99) - 1e-3] if tt > -2.49 else [-2.5, -2.49], color=C_W, stroke_width=3.5)

        gs = always_redraw(gshift)
        gsl = always_redraw(lambda: MT(r"g(t-\tau)", 24, C_W).next_to(A.c2p(min(t.get_value(), 4.5), 1.0), UP, buff=0.05))
        self.add(gs, gsl)
        self.play(FadeIn(steps[1]))
        self.play(t.animate.set_value(-1.5), run_time=1.2)
        self.play(FadeIn(steps[2]), FadeIn(steps[3]))

        def shade():
            tt = t.get_value()
            if tt <= 0.02:
                return VMobject()
            xs = np.linspace(0, tt, 80)
            pts = [A.c2p(0, 0)] + [A.c2p(x, float(f(x) * g(tt - x))) for x in xs] + [A.c2p(tt, 0)]
            return Polygon(*pts, color=C_HL, fill_color=C_HL, fill_opacity=0.6, stroke_width=1.5)

        sh = always_redraw(shade)
        tr = always_redraw(lambda: B.plot(lambda x: float(h(x)), x_range=[-2.5, max(t.get_value(), -2.49)], color=C_HL, stroke_width=3.5))
        dot = always_redraw(lambda: Dot(B.c2p(t.get_value(), float(h(t.get_value()))), color=C_HL, radius=0.06))
        vl = always_redraw(lambda: DashedLine(A.c2p(t.get_value(), 0), B.c2p(t.get_value(), 0), color=GREY_B, stroke_width=1.5))
        self.add(sh, tr, dot, vl)
        self.play(t.animate.set_value(0.0), run_time=1.5, rate_func=linear)
        z0 = CT(r"$t<0$：两者不重叠，$f*g=0$", 26, GREY_B).next_to(steps, DOWN, buff=0.45, aligned_edge=LEFT)
        self.play(FadeIn(z0))
        self.play(t.animate.set_value(4.8), run_time=7, rate_func=linear)
        rl = MT(r"f*g(t)=\left(e^{-t}-e^{-2t}\right)u(t)", 28, C_HL).next_to(z0, DOWN, buff=0.35, aligned_edge=LEFT)
        fitw(rl, 4.8)
        ct = CT(r"卷积定理：$\mathcal F[f*g]=\dfrac{1}{(1+j\omega)(2+j\omega)}$", 26, C_Z)
        fitw(ct, 4.8).next_to(rl, DOWN, buff=0.3, aligned_edge=LEFT)
        self.play(FadeIn(rl), FadeIn(ct))
        self.wait(3)
