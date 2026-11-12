# 公式速查

随各章完成逐步补充。

## 第 1 章 复数与复变函数

| 内容 | 公式 |
|---|---|
| 代数式 / 三角式 / 指数式 | $z=x+iy=r(\cos\theta+i\sin\theta)=re^{i\theta}$ |
| 模 | $\lvert z\rvert=\sqrt{x^2+y^2}=\sqrt{z\bar z}$ |
| 主辐角 | $\arg z\in(-\pi,\pi]$，$\operatorname{Arg}z=\arg z+2k\pi$ |
| 实部虚部 | $\operatorname{Re}z=\frac{z+\bar z}{2}$，$\operatorname{Im}z=\frac{z-\bar z}{2i}$ |
| 欧拉公式 | $e^{i\theta}=\cos\theta+i\sin\theta$ |
| 乘法 | $z_1z_2=r_1r_2e^{i(\theta_1+\theta_2)}$ |
| 除法 | $\frac{z_1}{z_2}=\frac{r_1}{r_2}e^{i(\theta_1-\theta_2)}$ |
| 棣莫弗 | $(\cos\theta+i\sin\theta)^n=\cos n\theta+i\sin n\theta$ |
| $n$ 次方根 | $\sqrt[n]{z}=\sqrt[n]{r}\,e^{i\frac{\theta+2k\pi}{n}},\ k=0,1,\dots,n-1$ |
| 三角不等式 | $\big\lvert\lvert z_1\rvert-\lvert z_2\rvert\big\rvert\le\lvert z_1+z_2\rvert\le\lvert z_1\rvert+\lvert z_2\rvert$ |
| 圆 | $\lvert z-z_0\rvert=R$；参数式 $z=z_0+Re^{i\theta}$ |
| 线段 | $z=z_1+t(z_2-z_1),\ t\in[0,1]$ |
| 一般圆 / 直线 | $Az\bar z+\bar\beta z+\beta\bar z+C=0$（$A,C\in\mathbb R$） |
| 复变函数 | $f(z)=u(x,y)+iv(x,y)$；$\lim f=A\iff\lim u=u_0,\ \lim v=v_0$ |

## 第 2 章 解析函数

| 内容 | 公式 |
|---|---|
| 导数 | $f'(z_0)=\lim\limits_{\Delta z\to0}\frac{f(z_0+\Delta z)-f(z_0)}{\Delta z}$ |
| C-R 方程 | $u_x=v_y,\quad u_y=-v_x$ |
| 导数公式 | $f'(z)=u_x+iv_x=v_y-iu_y$ |
| 极坐标 C-R | $u_r=\frac1r v_\theta,\quad v_r=-\frac1r u_\theta$ |
| 调和 | $\Delta u=u_{xx}+u_{yy}=0$ |
| 构造解析函数 | $f'(z)=u_x-iu_y$，令 $x=z,\ y=0$ 后积分 |
| 指数 | $e^z=e^x(\cos y+i\sin y)$，周期 $2\pi i$ |
| 对数 | $\operatorname{Ln}z=\ln\lvert z\rvert+i\arg z+2k\pi i$，$(\ln z)'=\frac1z$ |
| 乘幂 | $a^b=e^{b\operatorname{Ln}a}$ |
| 三角 | $\sin z=\frac{e^{iz}-e^{-iz}}{2i}$，$\cos z=\frac{e^{iz}+e^{-iz}}{2}$ |
| 双曲 | $\sinh z=\frac{e^z-e^{-z}}{2}$，$\cosh z=\frac{e^z+e^{-z}}{2}$，$\cos iy=\cosh y$，$\sin iy=i\sinh y$ |
| 反三角 | $\operatorname{Arcsin}z=-i\operatorname{Ln}\big(iz+\sqrt{1-z^2}\big)$，$\operatorname{Arctan}z=-\frac i2\operatorname{Ln}\frac{1+iz}{1-iz}$ |

## 第 3 章 复变函数的积分

| 内容 | 公式 |
|---|---|
| 参数法 | $\int_Cf(z)\,\mathrm dz=\int_\alpha^\beta f(z(t))\,z'(t)\,\mathrm dt$ |
| 实部虚部 | $\int_Cf\,\mathrm dz=\int_Cu\,\mathrm dx-v\,\mathrm dy+i\int_Cv\,\mathrm dx+u\,\mathrm dy$ |
| 长大不等式 | $\left\lvert\int_Cf\,\mathrm dz\right\rvert\le\int_C\lvert f\rvert\,\mathrm ds\le ML$ |
| 关键积分 | $\oint_{\lvert z-z_0\rvert=r}\frac{\mathrm dz}{(z-z_0)^n}=2\pi i\ (n=1)$，$0\ (n\ne1)$ |
| 柯西定理 | 单连通域内解析 ⇒ $\oint_Cf\,\mathrm dz=0$ |
| 复合闭路 | $\oint_Cf\,\mathrm dz=\sum_k\oint_{C_k}f\,\mathrm dz$ |
| 牛顿-莱布尼茨 | $\int_{z_0}^{z_1}f\,\mathrm dz=G(z_1)-G(z_0)$，$G'=f$ |
| 柯西积分公式 | $f(z_0)=\frac{1}{2\pi i}\oint_C\frac{f(z)}{z-z_0}\,\mathrm dz$ |
| 高阶导数 | $f^{(n)}(z_0)=\frac{n!}{2\pi i}\oint_C\frac{f(z)}{(z-z_0)^{n+1}}\,\mathrm dz$ |
| 平均值 | $f(z_0)=\frac1{2\pi}\int_0^{2\pi}f(z_0+Re^{i\theta})\,\mathrm d\theta$ |
| 柯西不等式 | $\lvert f^{(n)}(z_0)\rvert\le\frac{n!M}{R^n}$ |
| 刘维尔 | 全平面解析且有界 ⇒ 常数 |

## 第 4 章 解析函数的级数表示

| 内容 | 公式 |
|---|---|
| 几何级数 | $\sum_{n=0}^\infty z^n=\frac{1}{1-z}$，$\lvert z\rvert<1$ |
| 收敛半径 | $R=\frac1\lambda$，$\lambda=\lim\left\lvert\frac{c_{n+1}}{c_n}\right\rvert=\lim\sqrt[n]{\lvert c_n\rvert}$ |
| 泰勒系数 | $c_n=\frac{f^{(n)}(z_0)}{n!}=\frac{1}{2\pi i}\oint_C\frac{f(\zeta)}{(\zeta-z_0)^{n+1}}\,\mathrm d\zeta$ |
| 收敛半径（解析函数） | $R=$ 展开中心到最近奇点的距离 |
| $e^z$ | $\sum\frac{z^n}{n!}$，$\lvert z\rvert<\infty$ |
| $\sin z,\ \cos z$ | $\sum\frac{(-1)^nz^{2n+1}}{(2n+1)!}$，$\sum\frac{(-1)^nz^{2n}}{(2n)!}$ |
| $\ln(1+z)$ | $\sum_{n\ge0}\frac{(-1)^nz^{n+1}}{n+1}$，$\lvert z\rvert<1$ |
| $(1+z)^\alpha$ | $1+\alpha z+\frac{\alpha(\alpha-1)}{2!}z^2+\cdots$，$\lvert z\rvert<1$ |
| 洛朗级数 | $f=\sum_{n=-\infty}^{\infty}c_n(z-z_0)^n$，$c_n=\frac1{2\pi i}\oint_C\frac{f(\zeta)}{(\zeta-z_0)^{n+1}}\,\mathrm d\zeta$ |
| 负幂展开 | $\frac{1}{z-a}=\frac{1}{z}\sum_{n\ge0}\left(\frac az\right)^n$，$\lvert z\rvert>\lvert a\rvert$ |
| 积分与 $c_{-1}$ | $\oint_Cf\,\mathrm dz=2\pi i\,c_{-1}$ |

## 第 5 章 留数及其应用

| 内容 | 公式 |
|---|---|
| 奇点分类（洛朗展开） | 可去：无负幂；$m$ 阶极点：负幂到 $(z-z_0)^{-m}$ 为止；本性：无穷多负幂 |
| 奇点分类（极限） | $\lim f$ 有限 ⇔ 可去；$=\infty$ ⇔ 极点；不存在且不为 $\infty$ ⇔ 本性 |
| 零点与极点 | $z_0$ 是 $f$ 的 $m$ 阶零点 ⇔ $z_0$ 是 $\frac1f$ 的 $m$ 阶极点 |
| 留数 | $\operatorname{Res}[f,z_0]=c_{-1}=\frac{1}{2\pi i}\oint_Cf\,\mathrm dz$ |
| 留数定理 | $\oint_Cf\,\mathrm dz=2\pi i\sum_{k=1}^n\operatorname{Res}[f,z_k]$ |
| 规则 I（一阶极点） | $\operatorname{Res}[f,z_0]=\lim_{z\to z_0}(z-z_0)f(z)$ |
| 规则 III（$m$ 阶极点） | $\operatorname{Res}[f,z_0]=\frac{1}{(m-1)!}\lim_{z\to z_0}\frac{\mathrm d^{m-1}}{\mathrm dz^{m-1}}\left[(z-z_0)^mf(z)\right]$ |
| 规则 II（$f=\frac PQ$） | $P(z_0)\ne0$，$Q(z_0)=0$，$Q'(z_0)\ne0$ ⇒ $\operatorname{Res}[f,z_0]=\frac{P(z_0)}{Q'(z_0)}$ |
| 无穷远点的留数 | $\operatorname{Res}[f,\infty]=-c_{-1}=-\operatorname{Res}\left[f\left(\tfrac1z\right)\tfrac1{z^2},0\right]$ |
| 留数和为零 | $\sum_{k}\operatorname{Res}[f,z_k]+\operatorname{Res}[f,\infty]=0$ |
| 公式 I：$R(\cos\theta,\sin\theta)$ | $z=e^{i\theta}$，$\cos\theta=\frac{z^2+1}{2z}$，$\sin\theta=\frac{z^2-1}{2iz}$，$\mathrm d\theta=\frac{\mathrm dz}{iz}$；$\int_0^{2\pi}=2\pi i\sum_{\lvert z_k\rvert<1}\operatorname{Res}$ |
| 公式 II：$\frac{P}{Q}$ | $\deg Q\ge\deg P+2$，$Q$ 在实轴上无零点 ⇒ $\int_{-\infty}^{\infty}R(x)\,\mathrm dx=2\pi i\sum_{\operatorname{Im}z_k>0}\operatorname{Res}[R,z_k]$ |
| 公式 III：$R(x)e^{iax}$，$a>0$ | $\deg Q\ge\deg P+1$ ⇒ $\int_{-\infty}^{\infty}R(x)e^{iax}\,\mathrm dx=2\pi i\sum_{\operatorname{Im}z_k>0}\operatorname{Res}[R(z)e^{iaz},z_k]$，再取实部、虚部 |
| 约当引理 | $\left\lvert\int_{C_R}R(z)e^{iaz}\,\mathrm dz\right\rvert\to0$（$R(z)\to0$ 一致） |
| 狄利克雷积分 | $\int_0^{+\infty}\frac{\sin x}{x}\,\mathrm dx=\frac\pi2$ |

## 第 6 章 共形映射

| 内容 | 公式 |
|---|---|
| 保角性 | $f'(z_0)\ne0$：所有曲线转过 $\arg f'(z_0)$，伸缩率 $\lvert f'(z_0)\rvert$ |
| 分式线性映射 | $w=\frac{az+b}{cz+d}$，$ad-bc\ne0$ |
| 分解 | 平移 $z+b$、旋转伸缩 $az$、反演 $\frac1z$ 的复合 |
| 交比（三点定映射） | $\frac{w-w_1}{w-w_2}\cdot\frac{w_3-w_2}{w_3-w_1}=\frac{z-z_1}{z-z_2}\cdot\frac{z_3-z_2}{z_3-z_1}$ |
| 性质 | 保角、保圆（直线是过 $\infty$ 的圆）、保对称点 |
| 上半平面 → 单位圆盘 | $w=e^{i\theta}\frac{z-\lambda}{z-\bar\lambda}$，$\operatorname{Im}\lambda>0$ |
| 单位圆盘 → 单位圆盘 | $w=e^{i\theta}\frac{z-\alpha}{1-\bar\alpha z}$，$\lvert\alpha\rvert<1$ |
| 上半平面 → 上半平面 | $w=\frac{az+b}{cz+d}$，$a,b,c,d\in\mathbb R$，$ad-bc>0$ |
| 幂函数 | $w=z^n$：角形 $0<\arg z<\alpha$ → $0<\arg w<n\alpha$（$n\alpha\le2\pi$） |
| 指数函数 | $w=e^z$：带形 $0<\operatorname{Im}z<a$ → 角形 $0<\arg w<a$（$a\le2\pi$） |
| 两圆弧围成的区域 | 用 $\frac{z-z_1}{z-z_2}$ 把交点 $z_1,z_2$ 送到 $0,\infty$，得到角形域 |
| 构造步骤 | 角形 / 带形 → 上半平面（$z^n$ 或 $e^z$）→ 单位圆盘（$\frac{\eta-\lambda}{\eta-\bar\lambda}$） |

## 第 8 章 傅里叶变换

| 内容 | 公式 |
|---|---|
| 傅里叶级数（复指数形式） | $f_T(t)=\sum_{n=-\infty}^{+\infty}c_ne^{jn\omega_0t}$，$c_n=\frac1T\int_{-T/2}^{T/2}f_T(t)e^{-jn\omega_0t}\,\mathrm dt$，$\omega_0=\frac{2\pi}T$ |
| 傅里叶变换 | $F(\omega)=\int_{-\infty}^{+\infty}f(t)e^{-j\omega t}\,\mathrm dt$ |
| 逆变换 | $f(t)=\frac1{2\pi}\int_{-\infty}^{+\infty}F(\omega)e^{j\omega t}\,\mathrm d\omega$ |
| 筛选性质 | $\int_{-\infty}^{+\infty}\delta(t-t_0)f(t)\,\mathrm dt=f(t_0)$ |
| $\delta$ 与 $u$ | $\int_{-\infty}^t\delta(\tau)\,\mathrm d\tau=u(t)$，$u'(t)=\delta(t)$，$\delta(-t)=\delta(t)$ |
| 重要公式 | $\int_{-\infty}^{+\infty}e^{-j\omega t}\,\mathrm dt=2\pi\delta(\omega)$ |
| 常用变换对 | $\delta\leftrightarrow1$，$1\leftrightarrow2\pi\delta(\omega)$，$e^{j\omega_0t}\leftrightarrow2\pi\delta(\omega-\omega_0)$，$u(t)\leftrightarrow\frac1{j\omega}+\pi\delta(\omega)$ |
| 指数信号 | $e^{-\beta t}u(t)\leftrightarrow\frac1{\beta+j\omega}$，$e^{-\beta\lvert t\rvert}\leftrightarrow\frac{2\beta}{\beta^2+\omega^2}$（$\beta>0$） |
| 矩形 / 抽样 | $1\ (\lvert t\rvert\le\delta)\leftrightarrow\frac{2\sin\delta\omega}{\omega}$；$\frac{\sin\alpha t}{\pi t}\leftrightarrow1\ (\lvert\omega\rvert\le\alpha)$ |
| 正弦、余弦 | $\cos\omega_0t\leftrightarrow\pi[\delta(\omega+\omega_0)+\delta(\omega-\omega_0)]$，$\sin\omega_0t\leftrightarrow j\pi[\delta(\omega+\omega_0)-\delta(\omega-\omega_0)]$ |
| 高斯 | $e^{-t^2}\leftrightarrow\sqrt\pi\,e^{-\omega^2/4}$ |
| 位移 | $f(t-t_0)\leftrightarrow e^{-j\omega t_0}F(\omega)$，$e^{j\omega_0t}f(t)\leftrightarrow F(\omega-\omega_0)$ |
| 相似 | $f(at)\leftrightarrow\frac1{\lvert a\rvert}F\left(\frac\omega a\right)$ |
| 微分 | $f^{(n)}(t)\leftrightarrow(j\omega)^nF(\omega)$；$t^nf(t)\leftrightarrow j^nF^{(n)}(\omega)$ |
| 积分 | $\int_{-\infty}^tf\,\mathrm d\tau\leftrightarrow\frac{F(\omega)}{j\omega}$（$F(0)=0$ 时） |
| 卷积 | $f_1*f_2=\int_{-\infty}^{+\infty}f_1(\tau)f_2(t-\tau)\,\mathrm d\tau$ |
| 卷积定理 | $f_1*f_2\leftrightarrow F_1F_2$，$f_1f_2\leftrightarrow\frac1{2\pi}F_1*F_2$ |
| 帕塞瓦尔 | $\int_{-\infty}^{+\infty}f^2(t)\,\mathrm dt=\frac1{2\pi}\int_{-\infty}^{+\infty}\lvert F(\omega)\rvert^2\,\mathrm d\omega$ |

## 第 9 章 拉普拉斯变换

| 内容 | 公式 |
|---|---|
| 定义 | $F(s)=\mathcal L[f(t)]=\int_0^{+\infty}f(t)e^{-st}\,\mathrm dt$，$s=\beta+j\omega$ |
| 与傅里叶变换的关系 | $\mathcal L[f]=\mathcal F\left[f(t)u(t)e^{-\beta t}\right]$ |
| 存在定理 | $\lvert f(t)\rvert\le Me^{ct}$ ⇒ $F(s)$ 在 $\operatorname{Re}s>c$ 内存在且解析 |
| 常用变换 | $\mathcal L[1]=\frac1s$，$\mathcal L[\delta]=1$，$\mathcal L[e^{at}]=\frac1{s-a}$，$\mathcal L[t^m]=\frac{m!}{s^{m+1}}$ |
| 正弦、余弦 | $\mathcal L[\sin bt]=\frac{b}{s^2+b^2}$，$\mathcal L[\cos bt]=\frac{s}{s^2+b^2}$ |
| 位移 | $\mathcal L[e^{at}f(t)]=F(s-a)$ |
| 延迟 | $\mathcal L[f(t-\tau)u(t-\tau)]=e^{-s\tau}F(s)$（$\tau\ge0$） |
| 相似 | $\mathcal L[f(at)]=\frac1aF\left(\frac sa\right)$（$a>0$） |
| 微分 | $\mathcal L[f^{(n)}]=s^nF(s)-s^{n-1}f(0)-s^{n-2}f'(0)-\cdots-f^{(n-1)}(0)$ |
| 象函数的导数 | $F^{(n)}(s)=(-1)^n\mathcal L[t^nf(t)]$ |
| 积分 | $\mathcal L\left[\int_0^tf\,\mathrm d\tau\right]=\frac{F(s)}s$；$\int_s^{\infty}F(s)\,\mathrm ds=\mathcal L\left[\frac{f(t)}t\right]$ |
| 周期函数 | $\mathcal L[f]=\frac1{1-e^{-sT}}\int_0^Tf(t)e^{-st}\,\mathrm dt$ |
| 卷积定理 | $\mathcal L\left[\int_0^tf_1(\tau)f_2(t-\tau)\,\mathrm d\tau\right]=F_1(s)F_2(s)$ |
| 反演积分 | $f(t)=\frac1{2\pi j}\int_{\beta-j\infty}^{\beta+j\infty}F(s)e^{st}\,\mathrm ds=\sum_k\operatorname{Res}[F(s)e^{st},s_k]$（$t>0$） |
| 解常微分方程 | 两边取 $\mathcal L$（代入初值）→ 解代数方程得 $Y(s)$ → 求 $\mathcal L^{-1}$ |
