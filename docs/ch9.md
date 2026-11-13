# 第 9 章 拉普拉斯变换

!!! tip "本章目标"

    1. 理解拉普拉斯变换是“乘上衰减因子 $e^{-\sigma t}$ 后的傅里叶变换”；掌握定义、**存在定理**，知道象函数在半平面 $\operatorname{Re}s>c$（$c$ 由增长指数决定，可以大于 $0$）内**解析**；熟记常用函数的变换（$1,\ \delta(t),\ t^m,\ e^{at},\ \sin\omega t,\ \cos\omega t$）；
    2. 熟练运用**性质**求变换和逆变换：线性、相似、延迟、位移、微分（原函数的导数、象函数的导数）、积分（原函数的积分、象函数的积分）、初值与终值定理、**卷积定理**——期末必考；
    3. 会用两种方法求**逆变换**：部分分式 + 查表，以及**留数法** $f(t)=\sum\operatorname{Res}\!\left[F(s)e^{st},s_k\right]$，理解反演积分（Bromwich 积分）为什么能化成留数之和；
    4. 会用拉普拉斯变换解**常微分方程（组）的初值问题**、卷积型积分方程和电路方程——期末必考；理解**传递函数** $H(s)$、极点位置与时间响应、稳定性的关系。

**本章主线。** 第 8 章的傅里叶变换要求 $f(t)$ 在整个实轴上绝对可积，$1$、$t$、$\sin t$、$e^{2t}$ 这些最常见的信号都不满足。拉普拉斯变换做了两件事：只看 $t\ge0$（乘 $u(t)$），再乘一个衰减因子 $e^{-\sigma t}$。把 $\sigma+j\omega$ 合写成一个复数 $s$，变换的结果 $F(s)$ 就成了 $s$ 平面右半部分的**解析函数**，前八章的复变工具全部可以用上：

- 求逆变换的反演积分是一个复积分，补上左半圆弧后用**约当引理**和**留数定理**（第 5 章）算出来；
- 求导变成乘 $s$，积分变成除以 $s$，于是**微分方程变成代数方程**；
- 卷积变成乘积，系统的输入输出关系变成 $Y(s)=H(s)X(s)$，$H(s)$ 的**极点**在 $s$ 平面上的位置直接决定系统是衰减、振荡还是发散。

!!! note "记号约定"

    积分变换部分沿用工程中的习惯，用 $j$ 表示虚数单位（$j^2=-1$，与前六章的 $i$ 相同），复变量记为 $s=\sigma+j\omega$。$\mathcal L$ 表示拉普拉斯变换，$\mathcal L^{-1}$ 表示逆变换，$\mathcal F$ 表示傅里叶变换。$u(t)$ 为单位阶跃函数，$\delta(t)$ 为单位冲激函数（第 8 章）。

---

## 9.1 拉普拉斯变换的概念

### 9.1.1 从傅里叶变换到拉普拉斯变换

古典傅里叶变换 $\mathcal F[f]=\int_{-\infty}^{+\infty}f(t)e^{-j\omega t}\,\mathrm dt$ 有两个不方便的地方：

1. **绝对可积是很强的条件。** 常数、$t$、$\sin\omega t$ 都不绝对可积，只能借助 $\delta$ 函数得到广义傅里叶变换；$e^{2t}$ 这种增长的信号连广义变换都没有。
2. **很多实际问题只关心 $t\ge0$。** 电路在 $t=0$ 时接通，之前的情况由初始条件概括，没有必要在整个实轴上积分。

改造的办法很直接：先乘 $u(t)$，把 $t<0$ 的部分置零；再乘衰减因子 $e^{-\sigma t}$（$\sigma>0$ 足够大），把 $t>0$ 的部分压下来。对 $f(t)u(t)e^{-\sigma t}$ 做傅里叶变换：

$$\mathcal F\!\left[f(t)u(t)e^{-\sigma t}\right]=\int_{0}^{+\infty}f(t)e^{-\sigma t}e^{-j\omega t}\,\mathrm dt=\int_{0}^{+\infty}f(t)e^{-(\sigma+j\omega)t}\,\mathrm dt .$$

把 $\sigma+j\omega$ 记成一个复数 $s$，就得到拉普拉斯变换。**拉普拉斯变换就是“乘衰减因子后的傅里叶变换”**，$\sigma$ 管衰减，$\omega$ 管频率。

<figure class="video" markdown>
<video controls preload="metadata" width="100%" poster="assets/posters/C9_DampingROC.jpg">
<source src="assets/videos/C9_DampingROC.mp4" type="video/mp4">
</video>
<figcaption>动画 9-1　$f(t)=e^{0.5t}$ 没有古典傅里叶变换。乘上 $e^{-\sigma t}$ 后：$\sigma\le0.5$ 时仍不衰减，积分发散；$\sigma>0.5$ 时指数衰减，积分收敛。所有使积分收敛的 $s$ 组成右半平面 $\operatorname{Re}s>0.5$，象函数 $F(s)=\frac1{s-0.5}$ 的极点恰好在这个半平面的边界上。</figcaption>
</figure>

### 9.1.2 定义

!!! definition "定义 9.1（拉普拉斯变换）"

    设 $f(t)$ 在 $t\ge0$ 时有定义。若对复参数 $s=\sigma+j\omega$，积分

    $$F(s)=\int_0^{+\infty}f(t)e^{-st}\,\mathrm dt$$

    在 $s$ 平面的某一区域内收敛，则称 $F(s)$ 为 $f(t)$ 的**拉普拉斯变换**（**象函数**），记为 $F(s)=\mathcal L[f(t)]$；称 $f(t)$ 为 $F(s)$ 的**拉普拉斯逆变换**（**象原函数**），记为 $f(t)=\mathcal L^{-1}[F(s)]$。

两点约定，全章通用：

- **默认 $t<0$ 时 $f(t)=0$。** 拉普拉斯变换只用到 $t\ge0$ 的值，所以 $f(t)$ 与 $f(t)u(t)$ 的象函数相同。比如 $\mathcal L[1]=\mathcal L[u(t)]=\frac1s$，逆变换 $\mathcal L^{-1}\!\left[\frac1s\right]$ 写成 $1$ 或 $u(t)$ 都可以，指的都是 $t\ge0$ 时等于 $1$ 的函数。
- **积分下限理解为 $0^-$。** 当 $f(t)$ 在 $t=0$ 处含有冲激 $\delta(t)$ 时，积分从 $0^-$ 开始，把冲激完整地包进来。不含冲激时 $0^-$ 和 $0$ 没有区别。

!!! example "例题 9.1"

    求 $\mathcal L[1]$、$\mathcal L[e^{at}]$（$a$ 为复常数），并写出积分收敛的范围。

    **解。**

    $$\mathcal L[1]=\int_0^{+\infty}e^{-st}\,\mathrm dt=\left.-\frac{e^{-st}}{s}\right|_0^{+\infty}=\frac1s\qquad(\operatorname{Re}s>0),$$

    因为 $|e^{-st}|=e^{-t\operatorname{Re}s}$，只有 $\operatorname{Re}s>0$ 时 $t\to+\infty$ 的极限才为 $0$。同理

    $$\mathcal L[e^{at}]=\int_0^{+\infty}e^{-(s-a)t}\,\mathrm dt=\frac{1}{s-a}\qquad(\operatorname{Re}s>\operatorname{Re}a).$$

    $a=0$ 时就是 $\mathcal L[1]$。注意 $\mathcal L[e^{2t}]=\frac1{s-2}$ 存在，而 $e^{2t}$ 连广义傅里叶变换都没有。

### 9.1.3 存在定理与象函数的解析性

!!! definition "定义 9.2（指数级函数）"

    若存在常数 $M>0$ 和 $c\ge0$，使得 $t\ge0$ 时 $|f(t)|\le Me^{ct}$，则称 $f(t)$ 是**指数级函数**（增长不超过指数级），$c$ 称为它的**增长指数**。

$1$、$t^m$、$\sin\omega t$、$e^{at}$ 都是指数级函数；$e^{t^2}$ 增长得比任何 $e^{ct}$ 都快，不是指数级函数，它的拉普拉斯变换不存在。

!!! theorem "定理 9.1（存在定理）"

    设 $f(t)$ 满足：

    1. 在 $t\ge0$ 的任何有限区间上分段连续；
    2. 是指数级函数，增长指数为 $c$。

    则象函数 $F(s)=\int_0^{+\infty}f(t)e^{-st}\,\mathrm dt$ 在半平面 $\operatorname{Re}s>c$ 上存在，并且在该半平面内**解析**，

    $$F'(s)=\int_0^{+\infty}(-t)f(t)e^{-st}\,\mathrm dt=\mathcal L[-tf(t)] .$$

??? proof "证明要点"

    **收敛。** 设 $\operatorname{Re}s=\sigma>c$，则 $|f(t)e^{-st}|\le Me^{ct}e^{-\sigma t}=Me^{-(\sigma-c)t}$，而 $\int_0^{+\infty}Me^{-(\sigma-c)t}\,\mathrm dt=\frac{M}{\sigma-c}<+\infty$，所以积分绝对收敛。

    **解析。** 在任意闭半平面 $\operatorname{Re}s\ge c_1>c$ 上，$|f(t)e^{-st}|\le Me^{-(c_1-c)t}$，这个控制函数与 $s$ 无关且可积，所以积分**一致收敛**。同样 $|tf(t)e^{-st}|\le Mte^{-(c_1-c)t}$ 也可积，可以在积分号下对 $s$ 求导，得到上面的 $F'(s)$。$F'(s)$ 存在，说明 $F$ 在 $\operatorname{Re}s>c$ 内解析（第 2 章）。$\blacksquare$

这个定理有三个推论，后面反复用到：

- **收敛域是一个右半平面 $\operatorname{Re}s>c$。** 这里的 $c$ 由 $f$ 的增长速度决定，可正可负（如 $e^{2t}$ 的收敛域是 $\operatorname{Re}s>2$，不含虚轴），并不一定是 $\operatorname{Re}s>0$。只要 $\operatorname{Re}s$ 足够大，积分就收敛，所以做变换时通常略去收敛域，只在需要时注明。
- **奇点都在收敛域的左边。** $F(s)$ 在 $\operatorname{Re}s>c$ 内解析，所以它的奇点（对有理函数就是极点）全在某条竖直线 $\operatorname{Re}s=c$ 的左侧。求逆变换时，积分路径就选在这条线右边（9.3 节）。
- **$\lim\limits_{\operatorname{Re}s\to+\infty}F(s)=0$。** 由 $|F(s)|\le\frac{M}{\operatorname{Re}s-c}$ 立即可得。所以 $1$、$s$、$s^2$ 这类多项式**不是**普通函数的象函数，它们对应 $\delta(t)$、$\delta'(t)$ 等广义函数；有理函数的象函数在普通函数范围内一定是**真分式**。

存在定理的条件是**充分**的，不是必要的。$f(t)=\frac{1}{\sqrt t}$ 在 $t=0$ 附近无界，不满足条件 1，但 $\mathcal L\!\left[\frac1{\sqrt t}\right]=\sqrt{\frac{\pi}{s}}$ 照样存在（见例题 9.3）。

### 9.1.4 常用函数的拉普拉斯变换

!!! example "例题 9.2"

    求 $\mathcal L[\delta(t)]$ 和 $\mathcal L[\delta(t-\tau)]$（$\tau>0$）。

    **解。** 由 $\delta$ 函数的筛选性质（积分下限取 $0^-$）：

    $$\mathcal L[\delta(t)]=\int_{0^-}^{+\infty}\delta(t)e^{-st}\,\mathrm dt=e^{-s\cdot0}=1,\qquad \mathcal L[\delta(t-\tau)]=e^{-s\tau}.$$

    对整个 $s$ 平面都成立。象函数 $1$ 不趋于 $0$，再次说明 $\delta(t)$ 不是普通函数。

!!! example "例题 9.3"

    求 $\mathcal L[t^m]$（$m$ 为正整数），以及 $\mathcal L[t^{\alpha}]$（$\alpha>-1$）。

    **解。** 分部积分（$\operatorname{Re}s>0$）：

    $$\mathcal L[t^m]=\int_0^{+\infty}t^me^{-st}\,\mathrm dt=\left.-\frac{t^me^{-st}}{s}\right|_0^{+\infty}+\frac ms\int_0^{+\infty}t^{m-1}e^{-st}\,\mathrm dt=\frac ms\,\mathcal L[t^{m-1}] .$$

    反复递推，$\mathcal L[t^m]=\frac{m}{s}\cdot\frac{m-1}{s}\cdots\frac1s\,\mathcal L[1]=\dfrac{m!}{s^{m+1}}$。

    一般地，先取 $s>0$ 为实数，令 $x=st$：

    $$\mathcal L[t^\alpha]=\int_0^{+\infty}\frac{x^\alpha}{s^\alpha}e^{-x}\,\frac{\mathrm dx}{s}=\frac{\Gamma(\alpha+1)}{s^{\alpha+1}},\qquad \Gamma(x)=\int_0^{+\infty}t^{x-1}e^{-t}\,\mathrm dt .$$

    两边都在 $\operatorname{Re}s>0$ 内解析（右边取 $s^{\alpha+1}$ 的主值），在正实轴上相等，由解析函数的唯一性定理（两个在区域内解析的函数若在一段线段上相等，则在整个区域内相等）在整个右半平面相等。

    $\Gamma$ 函数满足 $\Gamma(1)=1$，$\Gamma(x+1)=x\Gamma(x)$，所以 $\Gamma(m+1)=m!$；另外 $\Gamma\!\left(\frac12\right)=\sqrt\pi$，于是 $\mathcal L\!\left[t^{-1/2}\right]=\sqrt{\frac\pi s}$。

!!! example "例题 9.4"

    求 $\mathcal L[\cos\omega t]$、$\mathcal L[\sin\omega t]$（$\omega$ 为实数）。

    **解。** 用欧拉公式和例题 9.1（$\operatorname{Re}s>0$）：

    $$\mathcal L[\cos\omega t]=\frac12\left(\mathcal L[e^{j\omega t}]+\mathcal L[e^{-j\omega t}]\right)=\frac12\left(\frac1{s-j\omega}+\frac1{s+j\omega}\right)=\frac{s}{s^2+\omega^2},$$

    $$\mathcal L[\sin\omega t]=\frac1{2j}\left(\frac1{s-j\omega}-\frac1{s+j\omega}\right)=\frac{\omega}{s^2+\omega^2}.$$

    同理 $\mathcal L[\cosh at]=\frac{s}{s^2-a^2}$，$\mathcal L[\sinh at]=\frac{a}{s^2-a^2}$（$\operatorname{Re}s>|a|$）。

**规律。** 常用函数的象函数几乎都是 $s$ 的**有理函数**，并且**极点的位置对应时间函数的“指数”**：$e^{at}\leftrightarrow$ 极点 $a$；$\sin\omega t,\cos\omega t\leftrightarrow$ 极点 $\pm j\omega$；$t^m\leftrightarrow$ $0$ 处的 $m+1$ 阶极点。这条对应在 9.4 节会变成判断稳定性的工具。

!!! theorem "定理 9.2（周期函数的象函数）"

    设 $f(t)$ 是 $[0,+\infty)$ 上以 $T$ 为周期的函数，且在一个周期内分段连续，则

    $$\mathcal L[f(t)]=\frac{1}{1-e^{-sT}}\int_0^Tf(t)e^{-st}\,\mathrm dt\qquad(\operatorname{Re}s>0).$$

??? proof "证明"

    把积分按周期拆开，第 $k$ 段令 $t=\tau+kT$：

    $$\int_0^{+\infty}f(t)e^{-st}\,\mathrm dt=\sum_{k=0}^{\infty}\int_{kT}^{(k+1)T}f(t)e^{-st}\,\mathrm dt=\sum_{k=0}^\infty e^{-skT}\int_0^Tf(\tau)e^{-s\tau}\,\mathrm d\tau .$$

    $\operatorname{Re}s>0$ 时 $|e^{-sT}|<1$，等比级数求和得 $\sum_k e^{-skT}=\frac1{1-e^{-sT}}$。$\blacksquare$

!!! example "例题 9.5"

    求周期为 $2a$ 的方波 $f(t)=\begin{cases}1,&0\le t<a,\\-1,&a\le t<2a\end{cases}$ 的象函数。

    **解。**

    $$\int_0^{2a}f(t)e^{-st}\,\mathrm dt=\frac{1-e^{-as}}{s}-\frac{e^{-as}-e^{-2as}}{s}=\frac{(1-e^{-as})^2}{s},$$

    $$\mathcal L[f]=\frac{(1-e^{-as})^2}{s(1-e^{-2as})}=\frac{1-e^{-as}}{s(1+e^{-as})}=\frac1s\tanh\frac{as}{2}.$$

!!! warning "易错点：收敛域不能乱丢"

    “通常略去收敛域”只是书写习惯，不是说收敛域不存在。$\mathcal L[e^{at}]=\frac1{s-a}$ 只在 $\operatorname{Re}s>\operatorname{Re}a$ 时成立。9.2 节用 $F(0)$ 计算 $\int_0^{+\infty}f(t)\,\mathrm dt$、用终值定理求 $f(+\infty)$ 时，都必须先确认 $s=0$ 在收敛域内（或在边界上且极限存在），否则会得出荒谬的结果。

---

## 9.2 拉普拉斯变换的性质

本节总假定所涉及函数的象函数都存在，$F(s)=\mathcal L[f(t)]$，$G(s)=\mathcal L[g(t)]$，并且 $t<0$ 时 $f(t)=g(t)=0$。求导、积分、求极限的次序交换不再一一说明。

**一句话概括。** 每条性质都是“时间域的一种操作 $\leftrightarrow$ $s$ 域的另一种操作”。记住下面这张对照，性质就不用死记：

| 时间域 $t$ | $s$ 域 |
|---|---|
| 平移 $f(t-\tau)$（延迟） | 乘 $e^{-s\tau}$ |
| 乘 $e^{at}$ | 平移 $F(s-a)$ |
| 求导 $f'(t)$ | 乘 $s$（再减初值） |
| 乘 $-t$ | 求导 $F'(s)$ |
| 积分 $\int_0^tf$ | 除以 $s$ |
| 除以 $t$ | 积分 $\int_s^\infty F$ |
| 卷积 $f*g$ | 乘积 $F\cdot G$ |

左右两列是对称的，这和第 8 章傅里叶变换的性质表一模一样，只是把 $j\omega$ 换成了 $s$。

### 9.2.1 线性性质与相似性质

!!! theorem "定理 9.3（线性性质）"

    对任意常数 $\alpha,\beta$，

    $$\mathcal L[\alpha f(t)+\beta g(t)]=\alpha F(s)+\beta G(s),\qquad \mathcal L^{-1}[\alpha F(s)+\beta G(s)]=\alpha f(t)+\beta g(t).$$

（积分的线性性。例题 9.4 已经用过。）

!!! example "例题 9.6"

    求 $\mathcal L^{-1}\!\left[\dfrac{1}{s^2-3s+2}\right]$。

    **解。** $\dfrac1{s^2-3s+2}=\dfrac1{(s-1)(s-2)}=\dfrac1{s-2}-\dfrac1{s-1}$，所以 $\mathcal L^{-1}=e^{2t}-e^{t}$。

!!! theorem "定理 9.4（相似性质 / 尺度变换）"

    设 $a>0$，则 $\mathcal L[f(at)]=\dfrac1aF\!\left(\dfrac sa\right)$。

??? proof "证明"

    令 $x=at$：$\displaystyle\int_0^{+\infty}f(at)e^{-st}\,\mathrm dt=\frac1a\int_0^{+\infty}f(x)e^{-\frac sax}\,\mathrm dx=\frac1aF\!\left(\frac sa\right)$。$\blacksquare$

时间压缩 $a$ 倍，$s$ 域就拉伸 $a$ 倍，这和傅里叶变换的尺度性质是同一个现象。

### 9.2.2 延迟性质与位移性质

!!! theorem "定理 9.5（延迟性质 / 时移性质）"

    设 $\tau\ge0$，$t<0$ 时 $f(t)=0$，则

    $$\mathcal L[f(t-\tau)]=e^{-s\tau}F(s),\qquad\text{即}\qquad \mathcal L[f(t-\tau)u(t-\tau)]=e^{-s\tau}F(s);$$

    反过来，$\mathcal L^{-1}[e^{-s\tau}F(s)]=f(t-\tau)u(t-\tau)$。

??? proof "证明"

    因为 $t<\tau$ 时 $f(t-\tau)=0$，令 $x=t-\tau$：

    $$\int_0^{+\infty}f(t-\tau)e^{-st}\,\mathrm dt=\int_\tau^{+\infty}f(t-\tau)e^{-st}\,\mathrm dt=\int_0^{+\infty}f(x)e^{-s(x+\tau)}\,\mathrm dx=e^{-s\tau}F(s).\ \blacksquare$$

!!! warning "易错点：延迟的是整个函数，不是只延迟 $u$"

    $f(t-\tau)u(t-\tau)$、$f(t)u(t-\tau)$、$f(t-\tau)u(t)$ 是三个不同的函数，只有第一个的象函数是 $e^{-s\tau}F(s)$。例如

    $$\mathcal L[t\,u(t-1)]=\mathcal L[(t-1)u(t-1)+u(t-1)]=e^{-s}\left(\frac1{s^2}+\frac1s\right),$$

    而不是 $\frac{e^{-s}}{s^2}$。做法是先把 $f(t)$ 改写成 $t-\tau$ 的函数，再用延迟性质。求逆变换时同样不能丢掉 $u(t-\tau)$：它表示 $t<\tau$ 时函数为零。

!!! example "例题 9.7"

    (1) 求宽为 $\tau$ 的矩形脉冲 $f(t)=u(t)-u(t-\tau)$ 的象函数；(2) 求 $\mathcal L^{-1}\!\left[\dfrac{e^{-2s}}{s-1}\right]$。

    **解。** (1) $\mathcal L[f]=\dfrac1s-\dfrac{e^{-s\tau}}{s}=\dfrac{1-e^{-s\tau}}{s}$。

    (2) $\mathcal L^{-1}\!\left[\dfrac1{s-1}\right]=e^t$，由延迟性质

    $$\mathcal L^{-1}\!\left[\frac{e^{-2s}}{s-1}\right]=e^{t-2}u(t-2)=\begin{cases}e^{t-2},&t>2,\\0,&t<2.\end{cases}$$

!!! theorem "定理 9.6（位移性质 / 频移性质）"

    $\mathcal L[e^{at}f(t)]=F(s-a)$（$\operatorname{Re}(s-a)>c$）。

（直接代入定义：$\int_0^{+\infty}f(t)e^{-(s-a)t}\,\mathrm dt=F(s-a)$。）

这条性质说明了**衰减振荡**的象函数为什么长成“$(s+\alpha)^2+\omega^2$”的样子：

$$\mathcal L[e^{-\alpha t}\cos\omega t]=\frac{s+\alpha}{(s+\alpha)^2+\omega^2},\qquad \mathcal L[e^{-\alpha t}\sin\omega t]=\frac{\omega}{(s+\alpha)^2+\omega^2}.$$

极点从 $\pm j\omega$ 平移到 $-\alpha\pm j\omega$：乘衰减因子 $e^{-\alpha t}$，极点就向左移动 $\alpha$。

### 9.2.3 微分性质

!!! theorem "定理 9.7（原函数的微分）"

    若 $f'(t)$ 的象函数存在，则

    $$\mathcal L[f'(t)]=sF(s)-f(0).$$

    一般地，

    $$\mathcal L[f^{(n)}(t)]=s^nF(s)-s^{n-1}f(0)-s^{n-2}f'(0)-\cdots-f^{(n-1)}(0),$$

    其中 $f^{(k)}(0)$ 理解为右极限 $\lim\limits_{t\to0^+}f^{(k)}(t)$。

??? proof "证明"

    分部积分：

    $$\int_0^{+\infty}f'(t)e^{-st}\,\mathrm dt=\Big[f(t)e^{-st}\Big]_0^{+\infty}+s\int_0^{+\infty}f(t)e^{-st}\,\mathrm dt .$$

    由 $|f(t)|\le Me^{ct}$，当 $\operatorname{Re}s>c$ 时 $|f(t)e^{-st}|\le Me^{-(\operatorname{Re}s-c)t}\to0$，所以第一项等于 $-f(0)$，得 $\mathcal L[f']=sF(s)-f(0)$。对 $f'$ 再用一次就得到 $\mathcal L[f'']=s\big(sF-f(0)\big)-f'(0)$，依此类推。$\blacksquare$

**这条性质是本章最重要的一条。** 求导变成乘 $s$，初始条件自动出现在公式里。于是线性常系数微分方程的初值问题，一做拉普拉斯变换就变成了关于 $F(s)$ 的**代数方程**（9.4 节）。在电路中，它就是“电感的阻抗是 $sL$，电容的阻抗是 $\frac1{sC}$”的来源。

!!! example "例题 9.8"

    用微分性质重新求 $\mathcal L[t^m]$。

    **解。** 设 $f(t)=t^m$，则 $f(0)=f'(0)=\cdots=f^{(m-1)}(0)=0$，$f^{(m)}(t)=m!$。由定理 9.7，

    $$\mathcal L[m!]=s^m\mathcal L[t^m]\quad\Longrightarrow\quad \mathcal L[t^m]=\frac{m!}{s^m}\cdot\frac1s=\frac{m!}{s^{m+1}}.$$

!!! theorem "定理 9.8（象函数的微分）"

    $F'(s)=\mathcal L[-tf(t)]$；一般地，$F^{(n)}(s)=\mathcal L[(-t)^nf(t)]$。

（就是定理 9.1 中的 $F'(s)$ 公式：在积分号下对 $s$ 求导。）

!!! example "例题 9.9"

    求：(1) $\mathcal L[t\sin\omega t]$；(2) $\mathcal L[te^{-3t}\sin2t]$；(3) $\mathcal L[t^2\cos^2t]$。

    **解。** (1) $\mathcal L[t\sin\omega t]=-\dfrac{\mathrm d}{\mathrm ds}\dfrac{\omega}{s^2+\omega^2}=\dfrac{2\omega s}{(s^2+\omega^2)^2}$。

    (2) 先用位移性质 $\mathcal L[e^{-3t}\sin2t]=\dfrac{2}{(s+3)^2+4}$，再求导：

    $$\mathcal L[te^{-3t}\sin2t]=-\frac{\mathrm d}{\mathrm ds}\frac{2}{(s+3)^2+4}=\frac{4(s+3)}{\left[(s+3)^2+4\right]^2}.$$

    (3) $t^2\cos^2t=\frac12t^2(1+\cos2t)$。$\mathcal L[t^2]=\frac2{s^3}$，$\mathcal L[t^2\cos2t]=\dfrac{\mathrm d^2}{\mathrm ds^2}\dfrac{s}{s^2+4}=\dfrac{2s(s^2-12)}{(s^2+4)^3}$，所以

    $$\mathcal L[t^2\cos^2t]=\frac1{s^3}+\frac{s(s^2-12)}{(s^2+4)^3}=\frac{2(s^6+24s^2+32)}{s^3(s^2+4)^3}.$$

### 9.2.4 积分性质

!!! theorem "定理 9.9（原函数的积分）"

    $$\mathcal L\!\left[\int_0^tf(\tau)\,\mathrm d\tau\right]=\frac{F(s)}{s};\qquad \mathcal L\Big[\underbrace{\int_0^t\mathrm dt\int_0^t\mathrm dt\cdots\int_0^t}_{n\text{ 次}}f(t)\,\mathrm dt\Big]=\frac{F(s)}{s^n}.$$

??? proof "证明"

    令 $g(t)=\int_0^tf(\tau)\,\mathrm d\tau$，则 $g'(t)=f(t)$，$g(0)=0$。由定理 9.7，$F(s)=\mathcal L[g']=sG(s)-0$，所以 $G(s)=\frac{F(s)}{s}$。$\blacksquare$

!!! theorem "定理 9.10（象函数的积分）"

    若 $\int_s^{\infty}F(s)\,\mathrm ds$ 收敛，则

    $$\int_s^{\infty}F(s)\,\mathrm ds=\mathcal L\!\left[\frac{f(t)}{t}\right];\qquad \underbrace{\int_s^\infty\mathrm ds\int_s^\infty\mathrm ds\cdots\int_s^\infty}_{n\text{ 次}}F(s)\,\mathrm ds=\mathcal L\!\left[\frac{f(t)}{t^n}\right].$$

    积分沿 $s$ 向右的水平射线进行（$F$ 解析，积分与路径无关）。

??? proof "证明"

    交换积分次序：$\displaystyle\int_s^\infty\!\!\int_0^{+\infty}f(t)e^{-ut}\,\mathrm dt\,\mathrm du=\int_0^{+\infty}f(t)\left[\int_s^\infty e^{-ut}\,\mathrm du\right]\mathrm dt=\int_0^{+\infty}\frac{f(t)}{t}e^{-st}\,\mathrm dt$。$\blacksquare$

!!! example "例题 9.10"

    求：(1) $\mathcal L\!\left[\displaystyle\int_0^t\tau\sin2\tau\,\mathrm d\tau\right]$；(2) $\mathcal L\!\left[\dfrac{\sin t}{t}\right]$，并由此求 $\displaystyle\int_0^{+\infty}\frac{\sin t}{t}\,\mathrm dt$。

    **解。** (1) 由例题 9.9(1)，$\mathcal L[t\sin2t]=\frac{4s}{(s^2+4)^2}$，再由积分性质

    $$\mathcal L\!\left[\int_0^t\tau\sin2\tau\,\mathrm d\tau\right]=\frac1s\cdot\frac{4s}{(s^2+4)^2}=\frac{4}{(s^2+4)^2}.$$

    (2) $\mathcal L\!\left[\dfrac{\sin t}{t}\right]=\displaystyle\int_s^\infty\frac{\mathrm du}{u^2+1}=\frac\pi2-\arctan s=\operatorname{arccot}s$，即

    $$\int_0^{+\infty}\frac{\sin t}{t}e^{-st}\,\mathrm dt=\operatorname{arccot}s .$$

    令 $s\to0^+$（左边在 $s=0$ 处的连续性可以严格证明），得 $\displaystyle\int_0^{+\infty}\frac{\sin t}{t}\,\mathrm dt=\frac\pi2$。这正是第 5 章用留数和“挖掉原点的半圆”算出的结果，这里用一行就得到了。

### 9.2.5 用拉普拉斯变换计算广义积分

例题 9.10(2) 的思路可以推广：在 $F(s)=\int_0^{+\infty}f(t)e^{-st}\,\mathrm dt$ 及其性质中代入特定的 $s$，就得到一批广义积分：

$$\int_0^{+\infty}f(t)e^{-s_0t}\,\mathrm dt=F(s_0),\qquad\int_0^{+\infty}tf(t)\,\mathrm dt=-F'(0),\qquad\int_0^{+\infty}\frac{f(t)}{t}\,\mathrm dt=\int_0^{\infty}F(s)\,\mathrm ds .$$

使用前要确认代入点在收敛域内，或者确认相应广义积分确实存在。

!!! example "例题 9.11"

    计算：(1) $\displaystyle\int_0^{+\infty}e^{-3t}\cos2t\,\mathrm dt$；(2) $\displaystyle\int_0^{+\infty}\frac{1-\cos t}{t}e^{-t}\,\mathrm dt$。

    **解。** (1) $\mathcal L[\cos2t]=\frac{s}{s^2+4}$，取 $s=3$（在收敛域 $\operatorname{Re}s>0$ 内）：积分 $=\dfrac{3}{9+4}=\dfrac3{13}$。

    (2) $\mathcal L[1-\cos t]=\frac1s-\frac{s}{s^2+1}=\frac{1}{s(s^2+1)}$，由象函数的积分性质

    $$\mathcal L\!\left[\frac{1-\cos t}{t}\right]=\int_s^\infty\left(\frac1u-\frac{u}{u^2+1}\right)\mathrm du=\left[\ln\frac{u}{\sqrt{u^2+1}}\right]_s^\infty=\frac12\ln\frac{s^2+1}{s^2}.$$

    取 $s=1$，积分 $=\dfrac12\ln2$。

### 9.2.6 初值定理与终值定理

!!! theorem "定理 9.11（初值定理与终值定理）"

    设 $f(t)$、$f'(t)$ 的象函数都存在。

    1. **初值定理**：若 $\lim\limits_{s\to\infty}sF(s)$ 存在，则 $f(0^+)=\lim\limits_{s\to\infty}sF(s)$（$s$ 沿正实轴，或在半平面 $\operatorname{Re}s>c$ 内趋于 $\infty$）；
    2. **终值定理**：若 $sF(s)$ 的所有奇点都在左半平面 $\operatorname{Re}s<0$ 内（从而 $f(+\infty)=\lim\limits_{t\to+\infty}f(t)$ 存在），则 $f(+\infty)=\lim\limits_{s\to0}sF(s)$。

??? proof "证明思路"

    由微分性质，$sF(s)-f(0)=\int_0^{+\infty}f'(t)e^{-st}\,\mathrm dt$。

    令 $s\to\infty$：$e^{-st}\to0$（$t>0$），右边 $\to0$，得 $\lim sF(s)=f(0)$。

    令 $s\to0$：右边 $\to\int_0^{+\infty}f'(t)\,\mathrm dt=f(+\infty)-f(0)$，消去 $f(0)$ 得 $\lim\limits_{s\to0}sF(s)=f(+\infty)$。这一步需要 $s=0$ 处能取极限，正是“$sF(s)$ 的奇点全在左半平面”保证的。$\blacksquare$

两个定理的用处是：**不求逆变换就能知道 $f$ 的起点和终点**。比如电路接通瞬间的电流、稳定后的电压，直接从象函数读出。

!!! example "例题 9.12"

    设 $F(s)=\dfrac{1}{s(s+1)}$，不求逆变换，求 $f(0^+)$ 和 $f(+\infty)$，并验证。

    **解。** $sF(s)=\frac1{s+1}$，唯一的极点 $-1$ 在左半平面。

    $$f(0^+)=\lim_{s\to\infty}\frac1{s+1}=0,\qquad f(+\infty)=\lim_{s\to0}\frac1{s+1}=1 .$$

    验证：$F=\frac1s-\frac1{s+1}$，$f(t)=1-e^{-t}$，确实 $f(0)=0$，$f(+\infty)=1$。（这是 RC 电路的充电曲线。）

!!! warning "易错点：终值定理不能乱用"

    $F(s)=\dfrac{1}{s^2+1}$ 时，$\lim\limits_{s\to0}sF(s)=0$，但 $f(t)=\sin t$ 根本没有极限。原因是 $sF(s)$ 的极点 $\pm j$ 在虚轴上。再如 $F(s)=\dfrac1{s-1}$，$\lim\limits_{s\to0}sF(s)=0$，可是 $f(t)=e^t\to+\infty$，极点 $1$ 在右半平面。**用终值定理前必须先检查 $sF(s)$ 的极点是否全在左半平面。**

### 9.2.7 卷积与卷积定理

第 8 章定义的卷积是 $f_1*f_2(t)=\int_{-\infty}^{+\infty}f_1(\tau)f_2(t-\tau)\,\mathrm d\tau$。当 $t<0$ 时 $f_1=f_2=0$，被积函数只在 $0\le\tau\le t$ 时可能非零，于是

!!! definition "定义 9.3（单边卷积）"

    $$f_1(t)*f_2(t)=\int_0^tf_1(\tau)f_2(t-\tau)\,\mathrm d\tau\qquad(t\ge0).$$

它仍满足交换律、结合律和对加法的分配律。

计算卷积的几何过程是：把 $f_2(\tau)$ **翻转**成 $f_2(-\tau)$，**平移** $t$ 得到 $f_2(t-\tau)$，与 $f_1(\tau)$ **相乘**，再**求面积**。

<figure class="video" markdown>
<video controls preload="metadata" width="100%" poster="assets/posters/C9_Convolution.jpg">
<source src="assets/videos/C9_Convolution.mp4" type="video/mp4">
</video>
<figcaption>动画 9-2　$f_1(\tau)=e^{-\tau}$ 与宽 $1.5$ 的矩形脉冲 $f_2$ 的卷积。窗口 $f_2(t-\tau)$ 从左向右滑过，阴影面积随 $t$ 变化就画出 $y(t)=f_1*f_2$：先“充电”到 $1-e^{-1.5}$，脉冲结束后“放电”。象函数就是两个象函数的乘积。</figcaption>
</figure>

!!! theorem "定理 9.12（卷积定理）"

    $$\mathcal L[f_1(t)*f_2(t)]=F_1(s)F_2(s),\qquad \mathcal L^{-1}[F_1(s)F_2(s)]=f_1(t)*f_2(t).$$

??? proof "证明"

    交换积分次序（积分区域 $0\le\tau\le t<+\infty$），再令 $x=t-\tau$：

    $$\int_0^{+\infty}\!\left[\int_0^tf_1(\tau)f_2(t-\tau)\,\mathrm d\tau\right]e^{-st}\,\mathrm dt=\int_0^{+\infty}f_1(\tau)e^{-s\tau}\left[\int_\tau^{+\infty}f_2(t-\tau)e^{-s(t-\tau)}\,\mathrm dt\right]\mathrm d\tau=F_1(s)F_2(s).\ \blacksquare$$

卷积定理在求逆变换时很有用：$F(s)$ 能拆成两个已知逆变换的因子的乘积时，逆变换就是两者的卷积。

!!! example "例题 9.13"

    求 $\mathcal L^{-1}\!\left[\dfrac{s^2}{(s^2+1)^2}\right]$。

    **解。** $\dfrac{s^2}{(s^2+1)^2}=\dfrac{s}{s^2+1}\cdot\dfrac{s}{s^2+1}$，$\mathcal L^{-1}\!\left[\frac{s}{s^2+1}\right]=\cos t$，所以

    $$f(t)=\cos t*\cos t=\int_0^t\cos\tau\cos(t-\tau)\,\mathrm d\tau=\frac12\int_0^t\left[\cos t+\cos(2\tau-t)\right]\mathrm d\tau=\frac12(t\cos t+\sin t).$$

!!! example "例题 9.14"

    用卷积定理求 $\mathcal L^{-1}\!\left[\dfrac{1}{s^2(s^2+1)}\right]$，再用部分分式验证。

    **解。** $\mathcal L^{-1}\!\left[\frac1{s^2}\right]=t$，$\mathcal L^{-1}\!\left[\frac1{s^2+1}\right]=\sin t$，

    $$f(t)=t*\sin t=\int_0^t\tau\sin(t-\tau)\,\mathrm d\tau=\Big[\tau\cos(t-\tau)\Big]_0^t-\int_0^t\cos(t-\tau)\,\mathrm d\tau=t-\sin t .$$

    验证：$\frac1{s^2(s^2+1)}=\frac1{s^2}-\frac1{s^2+1}$，逆变换为 $t-\sin t$，一致。

!!! warning "易错点：乘积的逆变换不是逆变换的乘积"

    $\mathcal L^{-1}[F_1F_2]\ne\mathcal L^{-1}[F_1]\cdot\mathcal L^{-1}[F_2]$。例如 $\mathcal L^{-1}\!\left[\frac1{s}\cdot\frac1s\right]=t$，而 $1\cdot1=1$。$s$ 域的乘积对应时间域的**卷积**；时间域的乘积（如乘 $e^{at}$）对应 $s$ 域的平移。

---

## 9.3 拉普拉斯逆变换

求逆变换有三条路：**查表 + 性质**（9.2 节已经在用）、**部分分式**、**留数**。后两种都针对 $F(s)$ 是有理函数（或含 $e^{-s\tau}$ 因子）的情形，是考试的重点。它们背后的依据是下面的反演积分公式。

### 9.3.1 反演积分公式

由 9.1.1 节，$F(\sigma+j\omega)$ 是 $f(t)u(t)e^{-\sigma t}$ 的傅里叶变换。在 $f$ 的连续点处用傅里叶逆变换（第 8 章）：

$$f(t)u(t)e^{-\sigma t}=\frac1{2\pi}\int_{-\infty}^{+\infty}F(\sigma+j\omega)e^{j\omega t}\,\mathrm d\omega .$$

两边乘 $e^{\sigma t}$，把积分变量换成 $s=\sigma+j\omega$（$\sigma$ 固定，$\mathrm ds=j\,\mathrm d\omega$），就得到：

!!! theorem "定理 9.13（反演积分公式）"

    设 $F(s)=\mathcal L[f(t)]$ 在 $\operatorname{Re}s>c$ 内解析，取 $\beta>c$。则在 $f(t)$ 的连续点 $t>0$ 处

    $$f(t)=\frac1{2\pi j}\int_{\beta-j\infty}^{\beta+j\infty}F(s)e^{st}\,\mathrm ds ;$$

    在间断点处，左边换成 $\frac12\left[f(t^+)+f(t^-)\right]$。积分路径是收敛域内的竖直线 $\operatorname{Re}s=\beta$（称为 **Bromwich 积分**）。

这个公式同时说明了**逆变换的唯一性**：除了间断点上的值以外，$f(t)$ 由 $F(s)$ 唯一确定。所以查表得到的结果就是唯一的答案。

直接沿竖直线积分很难。但被积函数 $F(s)e^{st}$ 是解析函数，可以用第 5 章的老办法：**补一段大圆弧，围成闭路，用留数定理**。

### 9.3.2 留数法

为了让大圆弧上的积分趋于零，需要第 5 章约当引理的一个“旋转”版本。

!!! theorem "引理 9.1（约当引理的另一形式）"

    设 $C_R$ 是以 $\beta$ 为圆心、半径为 $R$ 的**左**半圆弧 $s=\beta+Re^{j\theta}$，$\frac\pi2\le\theta\le\frac{3\pi}2$。若在 $\operatorname{Re}s\le\beta$ 上，当 $|s|\to\infty$ 时 $F(s)$ 一致地趋于 $0$，则对任意 $t>0$，

    $$\lim_{R\to\infty}\int_{C_R}F(s)e^{st}\,\mathrm ds=0 .$$

（令 $s=\beta+jz$，左半平面 $\operatorname{Re}s<\beta$ 变成上半平面 $\operatorname{Im}z>0$，$e^{st}=e^{\beta t}e^{jtz}$，正是第 5 章引理 5.2 中 $e^{iaz}$（$a=t>0$）的形式。）

!!! theorem "定理 9.14（用留数求逆变换）"

    设 $F(s)$ 除有限个孤立奇点 $s_1,s_2,\dots,s_n$ 外处处解析，这些奇点全在直线 $\operatorname{Re}s=\beta$ 的左侧，且当 $s\to\infty$ 时 $F(s)\to0$。则

    $$f(t)=\mathcal L^{-1}[F(s)]=\sum_{k=1}^n\operatorname{Res}\!\left[F(s)e^{st},s_k\right]\qquad(t>0).$$

??? proof "证明"

    取 $R$ 充分大，使所有奇点都落在由线段 $L$（从 $\beta-jR$ 到 $\beta+jR$）和左半圆弧 $C_R$ 围成的闭路内部。由留数定理（第 5 章），

    $$\frac1{2\pi j}\int_LF(s)e^{st}\,\mathrm ds+\frac1{2\pi j}\int_{C_R}F(s)e^{st}\,\mathrm ds=\sum_{k=1}^n\operatorname{Res}\!\left[F(s)e^{st},s_k\right].$$

    令 $R\to\infty$：由引理 9.1，第二项 $\to0$，第一项 $\to$ 反演积分 $f(t)$。$\blacksquare$

    $t<0$ 时，$e^{st}$ 在右半平面很小，应向**右**补圆弧；右边没有奇点，所以积分为 $0$，即 $f(t)=0$。这正好对应“$t<0$ 时 $f=0$”的约定。

<figure class="video" markdown>
<video controls preload="metadata" width="100%" poster="assets/posters/C9_Bromwich.jpg">
<source src="assets/videos/C9_Bromwich.mp4" type="video/mp4">
</video>
<figcaption>动画 9-3　$F(s)=\frac{1}{(s+1)(s+2)}$ 的反演积分。积分路径是所有奇点右侧的竖直线；$t>0$ 时向左补半圆，半径增大到包住两个极点，由约当引理圆弧上的积分趋于零，于是 $f(t)$ 等于两个留数之和 $e^{-t}-e^{-2t}$。$t<0$ 时向右补圆弧，里面没有奇点，$f(t)=0$。</figcaption>
</figure>

**有理函数的情形。** 设 $F(s)=\dfrac{P(s)}{Q(s)}$ 是真分式（$\deg P<\deg Q$，所以 $F\to0$），定理 9.14 的条件自动满足。留数用第 5 章的规则计算：

- $s_k$ 是 $Q$ 的**一阶零点**时：$\operatorname{Res}\!\left[\dfrac{P}{Q}e^{st},s_k\right]=\dfrac{P(s_k)}{Q'(s_k)}e^{s_kt}$；
- $s_k$ 是 $F$ 的 **$m$ 阶极点**时：$\operatorname{Res}\!\left[F(s)e^{st},s_k\right]=\dfrac{1}{(m-1)!}\lim\limits_{s\to s_k}\dfrac{\mathrm d^{m-1}}{\mathrm ds^{m-1}}\left[(s-s_k)^mF(s)e^{st}\right]$。

全是一阶极点时，

$$f(t)=\sum_{k=1}^n\frac{P(s_k)}{Q'(s_k)}e^{s_kt},$$

这在电路分析中叫**赫维赛德（Heaviside）展开定理**。注意求导是对 $s$ 求导，$t$ 当作常数；$m$ 阶极点的留数会带出 $t,t^2,\dots,t^{m-1}$ 这些因子。

### 9.3.3 部分分式法

把真分式 $\frac{P(s)}{Q(s)}$ 拆成简单分式之和，每一项查表。按 $Q(s)$ 的因子分三种情况：

1. **单重一次因子** $(s-a)$：对应项 $\dfrac{A}{s-a}$，$A=\left.\dfrac{P(s)}{Q(s)}(s-a)\right|_{s=a}$（**遮盖法**），逆变换 $Ae^{at}$；
2. **$m$ 重一次因子** $(s-a)^m$：对应项 $\dfrac{A_0}{(s-a)^m}+\dfrac{A_1}{(s-a)^{m-1}}+\cdots+\dfrac{A_{m-1}}{s-a}$，

    $$A_k=\frac1{k!}\left.\frac{\mathrm d^k}{\mathrm ds^k}\left[(s-a)^m\frac{P(s)}{Q(s)}\right]\right|_{s=a},$$

    逆变换用 $\mathcal L^{-1}\!\left[\dfrac1{(s-a)^{k+1}}\right]=\dfrac{t^k}{k!}e^{at}$；
3. **单重二次因子** $(s-a)^2+b^2$（一对共轭复根 $a\pm jb$）：对应项 $\dfrac{C(s-a)+bD}{(s-a)^2+b^2}$，逆变换 $e^{at}(C\cos bt+D\sin bt)$。系数可以比较系数求，也可以在 $s=a+jb$ 处取值：记 $Q(s)=[(s-a)^2+b^2]Q_3(s)$，则 $\dfrac{P(a+jb)}{Q_3(a+jb)}=jbC+bD$。

**两种方法的关系。** 部分分式中 $\frac{A}{s-a}$ 的系数 $A$ 正是 $F$ 在 $a$ 处的留数，$\frac{A}{s-a}$ 的逆变换 $Ae^{at}$ 正是 $F(s)e^{st}$ 在 $a$ 处的留数。所以两种方法本质上是一回事：留数法直接一步到位，部分分式法先拆再查表。复数极点成对出现时，部分分式法可以全程用实数，留数法则用“共轭两项之和 $=2\operatorname{Re}$”。

!!! example "例题 9.15（单极点）"

    求 $F(s)=\dfrac{5s-1}{(s+1)(s-2)}$ 的逆变换。

    **解法一（部分分式）。** 设 $F=\dfrac{A}{s+1}+\dfrac{B}{s-2}$，遮盖法：

    $$A=\left.\frac{5s-1}{s-2}\right|_{s=-1}=\frac{-6}{-3}=2,\qquad B=\left.\frac{5s-1}{s+1}\right|_{s=2}=\frac{9}{3}=3 .$$

    所以 $f(t)=2e^{-t}+3e^{2t}$。

    **解法二（留数）。** $s=-1$、$s=2$ 都是一阶极点，

    $$\operatorname{Res}[Fe^{st},-1]=\left.\frac{(5s-1)e^{st}}{s-2}\right|_{s=-1}=2e^{-t},\qquad \operatorname{Res}[Fe^{st},2]=\left.\frac{(5s-1)e^{st}}{s+1}\right|_{s=2}=3e^{2t}.$$

    $f(t)=2e^{-t}+3e^{2t}$。

!!! example "例题 9.16（重极点）"

    求 $F(s)=\dfrac{1}{(s-2)(s-1)^2}$ 的逆变换。

    **解法一（留数）。** $s=2$ 为一阶极点，$s=1$ 为二阶极点。

    $$\operatorname{Res}[Fe^{st},2]=\left.\frac{e^{st}}{(s-1)^2}\right|_{s=2}=e^{2t},$$

    $$\operatorname{Res}[Fe^{st},1]=\lim_{s\to1}\frac{\mathrm d}{\mathrm ds}\frac{e^{st}}{s-2}=\lim_{s\to1}\frac{te^{st}(s-2)-e^{st}}{(s-2)^2}=-te^t-e^t .$$

    所以 $f(t)=e^{2t}-e^t-te^t$。

    **解法二（部分分式）。** $F=\dfrac{A}{s-2}+\dfrac{A_0}{(s-1)^2}+\dfrac{A_1}{s-1}$，$A=1$，$A_0=\left.\frac1{s-2}\right|_{s=1}=-1$，$A_1=\left.\left(\frac1{s-2}\right)'\right|_{s=1}=\left.-\frac1{(s-2)^2}\right|_{s=1}=-1$。逆变换 $e^{2t}-te^t-e^t$，一致。

!!! example "例题 9.17（共轭复极点）"

    求 $F(s)=\dfrac{s+3}{s^2+2s+5}$ 的逆变换。

    **解法一（配方 + 位移性质）。** $s^2+2s+5=(s+1)^2+4$，

    $$F(s)=\frac{(s+1)+2}{(s+1)^2+2^2}=\frac{s+1}{(s+1)^2+2^2}+\frac{2}{(s+1)^2+2^2},$$

    所以 $f(t)=e^{-t}(\cos2t+\sin2t)=\sqrt2\,e^{-t}\sin\!\left(2t+\frac\pi4\right)$。

    **解法二（留数）。** 极点 $s_{1,2}=-1\pm2j$，都是一阶的。

    $$\operatorname{Res}[Fe^{st},-1+2j]=\left.\frac{(s+3)e^{st}}{2s+2}\right|_{s=-1+2j}=\frac{2+2j}{4j}e^{(-1+2j)t}=\frac{1-j}{2}e^{-t}e^{2jt}.$$

    另一个留数是它的共轭，两者之和是实部的 $2$ 倍：

    $$f(t)=2\operatorname{Re}\!\left[\frac{1-j}{2}e^{-t}(\cos2t+j\sin2t)\right]=e^{-t}(\cos2t+\sin2t).$$

    极点 $-1\pm2j$ 的实部 $-1$ 给出衰减 $e^{-t}$，虚部 $\pm2$ 给出振荡频率 $2$。

!!! example "例题 9.18（非有理象函数）"

    求：(1) $\mathcal L^{-1}\!\left[\dfrac{1-e^{-2s}}{s^2}\right]$；(2) $\mathcal L^{-1}\!\left[\ln\dfrac{s+1}{s-1}\right]$。

    **解。** (1) $\mathcal L^{-1}\!\left[\frac1{s^2}\right]=t$，由线性和延迟性质

    $$f(t)=t-(t-2)u(t-2)=\begin{cases}t,&0\le t<2,\\2,&t\ge2,\end{cases}$$

    一段斜坡加一段平台。$e^{-s\tau}$ 型因子不满足定理 9.14 的条件（左半平面里 $|e^{-2s}|$ 无界），要先用延迟性质处理掉。

    (2) 对数不能直接查表，先求导：$F'(s)=\dfrac1{s+1}-\dfrac1{s-1}$。由象函数的微分性质 $F'(s)=\mathcal L[-tf(t)]$，

    $$-tf(t)=e^{-t}-e^{t}\quad\Longrightarrow\quad f(t)=\frac{e^t-e^{-t}}{t}=\frac{2\sinh t}{t}.$$

!!! warning "易错点：求逆变换的几个常见错误"

    1. **用留数法前没检查真分式。** $\deg P\ge\deg Q$ 时 $F\not\to0$，先做多项式除法，整式部分对应 $\delta$ 及其导数（如 $\frac{s+2}{s+1}=1+\frac1{s+1}\leftrightarrow\delta(t)+e^{-t}$）。
    2. **重极点漏了 $t$ 的因子。** $m$ 阶极点处要对 $(s-s_k)^mF(s)e^{st}$ 求导，$e^{st}$ 对 $s$ 求导会出现 $t$，不能把 $e^{st}$ 当常数。
    3. **复极点只算了一个。** 实系数有理函数的复极点成对出现，两个留数互为共轭，结果是 $2\operatorname{Re}[\cdots]$，最后应得到实函数。如果答案里还有 $j$，一定是算错了。
    4. **含 $e^{-s\tau}$ 时直接套留数公式。** 应先用延迟性质，结果里带上 $u(t-\tau)$。

---

## 9.4 拉普拉斯变换的应用

### 9.4.1 求解常微分方程（组）

线性常系数微分方程的初值问题，用拉普拉斯变换解只需要三步：

```mermaid
flowchart LR
    A["微分方程（组）<br/>+ 初始条件"] -->|"取拉普拉斯变换<br/>求导 → 乘 s"| B["象函数的<br/>代数方程（组）"]
    B -->|"解代数方程"| C["象函数 Y(s)"]
    C -->|"逆变换<br/>部分分式 / 留数"| D["原方程的解 y(t)"]
```

和高等数学里“先求齐次通解、再找特解、最后代初值定常数”相比，这个方法的好处是：**初始条件在第一步就自动代进去了**，得到的直接是满足初值的解；右端是阶跃、冲激、分段函数时也一样好用。

!!! example "例题 9.19"

    求解 $y''-3y'+2y=2e^{-t}$，$y(0)=2$，$y'(0)=-1$。

    **解。** 设 $Y(s)=\mathcal L[y(t)]$。由微分性质，$\mathcal L[y']=sY-2$，$\mathcal L[y'']=s^2Y-2s+1$。方程两边取变换：

    $$s^2Y-2s+1-3(sY-2)+2Y=\frac{2}{s+1},$$

    $$(s^2-3s+2)Y=2s-7+\frac2{s+1}\quad\Longrightarrow\quad Y(s)=\frac{2s^2-5s-5}{(s+1)(s-1)(s-2)} .$$

    三个一阶极点，用留数（或遮盖法）：

    $$\operatorname{Res}_{s=-1}=\frac{2+5-5}{(-2)(-3)}e^{-t}=\frac13e^{-t},\quad \operatorname{Res}_{s=1}=\frac{2-5-5}{2\cdot(-1)}e^{t}=4e^{t},\quad \operatorname{Res}_{s=2}=\frac{8-10-5}{3\cdot1}e^{2t}=-\frac73e^{2t}.$$

    所以 $y(t)=\dfrac13e^{-t}+4e^t-\dfrac73e^{2t}$。检验：$y(0)=\frac13+4-\frac73=2$，$y'(0)=-\frac13+4-\frac{14}3=-1$。

!!! example "例题 9.20（重根与共振）"

    求解：(1) $x'''+3x''+3x'+x=6e^{-t}$，$x(0)=x'(0)=x''(0)=0$；(2) $x''+2x'+2x=2e^{-t}\cos t$，$x(0)=x'(0)=0$。

    **解。** (1) 零初值，取变换得 $(s^3+3s^2+3s+1)X=\frac6{s+1}$，即 $(s+1)^3X=\frac6{s+1}$，

    $$X(s)=\frac{6}{(s+1)^4}=\frac{3!}{(s+1)^4}\quad\Longrightarrow\quad x(t)=t^3e^{-t}.$$

    (2) $(s^2+2s+2)X=\dfrac{2(s+1)}{(s+1)^2+1}$，而 $s^2+2s+2=(s+1)^2+1$，

    $$X(s)=\frac{2(s+1)}{\left[(s+1)^2+1\right]^2}.$$

    由例题 9.9(1)，$\mathcal L[t\sin t]=\frac{2s}{(s^2+1)^2}$，再用位移性质，$x(t)=te^{-t}\sin t$。

    两小题的共同点是：右端的“指数”恰好是方程特征根（$-1$ 和 $-1\pm j$），象函数里极点阶数升高，解中出现因子 $t$。这就是**共振**：外力的频率与系统固有频率一致时，振幅随 $t$ 增长。

!!! example "例题 9.21（方程组）"

    求解 $\begin{cases}x'+x-y=e^t,\\y'+3x-2y=2e^t,\end{cases}$ $x(0)=y(0)=1$。

    **解。** 设 $X=\mathcal L[x]$，$Y=\mathcal L[y]$，两式取变换：

    $$\begin{cases}(s+1)X-Y=1+\dfrac1{s-1}=\dfrac{s}{s-1},\\[2mm]3X+(s-2)Y=1+\dfrac2{s-1}=\dfrac{s+1}{s-1}.\end{cases}$$

    由第一式 $Y=(s+1)X-\frac{s}{s-1}$，代入第二式：

    $$\left[3+(s-2)(s+1)\right]X=\frac{s+1}{s-1}+\frac{s(s-2)}{s-1}=\frac{s^2-s+1}{s-1},$$

    而 $3+(s-2)(s+1)=s^2-s+1$，所以 $X=\frac1{s-1}$，$Y=\frac{s+1}{s-1}-\frac{s}{s-1}=\frac1{s-1}$。故 $x(t)=y(t)=e^t$。

!!! example "例题 9.22（卷积型积分方程）"

    求解 $y(t)=t+\displaystyle\int_0^t\sin(t-\tau)\,y(\tau)\,\mathrm d\tau$。

    **解。** 积分项正是卷积 $\sin t*y(t)$。由卷积定理，

    $$Y=\frac1{s^2}+\frac{1}{s^2+1}Y\quad\Longrightarrow\quad Y\cdot\frac{s^2}{s^2+1}=\frac1{s^2}\quad\Longrightarrow\quad Y=\frac{s^2+1}{s^4}=\frac1{s^2}+\frac1{s^4},$$

    所以 $y(t)=t+\dfrac{t^3}{6}$。

!!! warning "易错点：解微分方程时的几个疏忽"

    1. $\mathcal L[y'']=s^2Y-sy(0)-y'(0)$，常见错误是写成 $s^2Y-y(0)-y'(0)$ 或漏掉某个初值。
    2. 右端函数也要变换，并且要用正确的性质（$e^{-t}\cos t$ 用位移性质，$u(t-a)$ 用延迟性质）。
    3. 解出 $Y(s)$ 后先**约分、检查是否为真分式**，再做部分分式；最后用 $y(0)$、$y'(0)$ 检验一下，能发现大部分计算错误。

### 9.4.2 电路的 $s$ 域分析

把元件的电压电流关系取拉普拉斯变换（零初始状态）：

| 元件 | 时间域 | $s$ 域 | 阻抗 $Z(s)=U(s)/I(s)$ |
|---|---|---|---|
| 电阻 $R$ | $u=Ri$ | $U=RI$ | $R$ |
| 电感 $L$ | $u=L\dfrac{\mathrm di}{\mathrm dt}$ | $U=sLI-Li(0)$ | $sL$ |
| 电容 $C$ | $u=\dfrac1C\displaystyle\int_0^ti\,\mathrm dt+u_C(0)$ | $U=\dfrac{I}{sC}+\dfrac{u_C(0)}{s}$ | $\dfrac1{sC}$ |

零初始状态下，电感、电容和电阻一样，都变成了“阻抗”，基尔霍夫定律照常使用，电路的微分方程就变成了**代数方程**。非零初值则以附加电源 $-Li(0)$、$\frac{u_C(0)}{s}$ 的形式出现。令 $s=j\omega$，$sL\to j\omega L$、$\frac1{sC}\to\frac1{j\omega C}$，就回到正弦稳态电路的相量法：相量法只管稳态，拉普拉斯变换连暂态一起算出来。

!!! example "例题 9.23（RLC 串联电路）"

    $R=2\,\Omega$、$L=1\,\mathrm H$、$C=0.2\,\mathrm F$ 串联，$t=0$ 时接入 $E=10\,\mathrm V$ 的直流电源，电感电流与电容电压初值均为零。求回路电流 $i(t)$。

    **解。** 回路方程 $L\dfrac{\mathrm di}{\mathrm dt}+Ri+\dfrac1C\displaystyle\int_0^ti\,\mathrm d\tau=10u(t)$。取变换（或直接写 $s$ 域阻抗）：

    $$\left(s+2+\frac5s\right)I(s)=\frac{10}{s}\quad\Longrightarrow\quad I(s)=\frac{10}{s^2+2s+5}=5\cdot\frac{2}{(s+1)^2+2^2}.$$

    所以 $i(t)=5e^{-t}\sin2t\ (\mathrm A)$。

    极点 $-1\pm2j$ 在左半平面且是复数：电流是**衰减振荡**（欠阻尼），衰减快慢由实部 $-1$ 决定，振荡角频率由虚部 $2$ 决定。$i(0^+)=0$，与初值定理 $\lim\limits_{s\to\infty}sI(s)=0$ 一致；稳态 $i(+\infty)=0$（电容隔直），与终值定理一致。

### 9.4.3 传递函数、极点与稳定性

考虑一般的线性时不变系统，输入 $x(t)$、输出 $y(t)$ 满足

$$a_ny^{(n)}+\cdots+a_1y'+a_0y=b_mx^{(m)}+\cdots+b_1x'+b_0x .$$

在**零初始状态**下两边取变换，$(a_ns^n+\cdots+a_0)Y(s)=(b_ms^m+\cdots+b_0)X(s)$。

!!! definition "定义 9.4（传递函数）"

    零初始状态下，输出与输入的象函数之比

    $$H(s)=\frac{Y(s)}{X(s)}=\frac{b_ms^m+\cdots+b_1s+b_0}{a_ns^n+\cdots+a_1s+a_0}$$

    称为系统的**传递函数**（系统函数）。$H(s)$ 的零点、极点称为系统的零点、极点。

$H(s)$ 只由系统本身决定，与输入无关。它有两个重要的解释：

1. **冲激响应。** 输入 $x(t)=\delta(t)$ 时 $X(s)=1$，输出 $Y(s)=H(s)$。所以 $h(t)=\mathcal L^{-1}[H(s)]$ 是系统的**冲激响应**；由卷积定理，任意输入下 $y(t)=h(t)*x(t)$。**一个系统由它的冲激响应完全刻画。**
2. **频率响应。** 若 $H(s)$ 的极点都在左半平面，令 $s=j\omega$，$H(j\omega)$ 就是 $h(t)$ 的傅里叶变换（第 8 章），称为**频率响应**：输入 $e^{j\omega t}$ 时，稳态输出是 $H(j\omega)e^{j\omega t}$，$|H(j\omega)|$ 是幅度放大倍数，$\arg H(j\omega)$ 是相移。

!!! example "例题 9.24"

    系统满足 $y''+3y'+2y=x(t)$。求传递函数、冲激响应和单位阶跃响应（零初始状态）。

    **解。** $H(s)=\dfrac{1}{s^2+3s+2}=\dfrac1{(s+1)(s+2)}$，冲激响应

    $$h(t)=\mathcal L^{-1}[H(s)]=e^{-t}-e^{-2t}$$

    （这正是动画 9-3 中的反演积分）。单位阶跃响应：$X=\frac1s$，

    $$Y(s)=\frac{1}{s(s+1)(s+2)}\quad\Longrightarrow\quad y(t)=\frac12-e^{-t}+\frac12e^{-2t}.$$

    由终值定理，$y(+\infty)=\lim\limits_{s\to0}sY(s)=H(0)=\frac12$：阶跃响应的稳态值就是“直流增益” $H(0)$。

**极点决定响应的形状。** 设 $H(s)$ 是真分式，用留数法，$h(t)$ 是各极点处留数之和，每个极点贡献一项：

| 极点位置 | $h(t)$ 中对应的项 | 随 $t\to+\infty$ |
|---|---|---|
| 负实轴 $s=-\alpha$（$\alpha>0$） | $Ae^{-\alpha t}$ | 单调衰减 |
| 左半平面共轭对 $-\alpha\pm j\omega_0$ | $e^{-\alpha t}(C\cos\omega_0t+D\sin\omega_0t)$ | 衰减振荡 |
| 原点 $s=0$（一阶） | 常数 $A$ | 不衰减 |
| 虚轴共轭对 $\pm j\omega_0$（一阶） | $C\cos\omega_0t+D\sin\omega_0t$ | 等幅振荡 |
| 虚轴上的重极点 | $t^k\cos\omega_0t$ 等 | 发散 |
| 右半平面 $\operatorname{Re}s>0$ | $e^{\alpha t}(\cdots)$，$\alpha>0$ | 指数发散 |

极点离虚轴越远，衰减越快；离实轴越远，振荡越快。

<figure class="video" markdown>
<video controls preload="metadata" width="100%" poster="assets/posters/C9_PoleResponse.jpg">
<source src="assets/videos/C9_PoleResponse.mp4" type="video/mp4">
</video>
<figcaption>动画 9-4　一对极点 $\sigma_0\pm j\omega_0$ 在 $s$ 平面上移动，右边同步画出对应的时间响应 $e^{\sigma_0t}\cos\omega_0t$：负实轴上单调衰减，左半平面衰减振荡，虚轴上等幅振荡，右半平面振幅指数增长。</figcaption>
</figure>

!!! theorem "定理 9.15（稳定性的极点判据）"

    设系统的传递函数 $H(s)$ 是有理真分式（分子分母无公因子）。则系统稳定（冲激响应 $h(t)\to0$，且 $\int_0^{+\infty}|h(t)|\,\mathrm dt<+\infty$，从而有界输入产生有界输出）的充要条件是：**$H(s)$ 的所有极点都在左半平面 $\operatorname{Re}s<0$ 内。**

（由上表：左半平面的极点贡献衰减指数项，可积；只要有一个极点在虚轴或右半平面，对应的项就不衰减。）

!!! example "例题 9.25"

    判断稳定性：(1) $H_1(s)=\dfrac{1}{s^2+s-2}$；(2) $H_2(s)=\dfrac{s}{s^2+4}$；(3) $H_3(s)=\dfrac{s+3}{s^2+2s+5}$。

    **解。** (1) $s^2+s-2=(s-1)(s+2)$，极点 $1$ 在右半平面，$h_1$ 含 $\frac13e^{t}$，**不稳定**。

    (2) 极点 $\pm2j$ 在虚轴上，$h_2(t)=\cos2t$ 等幅振荡，**不稳定**（临界）：输入 $\sin2t$ 时会发生共振，输出振幅随 $t$ 线性增长。

    (3) 极点 $-1\pm2j$ 在左半平面，$h_3(t)=e^{-t}(\cos2t+\sin2t)$（例题 9.17），**稳定**。

**曲面 $|H(s)|$ 与频率响应。** 把 $|H(s)|$ 画成 $s$ 平面上方的曲面，极点处曲面隆起到无穷高（像撑起帐篷的杆子），零点处曲面落到零。沿虚轴 $s=j\omega$ 切一刀，截线就是幅频特性 $|H(j\omega)|$。极点越靠近虚轴，截线在 $\omega$ 接近极点虚部时越“尖”，这就是谐振。

<figure class="video" markdown>
<video controls preload="metadata" width="100%" poster="assets/posters/C9_SurfaceFreqResponse.jpg">
<source src="assets/videos/C9_SurfaceFreqResponse.mp4" type="video/mp4">
</video>
<figcaption>动画 9-5　$H(s)=\frac{1}{(s+0.5)^2+4}$ 的模 $|H(s)|$ 曲面，在极点 $-0.5\pm2j$ 处隆起。沿虚轴切开得到的红色曲线就是频率响应 $|H(j\omega)|$，在 $\omega\approx2$ 处有谐振峰。拉普拉斯变换看整个 $s$ 平面，傅里叶变换只看虚轴这一条线。</figcaption>
</figure>

!!! warning "易错点：稳定性判据"

    1. 判断的是**分母**的根（极点），不是分子的根（零点）；并且要先约去公因子。
    2. 虚轴上的一阶极点**不算稳定**（冲激响应不衰减）；虚轴上的重极点和右半平面极点都会使响应发散。
    3. 传递函数是在**零初始状态**下定义的；有非零初值时，$Y(s)=H(s)X(s)+$（初值项），不能直接用 $H(s)X(s)$。

---

## 知识点归纳总结

!!! abstract "知识点归纳总结"

    ```mermaid
    flowchart TD
        FT["第8章：傅里叶变换<br/>要求绝对可积"] --> DEF["乘 u(t) 和 exp(−σt)<br/>F(s)=∫(0→∞) f(t) exp(−st) dt, s=σ+jω"]
        DEF --> EX["存在定理：指数级 + 分段连续<br/>⇒ Re s > c 内存在且解析"]
        EX --> TAB["常用变换对<br/>1, δ, tᵐ, exp(at), sin, cos"]
        DEF --> PROP["性质<br/>线性 · 相似 · 延迟 · 位移<br/>微分 · 积分 · 初值终值"]
        PROP --> CONV["卷积定理<br/>L[f₁*f₂]=F₁F₂"]
        TAB --> INV["逆变换"]
        PROP --> INV
        INV --> PF["部分分式 + 查表"]
        INV --> RES["反演积分 + 约当引理<br/>f(t)=Σ Res[F(s) exp(st), sₖ]"]
        R5["第5章：留数定理<br/>约当引理"] --> RES
        PROP --> ODE["解常微分方程（组）<br/>求导 → 乘 s"]
        CONV --> IE["卷积型积分方程"]
        ODE --> CIR["电路 s 域分析<br/>Z = R, sL, 1/(sC)"]
        CONV --> H["传递函数 H(s)=Y/X<br/>h(t)=L⁻¹[H], y=h*x"]
        H --> STAB["极点在左半平面 ⇔ 稳定<br/>H(jω) = 频率响应"]
    ```

    **常用变换对**（$\operatorname{Re}s$ 足够大）

    | $f(t)$（$t\ge0$） | $F(s)$ | $f(t)$ | $F(s)$ |
    |---|---|---|---|
    | $\delta(t)$ | $1$ | $\delta(t-\tau)$ | $e^{-s\tau}$ |
    | $1$ 或 $u(t)$ | $\dfrac1s$ | $u(t-\tau)$ | $\dfrac{e^{-s\tau}}{s}$ |
    | $t^m$ | $\dfrac{m!}{s^{m+1}}$ | $t^\alpha$（$\alpha>-1$） | $\dfrac{\Gamma(\alpha+1)}{s^{\alpha+1}}$ |
    | $e^{at}$ | $\dfrac1{s-a}$ | $t^me^{at}$ | $\dfrac{m!}{(s-a)^{m+1}}$ |
    | $\sin\omega t$ | $\dfrac{\omega}{s^2+\omega^2}$ | $\cos\omega t$ | $\dfrac{s}{s^2+\omega^2}$ |
    | $e^{-\alpha t}\sin\omega t$ | $\dfrac{\omega}{(s+\alpha)^2+\omega^2}$ | $e^{-\alpha t}\cos\omega t$ | $\dfrac{s+\alpha}{(s+\alpha)^2+\omega^2}$ |
    | $\sinh at$ | $\dfrac{a}{s^2-a^2}$ | $\cosh at$ | $\dfrac{s}{s^2-a^2}$ |
    | $t\sin\omega t$ | $\dfrac{2\omega s}{(s^2+\omega^2)^2}$ | $t\cos\omega t$ | $\dfrac{s^2-\omega^2}{(s^2+\omega^2)^2}$ |
    | $\dfrac{\sin t}{t}$ | $\operatorname{arccot}s$ | 周期 $T$ 的 $f$ | $\dfrac{\int_0^Tf e^{-st}\,\mathrm dt}{1-e^{-sT}}$ |

    **性质表**（$F=\mathcal L[f]$，$G=\mathcal L[g]$）

    | 性质 | 公式 |
    |---|---|
    | 线性 | $\mathcal L[\alpha f+\beta g]=\alpha F+\beta G$ |
    | 相似（$a>0$） | $\mathcal L[f(at)]=\frac1aF\!\left(\frac sa\right)$ |
    | 延迟（$\tau\ge0$） | $\mathcal L[f(t-\tau)u(t-\tau)]=e^{-s\tau}F(s)$ |
    | 位移 | $\mathcal L[e^{at}f(t)]=F(s-a)$ |
    | 原函数微分 | $\mathcal L[f']=sF-f(0)$，$\mathcal L[f'']=s^2F-sf(0)-f'(0)$ |
    | 象函数微分 | $F^{(n)}(s)=\mathcal L[(-t)^nf(t)]$ |
    | 原函数积分 | $\mathcal L\!\left[\int_0^tf\,\mathrm d\tau\right]=\frac{F(s)}{s}$ |
    | 象函数积分 | $\int_s^\infty F\,\mathrm ds=\mathcal L\!\left[\frac{f(t)}{t}\right]$ |
    | 初值定理 | $f(0^+)=\lim\limits_{s\to\infty}sF(s)$ |
    | 终值定理 | $f(+\infty)=\lim\limits_{s\to0}sF(s)$（$sF$ 的极点全在左半平面） |
    | 卷积定理 | $\mathcal L[f*g]=F\cdot G$，$f*g=\int_0^tf(\tau)g(t-\tau)\,\mathrm d\tau$ |
    | 反演积分 | $f(t)=\frac1{2\pi j}\int_{\beta-j\infty}^{\beta+j\infty}F(s)e^{st}\,\mathrm ds=\sum\operatorname{Res}[F(s)e^{st},s_k]$ |

    **方法清单**

    1. **求象函数**：先查表，再看是哪种“操作”——乘 $e^{at}$ 用位移，乘 $t$ 用象函数微分，除以 $t$ 用象函数积分，延迟用 $e^{-s\tau}$，积分上限是 $t$ 用除以 $s$；周期函数用定理 9.2；
    2. **求逆变换**：有理真分式 → 留数法或部分分式（复极点用配方 + 位移性质）；非真分式 → 先做除法；含 $e^{-s\tau}$ → 延迟性质；含 $\ln$、$\arctan$ → 先对 $s$ 求导；能拆成两个熟悉因子之积 → 卷积；
    3. **解微分方程**：取变换（带初值） → 解出 $Y(s)$ → 逆变换 → 用初值检验；
    4. **积分方程**：认出卷积 $\int_0^tk(t-\tau)y(\tau)\,\mathrm d\tau$，用卷积定理化成代数方程；
    5. **算广义积分**：$\int_0^{+\infty}f(t)e^{-s_0t}\,\mathrm dt=F(s_0)$，$\int_0^{+\infty}\frac{f}{t}\,\mathrm dt=\int_0^\infty F\,\mathrm ds$；
    6. **分析系统**：写 $H(s)$ → 看极点位置判稳定性和响应形状 → 冲激响应 $h=\mathcal L^{-1}[H]$，阶跃响应 $\mathcal L^{-1}[H/s]$，稳态值用终值定理。

    **易错点汇总**

    - 把收敛域丢了，在收敛域外代值，或终值定理用于虚轴、右半平面有极点的情形；
    - 延迟性质只延迟了 $u$ 没延迟 $f$：$\mathcal L[tu(t-1)]\ne\frac{e^{-s}}{s^2}$；
    - $\mathcal L[y'']$ 的初值项写错；
    - 以为 $\mathcal L^{-1}[F_1F_2]=\mathcal L^{-1}[F_1]\cdot\mathcal L^{-1}[F_2]$；
    - 重极点的留数忘了对 $e^{st}$ 求导（漏掉 $t$ 因子），复极点只算了一个；
    - 稳定性看成了零点，或把虚轴上的极点当成稳定。

---

## 练习题

!!! example "练习 1"

    求象函数：(1) $e^{-2t}\sin3t$；(2) $t^2e^{t}$；(3) $t\,u(t-1)$；(4) $\sin\!\left(t-\frac\pi3\right)u\!\left(t-\frac\pi3\right)$。

??? success "解答"

    (1) 位移性质：$\dfrac{3}{(s+2)^2+9}$。

    (2) $\mathcal L[t^2]=\frac2{s^3}$，位移：$\dfrac{2}{(s-1)^3}$。

    (3) $t=(t-1)+1$，$\mathcal L[(t-1)u(t-1)+u(t-1)]=e^{-s}\left(\dfrac1{s^2}+\dfrac1s\right)$。

    (4) 已经是 $f\!\left(t-\frac\pi3\right)u\!\left(t-\frac\pi3\right)$ 的形式，$f=\sin t$：$\dfrac{e^{-\pi s/3}}{s^2+1}$。

!!! example "练习 2"

    求 $\mathcal L\!\left[\dfrac{e^{-t}-e^{-2t}}{t}\right]$，并计算 $\displaystyle\int_0^{+\infty}\frac{e^{-t}-e^{-2t}}{t}\,\mathrm dt$。

??? success "解答"

    $\mathcal L[e^{-t}-e^{-2t}]=\frac1{s+1}-\frac1{s+2}$，由象函数的积分性质

    $$\mathcal L\!\left[\frac{e^{-t}-e^{-2t}}{t}\right]=\int_s^\infty\left(\frac1{u+1}-\frac1{u+2}\right)\mathrm du=\left[\ln\frac{u+1}{u+2}\right]_s^\infty=\ln\frac{s+2}{s+1}.$$

    收敛域为 $\operatorname{Re}s>-1$，包含 $s=0$。令 $s=0$，积分 $=\ln2$。

!!! example "练习 3"

    求 $\mathcal L\!\left[\displaystyle\int_0^t\tau e^{-3\tau}\sin2\tau\,\mathrm d\tau\right]$。

??? success "解答"

    由例题 9.9(2)，$\mathcal L[te^{-3t}\sin2t]=\frac{4(s+3)}{[(s+3)^2+4]^2}$。由积分性质，

    $$\mathcal L\!\left[\int_0^t\tau e^{-3\tau}\sin2\tau\,\mathrm d\tau\right]=\frac{4(s+3)}{s\left[(s+3)^2+4\right]^2}.$$

!!! example "练习 4"

    用部分分式求 $\mathcal L^{-1}\!\left[\dfrac{s+1}{s(s^2+4)}\right]$。

??? success "解答"

    设 $\dfrac{s+1}{s(s^2+4)}=\dfrac As+\dfrac{Bs+C}{s^2+4}$。遮盖法 $A=\frac14$；比较 $s^2$ 系数 $A+B=0$，$B=-\frac14$；比较 $s$ 系数 $C=1$。所以

    $$f(t)=\frac14-\frac14\cos2t+\frac12\sin2t .$$

!!! example "练习 5"

    用留数法求 $\mathcal L^{-1}\!\left[\dfrac{1}{s(s+1)^3}\right]$。

??? success "解答"

    $s=0$ 为一阶极点：$\operatorname{Res}=\left.\frac{e^{st}}{(s+1)^3}\right|_{s=0}=1$。

    $s=-1$ 为三阶极点：

    $$\operatorname{Res}=\frac1{2!}\lim_{s\to-1}\frac{\mathrm d^2}{\mathrm ds^2}\frac{e^{st}}{s}=\frac12\left.e^{st}\left(\frac{t^2}{s}-\frac{2t}{s^2}+\frac{2}{s^3}\right)\right|_{s=-1}=-e^{-t}\left(\frac{t^2}2+t+1\right).$$

    所以 $f(t)=1-e^{-t}\left(1+t+\dfrac{t^2}{2}\right)$。检验：$f(0)=0$，与初值定理 $\lim\limits_{s\to\infty}\frac{1}{(s+1)^3}=0$ 一致。

!!! example "练习 6"

    求解 $y''+4y=3\sin t$，$y(0)=y'(0)=0$。

??? success "解答"

    $(s^2+4)Y=\dfrac{3}{s^2+1}$，

    $$Y=\frac{3}{(s^2+1)(s^2+4)}=\frac1{s^2+1}-\frac1{s^2+4},$$

    所以 $y(t)=\sin t-\dfrac12\sin2t$。检验：$y'(0)=1-1=0$。

!!! example "练习 7"

    (1) 设 $F(s)=\dfrac{s+3}{s(s+1)(s+2)}$，不求逆变换，求 $f(0^+)$ 和 $f(+\infty)$，再求出 $f(t)$ 验证。(2) 对 $F(s)=\dfrac{s}{s^2-1}$ 能否用终值定理？

??? success "解答"

    (1) $sF(s)=\frac{s+3}{(s+1)(s+2)}$，极点 $-1,-2$ 都在左半平面，可以用终值定理。

    $$f(0^+)=\lim_{s\to\infty}\frac{s+3}{(s+1)(s+2)}=0,\qquad f(+\infty)=\lim_{s\to0}\frac{s+3}{(s+1)(s+2)}=\frac32 .$$

    留数法：$s=0$ 处 $\frac32$；$s=-1$ 处 $\frac{2}{(-1)(1)}e^{-t}=-2e^{-t}$；$s=-2$ 处 $\frac{1}{(-2)(-1)}e^{-2t}=\frac12e^{-2t}$。$f(t)=\frac32-2e^{-t}+\frac12e^{-2t}$，$f(0)=0$，$f(+\infty)=\frac32$，一致。

    (2) 不能。$sF(s)=\frac{s^2}{s^2-1}$ 有极点 $s=1$ 在右半平面。形式上 $\lim\limits_{s\to0}sF(s)=0$，但 $f(t)=\cosh t\to+\infty$。

!!! example "练习 8"

    系统满足 $y''+5y'+6y=x'+x$。(1) 求传递函数 $H(s)$ 及冲激响应 $h(t)$，判断稳定性；(2) 求零初始状态下的单位阶跃响应。

??? success "解答"

    (1) 零初值取变换：$(s^2+5s+6)Y=(s+1)X$，

    $$H(s)=\frac{s+1}{(s+2)(s+3)} .$$

    留数：$s=-2$ 处 $\frac{-1}{1}e^{-2t}$；$s=-3$ 处 $\frac{-2}{-1}e^{-3t}$。所以 $h(t)=2e^{-3t}-e^{-2t}$。极点 $-2,-3$ 都在左半平面，系统**稳定**。

    (2) $Y=\dfrac{s+1}{s(s+2)(s+3)}$。留数：$s=0$ 处 $\frac16$；$s=-2$ 处 $\frac{-1}{(-2)(1)}e^{-2t}=\frac12e^{-2t}$；$s=-3$ 处 $\frac{-2}{(-3)(-1)}e^{-3t}=-\frac23e^{-3t}$。

    $$y(t)=\frac16+\frac12e^{-2t}-\frac23e^{-3t},$$

    $y(0)=0$，稳态值 $\frac16=H(0)$，与终值定理一致。

---

**结语。** 到这里，整门课程连成了一条线：第 1–2 章的复数与解析函数提供了 $e^{st}$ 和“解析”这个核心概念；第 3–5 章的复积分、级数和留数，在本章变成了求逆变换的方法；第 8–9 章的积分变换把微分方程、电路和信号问题搬到复平面上解决。全书的知识地图和公式汇总见 [全课程总结](summary.md) 与 [公式速查](formulas.md)。
