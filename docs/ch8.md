# 第 8 章 傅里叶变换

!!! tip "本章目标"

    1. 会把周期函数展开成**傅里叶级数**的三角形式与复指数形式，理解离散频谱 $c_n=F(n\omega_0)$ 的含义；
    2. 理解从傅里叶级数到**傅里叶积分**的过渡（周期 $T\to+\infty$），掌握傅里叶变换 $F(\omega)=\int_{-\infty}^{+\infty}f(t)e^{-j\omega t}\,\mathrm dt$ 与逆变换的定义，会求振幅谱和相位谱；
    3. 理解**单位冲激函数** $\delta(t)$ 的定义与三条性质（筛选、偶函数、$\int\delta=u$），会用它求常数、阶跃、正弦、余弦等函数的**广义傅里叶变换**；
    4. 熟练运用傅里叶变换的**性质**（线性、位移、相似、微分、积分、帕塞瓦尔等式），理解**卷积**的计算与**卷积定理**；
    5. 会用第 3、5 章的工具（柯西定理、留数、约当引理）计算傅里叶变换和逆变换——期末考查“傅里叶变换的概念、$\delta$ 函数、卷积”。

**本章主线。** 第 2 章说 $e^{j\omega t}=\cos\omega t+j\sin\omega t$ 是一个旋转的单位向量。傅里叶的想法是：**任何（足够好的）信号都可以拆成许多个这样的旋转向量之和**。周期信号只需要频率为 $n\omega_0$ 的那些（傅里叶级数，离散谱）；非周期信号可以看成周期为 $\infty$ 的周期信号，需要所有频率（傅里叶积分，连续谱）。变换之后，时域里的求导变成乘 $j\omega$，时域里难算的卷积变成乘积，所以微分方程和线性系统的问题在频域里变成代数问题。

!!! info "记号约定"

    积分变换两章沿用教材和后续信号、电路课程的习惯：**虚数单位记作 $j$**（$j^2=-1$，与前六章的 $i$ 是同一个数），时间变量为 $t$，角频率为 $\omega$。傅里叶变换记作 $\mathcal F$，单位阶跃函数记作 $u(t)$。

---

## 8.1 傅里叶变换的概念

### 8.1.1 傅里叶级数

!!! theorem "定理 8.1（狄利克雷定理）"

    设 $f_T(t)$ 是以 $T$ 为周期的实值函数，且在 $\left[-\frac T2,\frac T2\right]$ 上满足**狄利克雷条件**：

    1. 连续或只有有限个第一类间断点；
    2. 只有有限个极值点。

    令 $\omega_0=\dfrac{2\pi}{T}$（**基频**），则在 $f_T(t)$ 的连续点处有

    $$f_T(t)=\frac{a_0}{2}+\sum_{n=1}^{\infty}\left(a_n\cos n\omega_0t+b_n\sin n\omega_0t\right),\tag{8.1}$$

    $$a_n=\frac2T\int_{-T/2}^{T/2}f_T(t)\cos n\omega_0t\,\mathrm dt\ (n=0,1,2,\dots),\qquad b_n=\frac2T\int_{-T/2}^{T/2}f_T(t)\sin n\omega_0t\,\mathrm dt\ (n=1,2,\dots).$$

    在间断点处，(8.1) 的左端应换成 $\frac12\left[f_T(t+0)+f_T(t-0)\right]$。

(8.1) 称为傅里叶级数的**三角形式**。系数公式来自三角函数系 $\{1,\cos n\omega_0t,\sin n\omega_0t\}$ 在一个周期上的正交性。

**复指数形式。** 由欧拉公式 $\cos n\omega_0t=\frac{e^{jn\omega_0t}+e^{-jn\omega_0t}}{2}$，$\sin n\omega_0t=\frac{e^{jn\omega_0t}-e^{-jn\omega_0t}}{2j}$，代入 (8.1) 整理得

$$f_T(t)=\frac{a_0}{2}+\sum_{n=1}^{\infty}\left(\frac{a_n-jb_n}{2}e^{jn\omega_0t}+\frac{a_n+jb_n}{2}e^{-jn\omega_0t}\right).$$

令 $c_0=\frac{a_0}{2}$，$c_n=\frac{a_n-jb_n}{2}$，$c_{-n}=\frac{a_n+jb_n}{2}=\overline{c_n}$（$n\ge1$），就得到更简洁的形式：

!!! success "傅里叶级数的复指数形式"

    $$f_T(t)=\sum_{n=-\infty}^{+\infty}c_n\,e^{jn\omega_0t},\qquad c_n=\frac1T\int_{-T/2}^{T/2}f_T(t)\,e^{-jn\omega_0t}\,\mathrm dt\quad(n=0,\pm1,\pm2,\dots).\tag{8.2}$$

    $c_n$ 也记作 $F(n\omega_0)$，称为 $f_T$ 的**离散频谱**；$|c_n|$ 为**振幅谱**，$\arg c_n$ 为**相位谱**。

**物理意义。** 记 $A_0=c_0$，$A_n=\sqrt{a_n^2+b_n^2}=2|c_n|$，则 (8.1) 可写为 $f_T(t)=A_0+\sum_{n\ge1}A_n\cos(n\omega_0t+\theta_n)$，其中 $\theta_n=\arg c_n$。若把 $f_T$ 看成一个信号，它就是频率为 $0,\omega_0,2\omega_0,\dots$ 的简谐波之和，第 $n$ 次谐波的振幅为 $A_n$、初相为 $\theta_n$。复指数形式里同一频率被拆成了 $e^{jn\omega_0t}$ 和 $e^{-jn\omega_0t}$ 两个反向旋转的向量，各占一半振幅 $|c_n|=|c_{-n}|=\frac{A_n}{2}$，合起来是实信号。

!!! note "前后联系"

    $\sum c_ne^{jn\omega_0t}$ 令 $z=e^{j\omega_0t}$ 就是 $\sum c_nz^n$，一个在单位圆上取值的**洛朗级数**（第 4 章）；系数公式 $c_n=\frac1T\int f_Te^{-jn\omega_0t}\,\mathrm dt$ 就是洛朗系数公式 $c_n=\frac{1}{2\pi j}\oint\frac{f(z)}{z^{n+1}}\,\mathrm dz$ 在 $|z|=1$ 上的形式。

!!! example "例题 8.1（教材例 8.1）"

    求以 $T$ 为周期的函数 $f_T(t)=\begin{cases}0,&-\frac T2<t<0,\\2,&0<t<\frac T2\end{cases}$ 的离散频谱和傅里叶级数的复指数形式。

    **解。** $\omega_0=\frac{2\pi}{T}$。$n=0$ 时，$c_0=\frac1T\int_0^{T/2}2\,\mathrm dt=1$。$n\ne0$ 时，

    $$c_n=\frac1T\int_0^{T/2}2e^{-jn\omega_0t}\,\mathrm dt=\frac{2}{T}\cdot\frac{1-e^{-jn\omega_0T/2}}{jn\omega_0}=\frac{1-e^{-jn\pi}}{jn\pi}=\begin{cases}0,&n\ \text{为非零偶数},\\[2pt]\dfrac{2}{jn\pi}=-\dfrac{2j}{n\pi},&n\ \text{为奇数}.\end{cases}$$

    所以

    $$f_T(t)=1+\sum_{n=-\infty}^{+\infty}\frac{-2j}{(2n-1)\pi}\,e^{j(2n-1)\omega_0t}.$$

    振幅谱：$|c_0|=1$，$n$ 为奇数时 $|c_n|=\frac{2}{|n|\pi}$，$n$ 为非零偶数时为 $0$。相位谱：$n=1,3,5,\dots$ 时 $\arg c_n=-\frac\pi2$；$n=-1,-3,\dots$ 时 $\arg c_n=\frac\pi2$。

    （合并 $\pm n$ 两项可还原三角形式：$f_T(t)=1+\frac4\pi\left(\sin\omega_0t+\frac13\sin3\omega_0t+\frac15\sin5\omega_0t+\cdots\right)$。）

### 8.1.2 傅里叶积分与傅里叶变换

**非周期函数 = 周期为无穷大的周期函数。** 设 $f(t)$ 定义在整个实轴上，不是周期函数。截取 $\left[-\frac T2,\frac T2\right]$ 上的一段向两边周期延拓，得到周期函数 $f_T(t)$，则 $f(t)=\lim\limits_{T\to+\infty}f_T(t)$。由 (8.2)，

$$f_T(t)=\sum_{n=-\infty}^{+\infty}\left[\frac1T\int_{-T/2}^{T/2}f_T(\tau)e^{-j\omega_n\tau}\,\mathrm d\tau\right]e^{j\omega_nt},\qquad\omega_n=n\omega_0 .$$

记谱线间隔 $\Delta\omega=\omega_0=\frac{2\pi}{T}$，即 $\frac1T=\frac{\Delta\omega}{2\pi}$，则

$$f_T(t)=\frac{1}{2\pi}\sum_{n=-\infty}^{+\infty}\left[\int_{-T/2}^{T/2}f_T(\tau)e^{-j\omega_n\tau}\,\mathrm d\tau\right]e^{j\omega_nt}\,\Delta\omega .$$

$T\to+\infty$ 时 $\Delta\omega\to0$：一方面谱线越来越密，**离散频谱变成连续频谱**；另一方面方括号里的积分趋于 $\int_{-\infty}^{+\infty}f(\tau)e^{-j\omega_n\tau}\,\mathrm d\tau$，求和 $\sum(\cdots)\Delta\omega$ 正是黎曼和，**级数求和变成积分**：

$$f(t)=\frac{1}{2\pi}\int_{-\infty}^{+\infty}\left[\int_{-\infty}^{+\infty}f(\tau)e^{-j\omega\tau}\,\mathrm d\tau\right]e^{j\omega t}\,\mathrm d\omega .\tag{8.3}$$

<figure class="video" markdown>
<video controls preload="metadata" width="100%" poster="assets/posters/C8_SeriesToIntegral.jpg">
<source src="assets/videos/C8_SeriesToIntegral.mp4" type="video/mp4">
</video>
<figcaption>动画 8-1　宽为 $2$ 的方波以 $T=4,8,16,32$ 为周期延拓。谱线 $Tc_n$ 的位置是 $n\omega_0$，间隔 $\omega_0=\frac{2\pi}{T}$ 越来越小，而顶点始终落在同一条包络 $F(\omega)=\frac{2\sin\omega}{\omega}$ 上。$T\to\infty$ 时离散谱线铺满成连续频谱，求和变成积分。</figcaption>
</figure>

上面的推导只是形式上的，严格的结论是：

!!! theorem "定理 8.2（傅里叶积分定理）"

    设 $f(t)$ 在 $(-\infty,+\infty)$ 上满足：

    1. 在任一有限区间上满足狄利克雷条件；
    2. 绝对可积，即 $\int_{-\infty}^{+\infty}|f(t)|\,\mathrm dt<+\infty$，

    则 (8.3) 式在 $f(t)$ 的连续点处成立；在间断点处，(8.3) 的左端应为 $\frac12\left[f(t+0)+f(t-0)\right]$。

(8.3) 称为 $f(t)$ 的**傅里叶积分公式**。把它拆成两步，就得到一对互逆的变换：

!!! definition "定义 8.1（傅里叶变换）"

    $$F(\omega)=\mathcal F[f(t)]=\int_{-\infty}^{+\infty}f(t)\,e^{-j\omega t}\,\mathrm dt,\tag{8.4}$$

    $$f(t)=\mathcal F^{-1}[F(\omega)]=\frac{1}{2\pi}\int_{-\infty}^{+\infty}F(\omega)\,e^{j\omega t}\,\mathrm d\omega .\tag{8.5}$$

    (8.4) 称为 $f(t)$ 的**傅里叶变换**（傅氏变换），$F(\omega)$ 称为**象函数**；(8.5) 称为 $F(\omega)$ 的**傅里叶逆变换**，$f(t)$ 称为**象原函数**。二者构成一个**傅里叶变换对**，记作 $f(t)\leftrightarrow F(\omega)$。

**频谱。** 与 (8.2) 对比可知 $F(\omega)$ 是离散频谱 $Tc_n$ 的连续版本，所以称 $F(\omega)$ 为 $f(t)$ 的**频谱密度函数**（简称**频谱**），$|F(\omega)|$ 为**振幅谱**，$\arg F(\omega)$ 为**相位谱**。(8.5) 说明：非周期信号是频率连续分布的简谐波 $e^{j\omega t}$ 的“积分和”，频率 $\omega$ 处的“密度”是 $\frac{F(\omega)}{2\pi}$。

当 $f$ 为实函数时，$F(-\omega)=\overline{F(\omega)}$，所以**振幅谱是偶函数、相位谱是奇函数**，画图时只需看 $\omega\ge0$ 的一半。

!!! example "例题 8.2（教材例 8.2：矩形脉冲）"

    求矩形脉冲 $f(t)=\begin{cases}1,&|t|\le\delta,\\0,&|t|>\delta\end{cases}$（$\delta>0$）的傅里叶变换及其积分表达式。

    **解。**

    $$F(\omega)=\int_{-\delta}^{\delta}e^{-j\omega t}\,\mathrm dt=\frac{e^{-j\omega t}}{-j\omega}\bigg|_{-\delta}^{\delta}=\frac{e^{j\omega\delta}-e^{-j\omega\delta}}{j\omega}=\frac{2\sin\delta\omega}{\omega}=2\delta\,\frac{\sin\delta\omega}{\delta\omega},$$

    $\omega=0$ 处取极限值 $F(0)=2\delta$。振幅谱 $|F(\omega)|=2\delta\left|\frac{\sin\delta\omega}{\delta\omega}\right|$；$F$ 为实数，相位谱在 $F>0$ 处为 $0$、$F<0$ 处为 $\pi$。

    由 (8.5)，并注意 $F$ 是偶函数（$\sin\omega t$ 部分积分为零）：

    $$f(t)=\frac1{2\pi}\int_{-\infty}^{+\infty}\frac{2\sin\delta\omega}{\omega}\cos\omega t\,\mathrm d\omega=\frac2\pi\int_0^{+\infty}\frac{\sin\delta\omega\cos\omega t}{\omega}\,\mathrm d\omega=\begin{cases}1,&|t|<\delta,\\\frac12,&|t|=\delta,\\0,&|t|>\delta.\end{cases}$$

    其中 $|t|=\delta$ 是间断点，取左右极限的平均值 $\frac12$。令 $t=0$、$\delta=1$，得到**狄利克雷积分** $\int_0^{+\infty}\frac{\sin x}{x}\,\mathrm dx=\frac\pi2$——和例题 5.13 用留数算出的结果一样。

!!! example "例题 8.3（教材例 8.3：理想低通频谱）"

    已知 $F(\omega)=\begin{cases}1,&|\omega|\le\alpha,\\0,&|\omega|>\alpha\end{cases}$（$\alpha>0$），求 $f(t)$。

    **解。**

    $$f(t)=\frac1{2\pi}\int_{-\alpha}^{\alpha}e^{j\omega t}\,\mathrm d\omega=\frac{1}{2\pi}\cdot\frac{e^{j\alpha t}-e^{-j\alpha t}}{jt}=\frac{\sin\alpha t}{\pi t}=\frac{\alpha}{\pi}\operatorname{Sa}(\alpha t),$$

    其中 $\operatorname{Sa}(x)=\frac{\sin x}{x}$ 称为**抽样函数**（$\operatorname{Sa}(0)=1$）。

    例题 8.2 与 8.3 互为“镜像”：时域的矩形对应频域的 $\operatorname{Sa}$ 函数，频域的矩形对应时域的 $\operatorname{Sa}$ 函数。这是傅里叶变换**对称性**的一个例子（见 8.3.1 节）。频谱只占 $|\omega|\le\alpha$ 的信号称为**带限信号**，$\frac{\sin\alpha t}{\pi t}$ 是理想低通滤波器的冲激响应，它在 $t<0$ 时不为零，所以理想低通滤波器物理上不可实现（非因果）。

!!! example "例题 8.4（教材例 8.4：单边指数衰减函数）"

    求 $f(t)=\begin{cases}e^{-\beta t},&t\ge0,\\0,&t<0\end{cases}$（$\beta>0$）的傅里叶变换，并求振幅谱和相位谱。

    **解。**

    $$F(\omega)=\int_0^{+\infty}e^{-\beta t}e^{-j\omega t}\,\mathrm dt=\frac{-1}{\beta+j\omega}e^{-(\beta+j\omega)t}\bigg|_0^{+\infty}=\frac{1}{\beta+j\omega}=\frac{\beta-j\omega}{\beta^2+\omega^2}.$$

    （上限处 $\left|e^{-(\beta+j\omega)t}\right|=e^{-\beta t}\to0$。）振幅谱 $|F(\omega)|=\dfrac{1}{\sqrt{\beta^2+\omega^2}}$，相位谱 $\arg F(\omega)=-\arctan\dfrac{\omega}{\beta}$。

    这正是 $RC$ 低通电路的频率特性的形状：低频几乎无衰减，$\omega=\beta$ 时振幅降到 $\frac{1}{\sqrt2}$ 倍（$-3\,\mathrm{dB}$），相位滞后 $\frac\pi4$。

!!! example "例题 8.5（双边指数函数；用留数求逆变换）"

    设 $\beta>0$。(1) 求 $f(t)=e^{-\beta|t|}$ 的傅里叶变换；(2) 用留数定理验证逆变换公式。

    **解。** (1) 把积分拆成 $t<0$、$t>0$ 两段，由例题 8.4，

    $$F(\omega)=\int_{-\infty}^{0}e^{\beta t}e^{-j\omega t}\,\mathrm dt+\int_0^{+\infty}e^{-\beta t}e^{-j\omega t}\,\mathrm dt=\frac{1}{\beta-j\omega}+\frac{1}{\beta+j\omega}=\frac{2\beta}{\beta^2+\omega^2}.$$

    (2) 计算 $\displaystyle\frac{1}{2\pi}\int_{-\infty}^{+\infty}\frac{2\beta}{\beta^2+\omega^2}e^{j\omega t}\,\mathrm d\omega$。把 $\omega$ 看成复变量，$R(\omega)=\frac{2\beta}{\beta^2+\omega^2}$ 的分母比分子高两次，极点为 $\omega=\pm j\beta$。

    - $t>0$：这是第 5 章公式 III 的情形（$a=t>0$），在**上半平面**补大圆弧，由约当引理弧上积分趋于零，只有 $\omega=j\beta$ 一个极点：

        $$\frac1{2\pi}\cdot2\pi j\operatorname{Res}\left[\frac{2\beta\,e^{j\omega t}}{(\omega-j\beta)(\omega+j\beta)},j\beta\right]=j\cdot\frac{2\beta e^{-\beta t}}{2j\beta}=e^{-\beta t}.$$

    - $t<0$：上半平面内 $|e^{j\omega t}|=e^{-t\operatorname{Im}\omega}$ 会变大，约当引理不能用，改在**下半平面**补大圆弧（顺时针，前面多一个负号），只有 $\omega=-j\beta$：

        $$\frac1{2\pi}\cdot(-2\pi j)\cdot\frac{2\beta e^{\beta t}}{-2j\beta}=e^{\beta t}.$$

    合起来正是 $e^{-\beta|t|}$。

!!! note "前后联系：傅里叶积分就是第 5 章的“第三类积分”"

    逆变换 $\frac1{2\pi}\int F(\omega)e^{j\omega t}\,\mathrm d\omega$ 正是第 5 章 5.3.3 节 $\int R(x)e^{iax}\,\mathrm dx$ 型的积分（$a$ 换成 $t$）。规律是：**$t>0$ 时在上半平面闭合，$t<0$ 时在下半平面闭合**。若 $F(\omega)$ 的极点全在上半平面（例如 $\frac{1}{\beta+j\omega}$ 的极点 $\omega=j\beta$），那么 $t<0$ 时下半平面里没有极点，$f(t)=0$：信号是**因果**的。

---

## 8.2 单位冲激函数（$\delta$ 函数）

### 8.2.1 单位冲激函数的概念及性质

**为什么需要 $\delta$ 函数。** 一方面，常数、单位阶跃、正弦、余弦这些最常用的函数都不是绝对可积的，按 (8.4) 它们没有傅里叶变换；另一方面，工程中有许多“瞬时”的量：冲击力、脉冲电压、集中在一点的质量或电荷，它们都没法用普通函数描述。

**引例。** 质量为 $m$ 的均匀细杆放在 $[0,\varepsilon]$ 上，线密度为 $\rho_\varepsilon(x)=\frac m\varepsilon$（$0\le x\le\varepsilon$），其余为 $0$。令 $\varepsilon\to0$ 得到原点处的质点，它的“密度”在 $x\ne0$ 处为 $0$、在 $x=0$ 处为 $\infty$，并且积分仍应等于质量 $m$。普通函数做不到“只在一点非零而积分不为零”，于是引入：

!!! definition "定义 8.2（单位冲激函数）"

    满足下列两个条件的 $\delta(t)$ 称为**单位冲激函数**（狄拉克函数，$\delta$ 函数）：

    1. 当 $t\ne0$ 时，$\delta(t)=0$；
    2. $\displaystyle\int_{-\infty}^{+\infty}\delta(t)\,\mathrm dt=1$。

    于是引例中质点的密度就是 $m\,\delta(x)$。

一种直观理解是把 $\delta(t)$ 看成宽 $\varepsilon$、高 $\frac1\varepsilon$ 的矩形脉冲 $\delta_\varepsilon(t)$ 当 $\varepsilon\to0$ 的极限：$\delta(t)=\lim\limits_{\varepsilon\to0}\delta_\varepsilon(t)$。

<figure class="video" markdown>
<video controls preload="metadata" width="100%" poster="assets/posters/C8_DeltaLimit.jpg">
<source src="assets/videos/C8_DeltaLimit.mp4" type="video/mp4">
</video>
<figcaption>动画 8-2　矩形脉冲 $\delta_\varepsilon(t)$ 越来越窄、越来越高，面积始终是 $1$；它的频谱 $\frac{\sin(\omega\varepsilon/2)}{\omega\varepsilon/2}$ 越来越平，趋于常数 $1$。和光滑函数 $f(t)$ 相乘再积分，结果是 $f$ 在小区间上的平均值，趋于 $f(0)$（筛选性质）。</figcaption>
</figure>

!!! warning "$\delta$ 函数不是普通函数"

    $\delta(t)$ 不是经典意义下的函数（没有哪个普通函数满足定义 8.2），而是**广义函数**。不能问“$\delta(0)$ 等于多少”，也不能像普通函数那样逐点使用它，而要**通过它的性质（在积分中的作用）来使用**。严格理论（分布论）中 $\delta$ 定义为把检验函数 $f$ 映为 $f(0)$ 的线性泛函，下面的筛选性质其实就是它的定义。

!!! theorem "性质 8.1（筛选性质）"

    若 $f(t)$ 在 $(-\infty,+\infty)$ 上有界，且在 $t=t_0$ 处连续，则

    $$\int_{-\infty}^{+\infty}\delta(t-t_0)\,f(t)\,\mathrm dt=f(t_0).\tag{8.6}$$

（用 $\delta_\varepsilon$ 理解：$\int\delta_\varepsilon(t)f(t)\,\mathrm dt=\frac1\varepsilon\int_{-\varepsilon/2}^{\varepsilon/2}f(t)\,\mathrm dt$ 是 $f$ 在 $\left[-\frac\varepsilon2,\frac\varepsilon2\right]$ 上的平均值，由 $f$ 的连续性，$\varepsilon\to0$ 时趋于 $f(0)$。）

!!! theorem "性质 8.2 与 8.3（偶函数；与单位阶跃函数的关系）"

    - $\delta(t)$ 为偶函数：$\delta(-t)=\delta(t)$。
    - 设 $u(t)=\begin{cases}1,&t>0,\\0,&t<0\end{cases}$ 为**单位阶跃函数**（Heaviside 函数），则

        $$\int_{-\infty}^{t}\delta(\tau)\,\mathrm d\tau=u(t),\qquad\frac{\mathrm d}{\mathrm dt}u(t)=\delta(t).$$

阶跃函数在 $t=0$ 处跳了 $1$，它的“导数”就是一个强度为 $1$ 的冲激。

**几个常用推论**（都可以由筛选性质验证）：

| 关系 | 说明 |
|---|---|
| $f(t)\,\delta(t-t_0)=f(t_0)\,\delta(t-t_0)$ | 乘积只“看到” $f$ 在 $t_0$ 的值 |
| $\delta(at)=\dfrac{1}{\lvert a\rvert}\delta(t)$（$a\ne0$） | 换元 $x=at$；时间压缩，冲激强度变为 $\frac1{\lvert a\rvert}$ |
| $\displaystyle\int_{-\infty}^{+\infty}\delta'(t)f(t)\,\mathrm dt=-f'(0)$ | 分部积分；更一般 $\int\delta^{(n)}f\,\mathrm dt=(-1)^nf^{(n)}(0)$ |
| $A\,\delta(t-t_0)$ | 发生在 $t_0$、**冲激强度**为 $A$ 的冲激，图上画成高为 $A$ 的箭头 |

### 8.2.2 $\delta$ 函数的傅里叶变换

由筛选性质立即得到

$$\mathcal F[\delta(t)]=\int_{-\infty}^{+\infty}\delta(t)e^{-j\omega t}\,\mathrm dt=e^{-j\omega t}\big|_{t=0}=1,\qquad\text{即}\quad\delta(t)\leftrightarrow1 .$$

单位冲激函数包含**所有频率成分，且幅度相等**，称为**均匀频谱**或**白色频谱**。反过来，由逆变换公式应有 $\mathcal F^{-1}[1]=\frac1{2\pi}\int_{-\infty}^{+\infty}e^{j\omega t}\,\mathrm d\omega=\delta(t)$，于是得到一个非常重要的公式：

!!! success "重要公式"

    $$\int_{-\infty}^{+\infty}e^{j\omega t}\,\mathrm d\omega=2\pi\,\delta(t),\qquad\text{同理}\quad\int_{-\infty}^{+\infty}e^{-j\omega t}\,\mathrm dt=2\pi\,\delta(\omega).\tag{8.7}$$

左边的积分在普通意义下是发散的。这里的积分以及由此得到的变换都是在 $\delta$ 函数（广义函数）意义下理解的，称为**广义傅里叶变换**。有了它，常数、阶跃、正弦、余弦和周期函数就都有傅里叶变换了。

!!! example "例题 8.6（教材例 8.7、8.9：常数与正弦型函数）"

    求 (1) $f_1(t)=1$；(2) $f_2(t)=e^{j\omega_0t}$；(3) $f_3(t)=\cos\omega_0t$；(4) $f_4(t)=\sin\omega_0t$ 的傅里叶变换。

    **解。** 都用 (8.7)：

    (1) $F_1(\omega)=\int_{-\infty}^{+\infty}e^{-j\omega t}\,\mathrm dt=2\pi\delta(\omega)$。

    (2) $F_2(\omega)=\int_{-\infty}^{+\infty}e^{-j(\omega-\omega_0)t}\,\mathrm dt=2\pi\delta(\omega-\omega_0)$。

    (3) $\cos\omega_0t=\frac12\left(e^{j\omega_0t}+e^{-j\omega_0t}\right)$，所以 $F_3(\omega)=\pi\left[\delta(\omega+\omega_0)+\delta(\omega-\omega_0)\right]$。

    (4) $\sin\omega_0t=\frac1{2j}\left(e^{j\omega_0t}-e^{-j\omega_0t}\right)$，所以 $F_4(\omega)=\frac{\pi}{j}\left[\delta(\omega-\omega_0)-\delta(\omega+\omega_0)\right]=j\pi\left[\delta(\omega+\omega_0)-\delta(\omega-\omega_0)\right]$。

    **解读。** 直流信号只含频率 $0$；单一频率的余弦只在 $\pm\omega_0$ 处有两根强度为 $\pi$ 的谱线。冲激函数可以把**离散频谱**和**连续频谱**统一在同一个框架里。

!!! example "例题 8.7（教材例 8.8：单位阶跃函数）"

    证明 $\mathcal F[u(t)]=\dfrac{1}{j\omega}+\pi\delta(\omega)$。

    **证。** 求 $F(\omega)=\frac{1}{j\omega}+\pi\delta(\omega)$ 的逆变换：

    $$\mathcal F^{-1}[F]=\frac1{2\pi}\int_{-\infty}^{+\infty}\pi\delta(\omega)e^{j\omega t}\,\mathrm d\omega+\frac1{2\pi}\int_{-\infty}^{+\infty}\frac{e^{j\omega t}}{j\omega}\,\mathrm d\omega=\frac12+\frac1\pi\int_0^{+\infty}\frac{\sin\omega t}{\omega}\,\mathrm d\omega .$$

    （$\frac{\cos\omega t}{j\omega}$ 是奇函数，按主值积分为零；$\frac{j\sin\omega t}{j\omega}$ 是偶函数。）由狄利克雷积分，$\int_0^{+\infty}\frac{\sin\omega t}{\omega}\,\mathrm d\omega$ 在 $t>0$ 时为 $\frac\pi2$、$t<0$ 时为 $-\frac\pi2$，所以逆变换在 $t>0$ 时为 $1$、$t<0$ 时为 $0$，正是 $u(t)$。$\blacksquare$

    **为什么多了 $\pi\delta(\omega)$。** $u(t)=\frac12+\frac12\operatorname{sgn}t$，直流分量 $\frac12$ 贡献 $\frac12\cdot2\pi\delta(\omega)=\pi\delta(\omega)$，而 $\mathcal F[\operatorname{sgn}t]=\frac{2}{j\omega}$（练习 2）贡献 $\frac1{j\omega}$。

!!! theorem "定理 8.3（周期函数的傅里叶变换）"

    设 $f(t)$ 是以 $T$ 为周期的函数，在 $\left[-\frac T2,\frac T2\right]$ 上满足狄利克雷条件，$\omega_0=\frac{2\pi}{T}$，$F(n\omega_0)=c_n$ 为它的离散频谱，则

    $$\mathcal F[f(t)]=2\pi\sum_{n=-\infty}^{+\infty}F(n\omega_0)\,\delta(\omega-n\omega_0).$$

（对 (8.2) 逐项变换，用例题 8.6(2) 即得。）周期信号的频谱是一串位于 $n\omega_0$ 处的冲激，强度是傅里叶系数的 $2\pi$ 倍：**离散频谱就是由冲激组成的连续频谱**。

---

## 8.3 傅里叶变换的性质

以下设所涉及函数的傅里叶变换都存在（可以是广义的），$F(\omega)=\mathcal F[f(t)]$，$G(\omega)=\mathcal F[g(t)]$；对求导、积分、求和与极限的交换次序不另作说明。

### 8.3.1 基本性质

!!! theorem "性质 8.4（线性性质）"

    设 $\alpha,\beta$ 为常数，则

    $$\mathcal F[\alpha f(t)+\beta g(t)]=\alpha F(\omega)+\beta G(\omega),\qquad\mathcal F^{-1}[\alpha F(\omega)+\beta G(\omega)]=\alpha f(t)+\beta g(t).$$

!!! theorem "性质 8.5（位移性质）"

    设 $t_0,\omega_0$ 为实常数，则

    $$\mathcal F[f(t-t_0)]=e^{-j\omega t_0}F(\omega)\quad\text{（时移）},\qquad\mathcal F^{-1}[F(\omega-\omega_0)]=e^{j\omega_0t}f(t)\quad\text{（频移）}.$$

??? proof "证明"

    令 $x=t-t_0$：$\mathcal F[f(t-t_0)]=\int_{-\infty}^{+\infty}f(x)e^{-j\omega(x+t_0)}\,\mathrm dx=e^{-j\omega t_0}F(\omega)$。频移性质同理。$\blacksquare$

**物理意义。** 信号延迟 $t_0$，各频率分量的**振幅不变**，只是相位各滞后 $\omega t_0$（与频率成正比，即“线性相位”）。频移性质是**调制**的数学基础：乘以 $e^{j\omega_0t}$ 把整个频谱搬到 $\omega_0$ 附近。

!!! success "推论（调制定理）"

    $$\mathcal F[f(t)\cos\omega_0t]=\frac12\left[F(\omega-\omega_0)+F(\omega+\omega_0)\right],\qquad\mathcal F[f(t)\sin\omega_0t]=\frac{1}{2j}\left[F(\omega-\omega_0)-F(\omega+\omega_0)\right].$$

    （由 $\cos\omega_0t=\frac{e^{j\omega_0t}+e^{-j\omega_0t}}{2}$ 与频移性质。）调幅广播就是把语音频谱搬到载频 $\pm\omega_0$ 附近。

!!! example "例题 8.8（教材例 8.10 改编）"

    已知 $G(\omega)=\dfrac{e^{-j2\omega}}{\beta+j\omega}$（$\beta>0$），求 $g(t)=\mathcal F^{-1}[G(\omega)]$。

    **解。** 由例题 8.4，$\mathcal F^{-1}\left[\frac{1}{\beta+j\omega}\right]=e^{-\beta t}u(t)$；$e^{-j2\omega}$ 对应延迟 $t_0=2$。由时移性质，

    $$g(t)=e^{-\beta(t-2)}u(t-2)=\begin{cases}e^{-\beta(t-2)},&t>2,\\0,&t<2.\end{cases}$$

!!! theorem "性质 8.6（相似性质 / 尺度变换）"

    设 $a$ 为非零实常数，则

    $$\mathcal F[f(at)]=\frac{1}{|a|}F\!\left(\frac{\omega}{a}\right).\tag{8.8}$$

??? proof "证明"

    令 $x=at$。$a>0$ 时，$\mathcal F[f(at)]=\frac1a\int_{-\infty}^{+\infty}f(x)e^{-j\frac{\omega}{a}x}\,\mathrm dx=\frac1aF\!\left(\frac\omega a\right)$；$a<0$ 时换元后积分上下限互换，多出一个负号，得 $-\frac1aF\!\left(\frac\omega a\right)$。两种情况合起来就是 (8.8)。$\blacksquare$

**物理意义。** $a>1$ 时信号在时间上**压缩**（快放），频谱**展宽**、幅度变小；$0<a<1$ 时信号**拉长**（慢放），频谱**变窄**。信号的持续时间与频带宽度成反比，二者不能同时任意小。这就是信号处理里的**时频测不准关系**，和量子力学的不确定性原理在数学上是同一件事。

<figure class="video" markdown>
<video controls preload="metadata" width="100%" poster="assets/posters/C8_Scaling.jpg">
<source src="assets/videos/C8_Scaling.mp4" type="video/mp4">
</video>
<figcaption>动画 8-3　高斯信号 $e^{-a^2t^2}$ 的频谱是 $\frac{\sqrt\pi}{a}e^{-\omega^2/(4a^2)}$（例题 8.15）。$a$ 增大时时域变窄、频域变宽、峰值降低；$a$ 减小时相反。橙色横线标出半高宽，两者的乘积不随 $a$ 改变。</figcaption>
</figure>

!!! example "例题 8.9（教材例 8.11）"

    已知抽样信号 $f(t)=\dfrac{\sin2t}{\pi t}$ 的频谱为 $F(\omega)=\begin{cases}1,&|\omega|\le2,\\0,&|\omega|>2,\end{cases}$ 求 $g(t)=f\!\left(\frac t2\right)$ 的频谱 $G(\omega)$。

    **解。** 由 (8.8)，$a=\frac12$：

    $$G(\omega)=2F(2\omega)=\begin{cases}2,&|\omega|\le1,\\0,&|\omega|>1.\end{cases}$$

    $g(t)$ 是 $f(t)$ 拉宽一倍的结果，频带从 $|\omega|<2$ 缩小到 $|\omega|<1$，高度加倍（面积 $\int G\,\mathrm d\omega=2\pi g(0)$ 也加倍）。

!!! theorem "性质 8.7（微分性质）"

    1. **导数的象函数**：若 $\lim\limits_{|t|\to+\infty}f(t)=0$，则 $\mathcal F[f'(t)]=j\omega\,F(\omega)$。一般地，若 $\lim\limits_{|t|\to+\infty}f^{(k)}(t)=0$（$k=0,1,\dots,n-1$），则

        $$\mathcal F[f^{(n)}(t)]=(j\omega)^nF(\omega).\tag{8.9}$$

    2. **象函数的导数**：

        $$F'(\omega)=\mathcal F[-jt\,f(t)],\qquad F^{(n)}(\omega)=\mathcal F\left[(-jt)^nf(t)\right].\tag{8.10}$$

??? proof "证明"

    分部积分：

    $$\mathcal F[f'(t)]=f(t)e^{-j\omega t}\Big|_{-\infty}^{+\infty}+j\omega\int_{-\infty}^{+\infty}f(t)e^{-j\omega t}\,\mathrm dt=j\omega F(\omega),$$

    因为 $|f(t)e^{-j\omega t}|=|f(t)|\to0$。反复使用即得 (8.9)。(8.10) 只需在 $F(\omega)=\int f(t)e^{-j\omega t}\,\mathrm dt$ 两边对 $\omega$ 求导（积分号下求导），被积函数多出因子 $-jt$。$\blacksquare$

**这是傅里叶变换最有用的性质。** 时域求导等于频域乘以 $j\omega$，于是线性常系数微分方程在频域里变成代数方程。电路里电感的阻抗 $j\omega L$、电容的阻抗 $\frac{1}{j\omega C}$ 就是这样来的：$u=L\frac{\mathrm di}{\mathrm dt}$ 变换后成为 $U=j\omega L\,I$。

**(8.10) 的记法。** 直接对 $F(\omega)=\int f(t)e^{-j\omega t}\,\mathrm dt$ 关于 $\omega$ 求导，每求一次导就“掉下来”一个 $-jt$。常用的形式是 $\mathcal F[t^nf(t)]=j^nF^{(n)}(\omega)$。

!!! example "例题 8.10（教材 PPT 例题）"

    设 $f(t)=t^2\cos t$，求 $\mathcal F[f(t)]$。

    **解。** 令 $g(t)=\cos t$，由例题 8.6，$G(\omega)=\pi\left[\delta(\omega+1)+\delta(\omega-1)\right]$。由 (8.10)，$G''(\omega)=\mathcal F[(-jt)^2g(t)]=-\mathcal F[t^2g(t)]$，所以

    $$\mathcal F[t^2\cos t]=-G''(\omega)=-\pi\left[\delta''(\omega+1)+\delta''(\omega-1)\right].$$

!!! theorem "性质 8.8（积分性质）"

    设 $g(t)=\displaystyle\int_{-\infty}^{t}f(\tau)\,\mathrm d\tau$。若 $\lim\limits_{t\to+\infty}g(t)=0$（即 $F(0)=\int_{-\infty}^{+\infty}f\,\mathrm dt=0$），则

    $$\mathcal F\left[\int_{-\infty}^{t}f(\tau)\,\mathrm d\tau\right]=\frac{1}{j\omega}F(\omega).$$

（$g'=f$，对 $g$ 用微分性质：$j\omega\,\mathcal F[g]=F$。）若 $F(0)\ne0$，$g(t)$ 在 $+\infty$ 处趋于常数 $F(0)$，结论修正为 $\mathcal F[g]=\frac{F(\omega)}{j\omega}+\pi F(0)\delta(\omega)$；取 $f=\delta$ 就回到 $\mathcal F[u(t)]=\frac1{j\omega}+\pi\delta(\omega)$。

!!! theorem "性质 8.9（帕塞瓦尔等式）"

    $$\int_{-\infty}^{+\infty}f^2(t)\,\mathrm dt=\frac{1}{2\pi}\int_{-\infty}^{+\infty}|F(\omega)|^2\,\mathrm d\omega\qquad(f\ \text{为实函数}).\tag{8.11}$$

??? proof "证明"

    $f$ 为实函数时 $\overline{F(\omega)}=\int f(t)e^{j\omega t}\,\mathrm dt$，所以

    $$\frac1{2\pi}\int|F|^2\,\mathrm d\omega=\frac1{2\pi}\int F(\omega)\left[\int f(t)e^{j\omega t}\,\mathrm dt\right]\mathrm d\omega=\int f(t)\left[\frac1{2\pi}\int F(\omega)e^{j\omega t}\,\mathrm d\omega\right]\mathrm dt=\int f^2(t)\,\mathrm dt.\ \blacksquare$$

**物理意义。** 左边是信号的总能量（把 $f$ 看成 $1\,\Omega$ 电阻上的电压），右边说总能量也可以按频率分开来算，$\frac{|F(\omega)|^2}{2\pi}$ 是单位频带内的能量，称为**能量谱密度**。

!!! example "例题 8.11（教材例 8.12）"

    求 $\displaystyle\int_0^{+\infty}\frac{\sin^2\omega}{\omega^2}\,\mathrm d\omega$。

    **解。** 由例题 8.2（$\delta=1$），$f(t)=\begin{cases}1,&|t|\le1\\0,&|t|>1\end{cases}$ 的频谱为 $F(\omega)=\frac{2\sin\omega}{\omega}$。由 (8.11)，

    $$\int_{-\infty}^{+\infty}\frac{4\sin^2\omega}{\omega^2}\,\mathrm d\omega=2\pi\int_{-1}^{1}1^2\,\mathrm dt=4\pi,$$

    被积函数为偶函数，故 $\displaystyle\int_0^{+\infty}\frac{\sin^2\omega}{\omega^2}\,\mathrm d\omega=\frac\pi2$。

!!! tip "补充：对称性（对偶性）"

    把 (8.4)、(8.5) 对照着看，两个公式只差一个 $2\pi$ 和指数的符号。由此可得：若 $f(t)\leftrightarrow F(\omega)$，则

    $$\mathcal F[F(t)]=2\pi f(-\omega).$$

    例如 $\delta(t)\leftrightarrow1$ 立即给出 $1\leftrightarrow2\pi\delta(\omega)$；矩形脉冲 $\leftrightarrow$ $\operatorname{Sa}$ 函数立即给出 $\operatorname{Sa}$ 函数 $\leftrightarrow$ 矩形（例题 8.2 与 8.3）。记住一半的变换对，另一半可以直接写出来。

### 8.3.2 卷积与卷积定理

!!! definition "定义 8.3（卷积）"

    设 $f_1(t),f_2(t)$ 在 $(-\infty,+\infty)$ 上有定义。若对每个实数 $t$，广义积分 $\int_{-\infty}^{+\infty}f_1(\tau)f_2(t-\tau)\,\mathrm d\tau$ 都收敛，则称它定义的关于 $t$ 的函数为 $f_1$ 与 $f_2$ 的**卷积**，记作

    $$f_1(t)*f_2(t)=\int_{-\infty}^{+\infty}f_1(\tau)\,f_2(t-\tau)\,\mathrm d\tau .$$

卷积满足**交换律** $f_1*f_2=f_2*f_1$（换元 $\tau'=t-\tau$）、**结合律** $f_1*(f_2*f_3)=(f_1*f_2)*f_3$ 与**分配律** $f_1*(f_2+f_3)=f_1*f_2+f_1*f_3$。

**几何意义：翻转、平移、相乘、积分。** 把 $f_2(\tau)$ 关于纵轴翻转得 $f_2(-\tau)$，再平移 $t$ 得 $f_2(t-\tau)$，与 $f_1(\tau)$ 相乘，乘积曲线下的面积就是 $t$ 时刻的卷积值。$t$ 从左到右滑动，就画出整条 $f_1*f_2$ 曲线，所以卷积也叫**褶积**。

<figure class="video" markdown>
<video controls preload="metadata" width="100%" poster="assets/posters/C8_Convolution.jpg">
<source src="assets/videos/C8_Convolution.mp4" type="video/mp4">
</video>
<figcaption>动画 8-4　$f(\tau)=e^{-\tau}u(\tau)$ 与翻转、平移后的 $g(t-\tau)=e^{-2(t-\tau)}u(t-\tau)$ 相乘，阴影面积随 $t$ 变化，画出 $f*g(t)=(e^{-t}-e^{-2t})u(t)$（例题 8.12，$\alpha=1,\beta=2$）。$t<0$ 时两者不重叠，卷积为零。</figcaption>
</figure>

!!! example "例题 8.12（教材例 8.13）"

    设 $f(t)=e^{-\alpha t}u(t)$，$g(t)=e^{-\beta t}u(t)$，$\alpha,\beta>0$，$\alpha\ne\beta$，求 $f*g$。

    **解。** $f(\tau)g(t-\tau)$ 只有在 $\tau\ge0$ 且 $t-\tau\ge0$ 时才不为零。

    - $t<0$ 时，不存在这样的 $\tau$，$f*g(t)=0$；
    - $t\ge0$ 时，积分区间是 $[0,t]$：

        $$f*g(t)=\int_0^te^{-\alpha\tau}e^{-\beta(t-\tau)}\,\mathrm d\tau=e^{-\beta t}\int_0^te^{(\beta-\alpha)\tau}\,\mathrm d\tau=\frac{e^{-\alpha t}-e^{-\beta t}}{\beta-\alpha}.$$

    即 $f*g(t)=\dfrac{1}{\beta-\alpha}\left(e^{-\alpha t}-e^{-\beta t}\right)u(t)$。

    **规律。** 两个**因果**函数（$t<0$ 时为零）的卷积仍是因果的，并且 $f_1*f_2(t)=\int_0^tf_1(\tau)f_2(t-\tau)\,\mathrm d\tau$。第 9 章的卷积就是这种形式。

!!! example "例题 8.13（教材例 8.14：分段函数的卷积）"

    求 $f(t)=t^2u(t)$ 与 $g(t)=\begin{cases}1,&|t|\le1,\\0,&|t|>1\end{cases}$ 的卷积。

    **解。** 用 $f*g=\int f(\tau)g(t-\tau)\,\mathrm d\tau$，被积函数非零要求 $\tau\ge0$ 且 $|t-\tau|\le1$，即 $\tau\in[\max(0,t-1),\,t+1]$。

    - $t<-1$：$t+1<0$，区间为空，$f*g=0$；
    - $-1\le t\le1$：$\tau\in[0,t+1]$，$f*g=\int_0^{t+1}\tau^2\,\mathrm d\tau=\frac13(t+1)^3$；
    - $t>1$：$\tau\in[t-1,t+1]$，$f*g=\frac13\left[(t+1)^3-(t-1)^3\right]=\frac13(6t^2+2)$。

    **方法。** 分段函数的卷积，关键是根据 $t$ 的范围定出 $\tau$ 的积分区间：画出 $f(\tau)$ 和滑动的 $g(t-\tau)$，看它们在哪一段重叠。

!!! theorem "定理 8.4（卷积定理）"

    设 $F_1(\omega)=\mathcal F[f_1(t)]$，$F_2(\omega)=\mathcal F[f_2(t)]$，则

    $$\mathcal F[f_1(t)*f_2(t)]=F_1(\omega)\cdot F_2(\omega),\qquad\mathcal F[f_1(t)\cdot f_2(t)]=\frac{1}{2\pi}F_1(\omega)*F_2(\omega).\tag{8.12}$$

??? proof "证明"

    交换积分次序，并在内层令 $x=t-\tau$：

    $$\mathcal F[f_1*f_2]=\int_{-\infty}^{+\infty}\left[\int_{-\infty}^{+\infty}f_1(\tau)f_2(t-\tau)\,\mathrm d\tau\right]e^{-j\omega t}\,\mathrm dt=\int_{-\infty}^{+\infty}f_1(\tau)e^{-j\omega\tau}\left[\int_{-\infty}^{+\infty}f_2(x)e^{-j\omega x}\,\mathrm dx\right]\mathrm d\tau=F_1(\omega)F_2(\omega).$$

    第二式同理（注意逆变换公式里的 $\frac1{2\pi}$）。$\blacksquare$

**时域卷积 = 频域相乘；时域相乘 = 频域卷积（再除以 $2\pi$）。** 卷积积分很难算，乘积很容易，所以卷积定理既能用来算卷积，也能用来算某些函数的傅里叶变换。

!!! example "例题 8.14（教材例 8.15、8.16）"

    (1) 设 $\alpha,\beta>0$，求 $f(t)=\dfrac{\sin\alpha t}{\pi t}$ 与 $g(t)=\dfrac{\sin\beta t}{\pi t}$ 的卷积；(2) 求 $f(t)=e^{-\beta t}u(t)\cos\omega_0t$（$\beta>0$）的傅里叶变换。

    **解。** (1) 由例题 8.3，$F(\omega)$、$G(\omega)$ 分别是 $|\omega|\le\alpha$、$|\omega|\le\beta$ 上取 $1$ 的矩形。乘积 $F\cdot G$ 是 $|\omega|\le\gamma$ 上取 $1$ 的矩形，$\gamma=\min(\alpha,\beta)$。由卷积定理与例题 8.3，

    $$f*g=\mathcal F^{-1}[FG]=\frac{\sin\gamma t}{\pi t}.$$

    直接计算这个卷积积分非常困难，用卷积定理一步就得到了。

    (2) 由 (8.12) 第二式，$\mathcal F[e^{-\beta t}u(t)]=\frac{1}{\beta+j\omega}$，$\mathcal F[\cos\omega_0t]=\pi[\delta(\omega+\omega_0)+\delta(\omega-\omega_0)]$，再用筛选性质：

    $$\mathcal F[f]=\frac{1}{2\pi}\int_{-\infty}^{+\infty}\frac{\pi\left[\delta(\tau+\omega_0)+\delta(\tau-\omega_0)\right]}{\beta+j(\omega-\tau)}\,\mathrm d\tau=\frac12\left[\frac{1}{\beta+j(\omega+\omega_0)}+\frac{1}{\beta+j(\omega-\omega_0)}\right]=\frac{\beta+j\omega}{(\beta+j\omega)^2+\omega_0^2}.$$

    （用调制定理也能直接得到同样结果。）

**$\delta$ 函数是卷积的“单位元”。** 由筛选性质，$h(t)*\delta(t)=\int h(\tau)\delta(t-\tau)\,\mathrm d\tau=h(t)$，更一般地 $h(t)*\delta(t-t_0)=h(t-t_0)$。用卷积定理也能看出来：$\mathcal F[\delta]=1$，乘 $1$ 不改变 $H(\omega)$。

!!! note "前后联系：线性时不变系统"

    信号与系统课程里，线性时不变系统由它的**冲激响应** $h(t)$（输入 $\delta(t)$ 时的输出）完全决定：任意输入 $x(t)$ 的输出是

    $$y(t)=x(t)*h(t)\qquad\Longleftrightarrow\qquad Y(\omega)=X(\omega)\,H(\omega).$$

    $H(\omega)=\mathcal F[h(t)]$ 称为**频率响应**。输入是 $e^{j\omega_0t}$ 时，输出是 $H(\omega_0)e^{j\omega_0t}$：系统只把每个频率分量乘上一个复数，也就是改变它的振幅和相位（第 1 章“乘复数 = 旋转 + 伸缩”）。滤波器设计就是设计 $|H(\omega)|$ 的形状。

### 8.3.3 综合举例

!!! example "例题 8.15（教材例 8.20：高斯函数，用柯西定理平移积分路径）"

    已知 $\displaystyle\int_{-\infty}^{+\infty}e^{-t^2}\,\mathrm dt=\sqrt\pi$，求 $f(t)=e^{-t^2}$ 的傅里叶变换。

    **解。** 配方：$-t^2-j\omega t=-\left(t+\frac{j\omega}{2}\right)^2-\frac{\omega^2}{4}$，所以

    $$F(\omega)=e^{-\omega^2/4}\int_{-\infty}^{+\infty}e^{-(t+j\omega/2)^2}\,\mathrm dt=e^{-\omega^2/4}\lim_{\beta\to+\infty}\int_{-\beta+j\omega/2}^{\beta+j\omega/2}e^{-z^2}\,\mathrm dz .$$

    后一个积分沿水平线 $\operatorname{Im}z=\frac\omega2$。取矩形闭路：实轴上的 $[-\beta,\beta]$、两条竖边、上边 $\operatorname{Im}z=\frac\omega2$。$e^{-z^2}$ 在全平面解析，由**柯西积分定理**（第 3 章），闭路积分为零。竖边 $z=\pm\beta+jy$（$y$ 在 $0$ 与 $\frac\omega2$ 之间）上 $|e^{-z^2}|=e^{-\beta^2+y^2}\le e^{-\beta^2+\omega^2/4}$，边长有限，由长大不等式，$\beta\to+\infty$ 时竖边上的积分趋于零。所以上边的积分等于实轴上的积分 $\sqrt\pi$：

    $$\mathcal F\left[e^{-t^2}\right]=\sqrt\pi\,e^{-\omega^2/4}.$$

    **结论。** 高斯函数的傅里叶变换仍是高斯函数。再由相似性质，$\mathcal F\left[e^{-a^2t^2}\right]=\frac{\sqrt\pi}{a}e^{-\omega^2/(4a^2)}$（动画 8-3）。

!!! example "例题 8.16（教材例 8.19：用傅里叶变换证明积分等式）"

    证明 $\displaystyle\int_0^{+\infty}\frac{\omega^2+2}{\omega^4+4}\cos\omega t\,\mathrm d\omega=\frac\pi2e^{-|t|}\cos t$。

    **证。** 令 $f(t)=e^{-|t|}\cos t$。由例题 8.5（$\beta=1$），$\mathcal F[e^{-|t|}]=\frac{2}{1+\omega^2}$；由调制定理，

    $$F(\omega)=\frac12\left[\frac{2}{1+(\omega-1)^2}+\frac{2}{1+(\omega+1)^2}\right]=\frac{1}{\omega^2-2\omega+2}+\frac{1}{\omega^2+2\omega+2}=\frac{2(\omega^2+2)}{\omega^4+4}.$$

    $F$ 为偶函数，由逆变换公式，

    $$f(t)=\frac1{2\pi}\int_{-\infty}^{+\infty}F(\omega)\cos\omega t\,\mathrm d\omega=\frac2\pi\int_0^{+\infty}\frac{\omega^2+2}{\omega^4+4}\cos\omega t\,\mathrm d\omega,$$

    移项即得所证。$\blacksquare$

    **方法。** 用傅里叶变换计算（或证明）一个含参数的积分：找一个函数 $f$，使它的频谱恰好是被积函数，再写出逆变换。

!!! example "例题 8.17（用留数求逆变换，判断因果性）"

    求 $F(\omega)=\dfrac{1}{(1+j\omega)(2+j\omega)}$ 的傅里叶逆变换。

    **解法一（卷积定理）。** $F=\frac{1}{1+j\omega}\cdot\frac{1}{2+j\omega}$，由例题 8.4 与 8.12，$f(t)=e^{-t}u(t)*e^{-2t}u(t)=\left(e^{-t}-e^{-2t}\right)u(t)$。

    **解法二（部分分式）。** $F=\frac{1}{1+j\omega}-\frac{1}{2+j\omega}$，逐项查表得同样结果。

    **解法三（留数）。** $(1+j\omega)(2+j\omega)=j^2(\omega-j)(\omega-2j)=-(\omega-j)(\omega-2j)$，两个极点 $\omega=j,2j$ 都在上半平面，分母比分子高两次。

    - $t>0$：上半平面闭合，

        $$f(t)=\frac{1}{2\pi}\cdot2\pi j\left\{\operatorname{Res}\left[\frac{-e^{j\omega t}}{(\omega-j)(\omega-2j)},j\right]+\operatorname{Res}\left[\frac{-e^{j\omega t}}{(\omega-j)(\omega-2j)},2j\right]\right\}=j\left(\frac{-e^{-t}}{-j}+\frac{-e^{-2t}}{j}\right)=e^{-t}-e^{-2t};$$

    - $t<0$：下半平面闭合，里面没有极点，$f(t)=0$。

    三种方法结果一致。解法三说明了一般规律：**频谱的极点全在上半 $\omega$ 平面 $\iff$ 信号是因果的**。换成 $s=j\omega$，上半 $\omega$ 平面就是左半 $s$ 平面，这正是第 9 章“极点在左半平面 $\iff$ 稳定”的来历。

!!! warning "易错点"

    - **$2\pi$ 放在哪里。** 本书约定正变换没有系数、逆变换带 $\frac1{2\pi}$。有的书用对称形式（两边各带 $\frac{1}{\sqrt{2\pi}}$），有的信号书用频率 $f$（$\omega=2\pi f$，正逆变换都没有系数）。换一本书做题时，先核对它的定义。
    - **指数的符号。** 正变换是 $e^{-j\omega t}$，逆变换是 $e^{+j\omega t}$；记反了会导致时移性质、因果性的结论全部反号。
    - **$\mathcal F[u(t)]$ 不是 $\frac1{j\omega}$**，还有 $\pi\delta(\omega)$ 这一项（直流分量）；$\mathcal F[1]=2\pi\delta(\omega)$，不是 $\delta(\omega)$。
    - **相似性质有 $\frac1{|a|}$**，$a<0$ 时取绝对值。同时有平移和伸缩时，先把 $f(at-b)$ 写成 $f\big(a(t-\frac ba)\big)$：$\mathcal F[f(at-b)]=\frac{1}{|a|}e^{-j\omega b/a}F\!\left(\frac\omega a\right)$。
    - **频域卷积要除以 $2\pi$**：$\mathcal F[f_1f_2]=\frac1{2\pi}F_1*F_2$，时域卷积对应的乘积则没有系数。
    - **卷积的积分区间**随 $t$ 改变，分段函数一定要分情况讨论（例题 8.13）。
    - 用留数求逆变换时，**$t>0$ 与 $t<0$ 要分别闭合**（上半平面 / 下半平面），下半平面闭合是顺时针方向，要多一个负号。

---

## 知识点归纳总结

!!! abstract "知识点归纳总结"

    ```mermaid
    flowchart TD
        FS["傅里叶级数<br/>f_T = Σ c_n exp(jnω₀t)<br/>离散频谱 c_n"] -- "T → ∞" --> FI["傅里叶积分公式<br/>f = (1/2π)∫[∫f exp(−jωτ) dτ] exp(jωt) dω"]
        FI --> FT["傅里叶变换对<br/>F(ω) = ∫ f exp(−jωt) dt<br/>f = (1/2π)∫ F exp(jωt) dω"]
        DL["δ 函数<br/>筛选性质 · 偶函数 · ∫δ = u"] --> GFT["广义傅里叶变换<br/>δ ↔ 1，1 ↔ 2πδ(ω)<br/>cos ω₀t、u(t)、周期函数"]
        FT --> GFT
        FT --> PR["性质<br/>线性 · 时移/频移 · 相似<br/>微分 jω · 积分 1/(jω) · 帕塞瓦尔"]
        FT --> CV["卷积 f₁*f₂<br/>卷积定理：F₁·F₂"]
        DL --> CV
        C3["第3章 柯西定理"] -.平移积分路径.-> FT
        C5["第5章 留数 + 约当引理"] -.逆变换：t>0 上半平面闭合.-> FT
        PR --> APP["应用：解微分方程、频谱分析<br/>调制、滤波 Y = HX"]
        CV --> APP
        FT -.乘衰减因子 exp(−σt).-> LAP["第9章 拉普拉斯变换"]
    ```

    **常用傅里叶变换对**（$\beta>0$）

    | $f(t)$ | $F(\omega)=\mathcal F[f(t)]$ | 出处 |
    |---|---|---|
    | 矩形脉冲 $1\ (\lvert t\rvert\le\delta)$ | $\dfrac{2\sin\delta\omega}{\omega}$ | 例题 8.2 |
    | $\dfrac{\sin\alpha t}{\pi t}$ | $1\ (\lvert\omega\rvert\le\alpha)$，否则 $0$ | 例题 8.3 |
    | $e^{-\beta t}u(t)$ | $\dfrac{1}{\beta+j\omega}$ | 例题 8.4 |
    | $te^{-\beta t}u(t)$ | $\dfrac{1}{(\beta+j\omega)^2}$ | 练习 5 |
    | $e^{-\beta\lvert t\rvert}$ | $\dfrac{2\beta}{\beta^2+\omega^2}$ | 例题 8.5 |
    | $e^{-t^2}$ | $\sqrt\pi\,e^{-\omega^2/4}$ | 例题 8.15 |
    | $\delta(t)$；$\delta(t-t_0)$ | $1$；$e^{-j\omega t_0}$ | 8.2.2 |
    | $1$ | $2\pi\delta(\omega)$ | 例题 8.6 |
    | $e^{j\omega_0t}$ | $2\pi\delta(\omega-\omega_0)$ | 例题 8.6 |
    | $\cos\omega_0t$ | $\pi[\delta(\omega+\omega_0)+\delta(\omega-\omega_0)]$ | 例题 8.6 |
    | $\sin\omega_0t$ | $j\pi[\delta(\omega+\omega_0)-\delta(\omega-\omega_0)]$ | 例题 8.6 |
    | $u(t)$ | $\dfrac{1}{j\omega}+\pi\delta(\omega)$ | 例题 8.7 |
    | $\operatorname{sgn}t$ | $\dfrac{2}{j\omega}$ | 练习 2 |
    | $e^{-\beta t}u(t)\cos\omega_0t$ | $\dfrac{\beta+j\omega}{(\beta+j\omega)^2+\omega_0^2}$ | 例题 8.14 |
    | $e^{-\beta t}u(t)\sin\omega_0t$ | $\dfrac{\omega_0}{(\beta+j\omega)^2+\omega_0^2}$ | 练习 4 |

    **性质一览**（$f\leftrightarrow F$，$g\leftrightarrow G$）

    | 性质 | 时域 | 频域 |
    |---|---|---|
    | 线性 | $\alpha f+\beta g$ | $\alpha F+\beta G$ |
    | 时移 | $f(t-t_0)$ | $e^{-j\omega t_0}F(\omega)$ |
    | 频移 | $e^{j\omega_0t}f(t)$ | $F(\omega-\omega_0)$ |
    | 调制 | $f(t)\cos\omega_0t$ | $\frac12[F(\omega-\omega_0)+F(\omega+\omega_0)]$ |
    | 相似 | $f(at)$ | $\frac1{\lvert a\rvert}F\!\left(\frac\omega a\right)$ |
    | 对称 | $F(t)$ | $2\pi f(-\omega)$ |
    | 时域微分 | $f^{(n)}(t)$ | $(j\omega)^nF(\omega)$ |
    | 频域微分 | $t^nf(t)$ | $j^nF^{(n)}(\omega)$ |
    | 积分（$F(0)=0$） | $\int_{-\infty}^tf(\tau)\,\mathrm d\tau$ | $\frac{1}{j\omega}F(\omega)$ |
    | 时域卷积 | $f*g$ | $F\cdot G$ |
    | 频域卷积 | $f\cdot g$ | $\frac1{2\pi}F*G$ |
    | 帕塞瓦尔 | $\int f^2\,\mathrm dt$ | $\frac1{2\pi}\int\lvert F\rvert^2\,\mathrm d\omega$ |

    **方法清单**

    1. **按定义求变换**：分段写出积分区间；含 $e^{-\beta|t|}$ 时拆成 $t<0$、$t>0$ 两段；验证上下限处指数项趋于零；
    2. **广义变换**：常数、$e^{j\omega_0t}$、正弦、余弦、阶跃一律化成 $e^{j\omega_0t}$，用 $\int e^{-j\omega t}\,\mathrm dt=2\pi\delta(\omega)$；
    3. **用性质求变换**：先找一个已知的“母函数”，再看是时移、频移（调制）、伸缩、乘 $t^n$ 还是求导，依次套用；
    4. **求逆变换**：有理函数优先部分分式 + 查表；也可用留数（$t>0$ 上半平面、$t<0$ 下半平面）；乘积形式用卷积定理；
    5. **卷积**：画图定积分区间，分段计算；或者变到频域相乘后再逆变换；
    6. **计算积分**：用逆变换公式（例题 8.2、8.16）或帕塞瓦尔等式（例题 8.11、练习 7）。

    **易错点汇总**

    - 正变换是 $e^{-j\omega t}$，逆变换带 $\frac1{2\pi}$ 和 $e^{+j\omega t}$；
    - $\mathcal F[1]=2\pi\delta(\omega)$，$\mathcal F[u(t)]=\frac{1}{j\omega}+\pi\delta(\omega)$；
    - 相似性质的 $\frac1{|a|}$；平移与伸缩同时出现时先提出 $a$；
    - 频域卷积定理的 $\frac1{2\pi}$；
    - $\delta$ 函数只能在积分中通过性质使用，不能逐点取值；
    - 用留数求逆变换时，$t$ 的正负决定在哪个半平面闭合。

---

## 练习题

!!! example "练习 1"

    求三角形脉冲 $f(t)=\begin{cases}1-|t|,&|t|\le1,\\0,&|t|>1\end{cases}$ 的傅里叶变换。

??? success "解答"

    $f$ 为偶函数，$F(\omega)=2\int_0^1(1-t)\cos\omega t\,\mathrm dt$。分部积分：

    $$F(\omega)=2\left[(1-t)\frac{\sin\omega t}{\omega}\bigg|_0^1+\int_0^1\frac{\sin\omega t}{\omega}\,\mathrm dt\right]=\frac{2(1-\cos\omega)}{\omega^2}=\frac{4\sin^2(\omega/2)}{\omega^2}.$$

    **另一种看法。** 三角形脉冲是宽为 $1$ 的矩形脉冲（$|t|\le\frac12$）与自身的卷积。由例题 8.2（$\delta=\frac12$），矩形脉冲的频谱为 $\frac{2\sin(\omega/2)}{\omega}$，由卷积定理，$F=\left(\frac{2\sin(\omega/2)}{\omega}\right)^2$，与上面一致。

!!! example "练习 2"

    求符号函数 $\operatorname{sgn}t=\begin{cases}-1,&t<0,\\1,&t>0\end{cases}$ 的傅里叶变换，并由此验证 $\mathcal F[u(t)]=\frac1{j\omega}+\pi\delta(\omega)$。

??? success "解答"

    $\operatorname{sgn}t=\lim\limits_{\beta\to0^+}\left[e^{-\beta t}u(t)-e^{\beta t}u(-t)\right]$。由例题 8.4、8.5 的计算，

    $$\mathcal F\left[e^{-\beta t}u(t)-e^{\beta t}u(-t)\right]=\frac{1}{\beta+j\omega}-\frac{1}{\beta-j\omega}=\frac{-2j\omega}{\beta^2+\omega^2}\xrightarrow{\beta\to0^+}\frac{-2j}{\omega}=\frac{2}{j\omega}.$$

    所以 $\mathcal F[\operatorname{sgn}t]=\frac{2}{j\omega}$。又 $u(t)=\frac12+\frac12\operatorname{sgn}t$，由线性性质，$\mathcal F[u(t)]=\frac12\cdot2\pi\delta(\omega)+\frac12\cdot\frac{2}{j\omega}=\pi\delta(\omega)+\frac{1}{j\omega}$。

!!! example "练习 3"

    求 $f(t)=\sin\left(5t+\dfrac\pi3\right)$ 的傅里叶变换。

??? success "解答"

    $f(t)=\frac12\sin5t+\frac{\sqrt3}{2}\cos5t$。由例题 8.6，

    $$F(\omega)=\frac{j\pi}{2}\left[\delta(\omega+5)-\delta(\omega-5)\right]+\frac{\sqrt3\pi}{2}\left[\delta(\omega+5)+\delta(\omega-5)\right]=\frac{\pi}{2}\left[(\sqrt3+j)\delta(\omega+5)+(\sqrt3-j)\delta(\omega-5)\right].$$

    （也可以用时移性质：$f(t)=\sin\left[5\left(t+\frac{\pi}{15}\right)\right]$，乘以 $e^{j\omega\pi/15}$，再用 $\delta$ 的筛选性质 $e^{j\omega\pi/15}\delta(\omega\mp5)=e^{\pm j\pi/3}\delta(\omega\mp5)$，结果相同。）

!!! example "练习 4"

    求 $f(t)=e^{-\beta t}u(t)\sin\omega_0t$（$\beta>0$）的傅里叶变换。

??? success "解答"

    由调制定理（$\sin$ 形式）与例题 8.4：

    $$F(\omega)=\frac{1}{2j}\left[\frac{1}{\beta+j(\omega-\omega_0)}-\frac{1}{\beta+j(\omega+\omega_0)}\right]=\frac1{2j}\cdot\frac{2j\omega_0}{(\beta+j\omega)^2+\omega_0^2}=\frac{\omega_0}{(\beta+j\omega)^2+\omega_0^2}.$$

!!! example "练习 5"

    求 $f(t)=te^{-\beta t}u(t)$（$\beta>0$）的傅里叶变换。

??? success "解答"

    由频域微分性质，$\mathcal F[tg(t)]=jG'(\omega)$，$g=e^{-\beta t}u(t)$，$G=\frac1{\beta+j\omega}$：

    $$F(\omega)=j\cdot\frac{-j}{(\beta+j\omega)^2}=\frac{1}{(\beta+j\omega)^2}.$$

    （也可以由卷积定理：$e^{-\beta t}u(t)*e^{-\beta t}u(t)=\int_0^te^{-\beta\tau}e^{-\beta(t-\tau)}\,\mathrm d\tau=te^{-\beta t}$，$t\ge0$，所以它的频谱是 $G^2$。）

!!! example "练习 6"

    求 $f_1(t)=u(t)$ 与 $f_2(t)=e^{-t}u(t)$ 的卷积，并用积分性质验证。

??? success "解答"

    两者都是因果函数，$t<0$ 时卷积为零；$t\ge0$ 时

    $$f_1*f_2(t)=\int_0^tf_2(\tau)f_1(t-\tau)\,\mathrm d\tau=\int_0^te^{-\tau}\,\mathrm d\tau=1-e^{-t}.$$

    即 $u(t)*e^{-t}u(t)=(1-e^{-t})u(t)$。一般地 $f*u=\int_{-\infty}^tf(\tau)\,\mathrm d\tau$，与 $u$ 卷积就是**求变上限积分**。频域验证：$\left(\frac{1}{j\omega}+\pi\delta(\omega)\right)\frac{1}{1+j\omega}=\frac{1}{j\omega(1+j\omega)}+\pi\delta(\omega)=\frac{1}{j\omega}+\pi\delta(\omega)-\frac{1}{1+j\omega}$，逆变换正是 $u(t)-e^{-t}u(t)$。

!!! example "练习 7"

    用帕塞瓦尔等式计算 $\displaystyle\int_{-\infty}^{+\infty}\frac{\mathrm d\omega}{(1+\omega^2)^2}$。

??? success "解答"

    取 $f(t)=e^{-|t|}$，$F(\omega)=\frac{2}{1+\omega^2}$（例题 8.5）。$\int_{-\infty}^{+\infty}f^2\,\mathrm dt=2\int_0^{+\infty}e^{-2t}\,\mathrm dt=1$，所以

    $$\frac1{2\pi}\int_{-\infty}^{+\infty}\frac{4}{(1+\omega^2)^2}\,\mathrm d\omega=1,\qquad\int_{-\infty}^{+\infty}\frac{\mathrm d\omega}{(1+\omega^2)^2}=\frac\pi2 .$$

    与例题 5.11(3) 用留数得到的结果相同。

!!! example "练习 8"

    求 $F(\omega)=\dfrac{1}{\omega^2+4\omega+5}$ 的傅里叶逆变换：(1) 用留数；(2) 用频移性质验证。

??? success "解答"

    (1) $\omega^2+4\omega+5=(\omega+2-j)(\omega+2+j)$，极点 $-2\pm j$，分母比分子高两次。

    $t>0$：上半平面闭合，只有 $\omega=-2+j$：

    $$f(t)=\frac{1}{2\pi}\cdot2\pi j\cdot\frac{e^{j(-2+j)t}}{2j}=\frac12e^{-t}e^{-j2t}.$$

    $t<0$：下半平面闭合（顺时针），只有 $\omega=-2-j$：

    $$f(t)=\frac{1}{2\pi}\cdot(-2\pi j)\cdot\frac{e^{j(-2-j)t}}{-2j}=\frac12e^{t}e^{-j2t}.$$

    合起来 $f(t)=\frac12e^{-|t|}e^{-j2t}$。

    (2) $F(\omega)=\frac{1}{1+(\omega+2)^2}$。由例题 8.5，$\mathcal F^{-1}\left[\frac{1}{1+\omega^2}\right]=\frac12e^{-|t|}$；$F(\omega+2)$ 对应乘以 $e^{-j2t}$（频移性质，$\omega_0=-2$），得 $f(t)=\frac12e^{-|t|}e^{-j2t}$，一致。

---

**下一章预告。** 傅里叶变换要求信号绝对可积，连 $u(t)$、$\sin t$ 都要借助 $\delta$ 函数，$e^{t}$ 这样增长的信号则根本没有傅里叶变换。第 9 章给信号乘上衰减因子 $e^{-\sigma t}$，只看 $t\ge0$ 的部分，再做傅里叶变换，就得到**拉普拉斯变换** $F(s)=\int_0^{+\infty}f(t)e^{-st}\,\mathrm dt$，其中 $s=\sigma+j\omega$。它在一个右半平面内解析，反演要用留数，正好用上第 2–5 章的全部工具。→ [第 9 章 拉普拉斯变换](ch9.md)
