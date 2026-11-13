# 第 3 章 复变函数的积分

!!! tip "本章目标"

    1. 理解复积分 $\int_C f(z)\,\mathrm dz$ 的定义（黎曼和的极限），会用**参数方程**和**实部虚部**两种方法计算，会用**长大不等式**估计积分；
    2. 牢记最重要的积分 $\displaystyle\oint_{|z-z_0|=r}\frac{\mathrm dz}{(z-z_0)^{n}}$；
    3. 掌握**柯西积分定理**、闭路变形原理、复合闭路定理，以及由此得到的原函数与牛顿-莱布尼茨公式；
    4. 熟练运用**柯西积分公式**与**高阶导数公式**计算闭路积分——这是第 5 章留数法的前身，也是期末必考；
    5. 了解平均值公式、最大模原理、柯西不等式、刘维尔定理，知道“解析”这一条件有多强。

**本章主线。** 第 2 章说解析函数“局部是旋转 + 伸缩”。本章从积分角度看同一件事：在没有奇点的地方，闭路积分为零（柯西定理）；回路的形状无关紧要，**只有它围住了哪些奇点**才要紧（闭路变形）；被围住的奇点又可以用一个公式精确算出贡献（柯西积分公式）。顺带得到一个惊人的结论：解析函数**无穷次可导**。

---

## 3.1 复积分的概念

### 3.1.1 复积分的定义

!!! definition "定义 3.1（复积分）"

    设 $C$ 是复平面上以 $A$ 为起点、$B$ 为终点的光滑（或分段光滑）有向曲线，函数 $f(z)$ 在 $C$ 上有定义。在 $C$ 上依次取分点

    $$A=z_0,\ z_1,\ \dots,\ z_{n-1},\ z_n=B,$$

    在每段弧 $\widehat{z_{k-1}z_k}$ 上任取一点 $\zeta_k$，记 $\Delta z_k=z_k-z_{k-1}$，作和式 $S_n=\sum_{k=1}^n f(\zeta_k)\,\Delta z_k$。记 $\delta=\max_k\lvert\Delta s_k\rvert$（$\Delta s_k$ 为弧长）。若当 $\delta\to0$ 时，无论如何分割、如何取 $\zeta_k$，$S_n$ 都趋于同一个复数 $I$，则称 $f$ 沿 $C$ **可积**，记

    $$\int_C f(z)\,\mathrm dz=\lim_{\delta\to0}\sum_{k=1}^n f(\zeta_k)\,\Delta z_k .$$

    当 $C$ 为闭曲线时记作 $\oint_C f(z)\,\mathrm dz$，约定**正向为逆时针**（区域在左手边）。

和实积分 $\int_a^b f(x)\,\mathrm dx$ 的区别：$\Delta z_k$ 是**有方向的小复数**，$f(\zeta_k)\Delta z_k$ 是把 $\Delta z_k$ 旋转、伸缩后的小向量；积分就是这些小向量**首尾相接**后的总位移。它更像物理里的**第二类曲线积分**（做功），而不是面积。

<figure class="video" markdown>
<video controls preload="metadata" width="100%" poster="assets/posters/C3_RiemannSum.jpg">
<source src="assets/videos/C3_RiemannSum.mp4" type="video/mp4">
</video>
<figcaption>动画 3-1　复积分是黎曼和的极限：$z^2$ 沿四分之一单位圆从 $1$ 到 $i$，小向量 $f(\zeta_k)\Delta z_k$ 首尾相接，分点越多，终点越接近 $\frac{-1-i}{3}$。</figcaption>
</figure>

### 3.1.2 存在性与计算方法

!!! theorem "定理 3.1（存在性与化为实积分）"

    若 $f(z)=u(x,y)+iv(x,y)$ 在分段光滑曲线 $C$ 上连续，则 $\int_C f(z)\,\mathrm dz$ 存在，且

    $$\int_C f(z)\,\mathrm dz=\int_C u\,\mathrm dx-v\,\mathrm dy+i\int_C v\,\mathrm dx+u\,\mathrm dy .$$

    若 $C$ 的参数方程为 $z=z(t)=x(t)+iy(t)$，$t$ 从 $\alpha$ 变到 $\beta$（$t=\alpha$ 对应起点），则

    $$\int_C f(z)\,\mathrm dz=\int_\alpha^\beta f\big(z(t)\big)\,z'(t)\,\mathrm dt .$$

**记忆方法。** 把 $f\,\mathrm dz=(u+iv)(\mathrm dx+i\,\mathrm dy)$ 形式地展开即得第一式；把 $z=z(t)$、$\mathrm dz=z'(t)\,\mathrm dt$ 代入即得第二式。实际计算几乎总用**参数法**。

!!! success "常用参数方程"

    | 曲线 | 参数方程 |
    |---|---|
    | 从 $z_1$ 到 $z_2$ 的直线段 | $z=z_1+t(z_2-z_1),\ t:0\to1$ |
    | 圆周 $\lvert z-z_0\rvert=r$（正向） | $z=z_0+re^{i\theta},\ \theta:0\to2\pi$ |
    | 曲线 $y=g(x)$，$x$ 从 $a$ 到 $b$ | $z=x+ig(x),\ x:a\to b$ |

!!! example "例题 3.1（教材例 3.1、例 3.3 改编）"

    设 $C_1$：从 $0$ 到 $1$ 的线段；$C_2$：从 $1$ 到 $1+i$ 的线段；$C_3$：从 $0$ 到 $1+i$ 的线段；$C_4$：抛物线 $x=y^2$ 上从 $0$ 到 $1+i$ 的一段。分别计算 $\int_C z\,\mathrm dz$ 与 $\int_C\bar z\,\mathrm dz$，其中 $C=C_1+C_2,\ C_3,\ C_4$。

    **解。** 参数方程：$C_1:z=x,\ x:0\to1$；$C_2:z=1+iy,\ y:0\to1$；$C_3:z=(1+i)t,\ t:0\to1$；$C_4:z=t^2+it,\ t:0\to1$。

    **(a) $f=z$。**

    $$\int_{C_3}z\,\mathrm dz=\int_0^1(1+i)t\cdot(1+i)\,\mathrm dt=\frac{(1+i)^2}{2}=i,$$

    $$\int_{C_1+C_2}z\,\mathrm dz=\int_0^1x\,\mathrm dx+\int_0^1(1+iy)\,i\,\mathrm dy=\frac12+\left(i-\frac12\right)=i .$$

    $C_4$ 上同样得 $i$。三条路径结果相同，都等于 $\frac{z^2}{2}\Big|_0^{1+i}=i$。

    **(b) $f=\bar z$。**

    $$\int_{C_1+C_2}\bar z\,\mathrm dz=\int_0^1x\,\mathrm dx+\int_0^1(1-iy)\,i\,\mathrm dy=\frac12+\left(i+\frac12\right)=1+i,$$

    $$\int_{C_3}\bar z\,\mathrm dz=\int_0^1(1-i)t\,(1+i)\,\mathrm dt=2\cdot\frac12=1,$$

    $$\int_{C_4}\bar z\,\mathrm dz=\int_0^1(t^2-it)(2t+i)\,\mathrm dt=\int_0^1(2t^3+t-it^2)\,\mathrm dt=1-\frac i3 .$$

    三条路径结果各不相同。

    **结论。** $z$ 解析，积分与路径无关；$\bar z$ 处处不解析，积分依赖路径。这正是 3.2 节柯西定理要回答的问题。

### 3.1.3 复积分的性质

!!! success "性质（设所涉及的积分都存在）"

    1. 线性：$\int_C[\alpha f+\beta g]\,\mathrm dz=\alpha\int_Cf\,\mathrm dz+\beta\int_Cg\,\mathrm dz$；
    2. 反向变号：$\int_{C^-}f\,\mathrm dz=-\int_Cf\,\mathrm dz$；
    3. 路径可加：$C=C_1+C_2\Rightarrow\int_C=\int_{C_1}+\int_{C_2}$；
    4. **估值不等式**：

        $$\left|\int_C f(z)\,\mathrm dz\right|\le\int_C|f(z)|\,\mathrm ds\le ML,$$

        其中 $M=\max_{z\in C}|f(z)|$，$L$ 为 $C$ 的长度。第二个不等式俗称**长大不等式**（“长度乘最大值”）。

!!! warning "易错点：$\mathrm dz$ 与 $\mathrm ds$"

    $|\mathrm dz|=\mathrm ds$ 是弧长元素（非负实数），而 $\mathrm dz$ 是复数。$\int_C|f|\,\mathrm ds$ 是实积分；不要写成 $\int_C|f|\,\mathrm dz$。

### 3.1.4 最重要的一个积分

!!! example "例题 3.2（教材例 3.2）"

    设 $n$ 为整数，计算 $\displaystyle I_n=\oint_{|z-z_0|=r}\frac{\mathrm dz}{(z-z_0)^{n}}$（$r>0$）。

    **解。** 令 $z=z_0+re^{i\theta}$，$\theta:0\to2\pi$，$\mathrm dz=ire^{i\theta}\mathrm d\theta$：

    $$I_n=\int_0^{2\pi}\frac{ire^{i\theta}}{r^ne^{in\theta}}\,\mathrm d\theta=\frac{i}{r^{n-1}}\int_0^{2\pi}e^{-i(n-1)\theta}\,\mathrm d\theta .$$

    若 $n=1$，被积函数为 $1$，$I_1=2\pi i$。若 $n\ne1$，$\int_0^{2\pi}\big[\cos(n-1)\theta-i\sin(n-1)\theta\big]\mathrm d\theta=0$，$I_n=0$。

!!! theorem "必背结论"

    $$\oint_{|z-z_0|=r}\frac{\mathrm dz}{(z-z_0)^{n}}=\begin{cases}2\pi i,& n=1,\\ 0,& n\ne1\ (n\in\mathbb Z),\end{cases}$$

    结果与圆心 $z_0$ 和半径 $r$ 都无关。

**为什么只有 $n=1$ 不为零？** 在圆上，$\dfrac{\mathrm dz}{z-z_0}=i\,\mathrm d\theta$：每一小段的贡献都是同一个方向 $i$，加起来 $2\pi i$。而当 $n=2$ 时，贡献 $ie^{-i\theta}\mathrm d\theta/r$ 的方向转了一整圈，首尾相接回到原点。$n=1$ 的这个 $2\pi i$ 就是第 5 章“留数”的来源。

<figure class="video" markdown>
<video controls preload="metadata" width="100%" poster="assets/posters/C3_OneOverZ.jpg">
<source src="assets/videos/C3_OneOverZ.mp4" type="video/mp4">
</video>
<figcaption>动画 3-2　$\oint\frac{\mathrm dz}{z}$：左边是小弦 $\Delta z_k$，右边是累加的 $\frac{\Delta z_k}{z_k}$，全都指向 $i$，总和 $2\pi i$；$\oint\frac{\mathrm dz}{z^2}$ 的小向量绕一圈，总和为 $0$。</figcaption>
</figure>

!!! example "例题 3.3（长大不等式）"

    设 $C$ 为上半单位圆周，从 $1$ 到 $-1$。证明 $\left|\displaystyle\int_C\frac{e^z}{z}\,\mathrm dz\right|\le e\pi$。

    **证。** 在 $C$ 上 $|z|=1$，$\left|\dfrac{e^z}{z}\right|=e^x\le e^1$（因为 $x\le1$），$C$ 的长度为 $\pi$，所以积分的模 $\le e\pi$。

!!! example "例题 3.4（教材例 3.4）"

    设 $C$ 为从 $0$ 到 $3+4i$ 的直线段，估计 $\left|\displaystyle\int_C\frac{\mathrm dz}{z-i}\right|$ 的一个上界。

    **解。** $C:z=(3+4i)t,\ t\in[0,1]$，长度 $L=5$。

    $$|z-i|^2=|3t+i(4t-1)|^2=25t^2-8t+1=25\left(t-\tfrac{4}{25}\right)^2+\tfrac{9}{25}\ge\tfrac{9}{25},$$

    故 $\left|\frac1{z-i}\right|\le\frac53$，从而

    $$\left|\int_C\frac{\mathrm dz}{z-i}\right|\le\frac53\cdot5=\frac{25}{3}.$$

    （从几何上看，$\frac35$ 就是点 $i$ 到这条线段的距离。）

!!! example "例题 3.5（教材例 3.5）"

    证明 $\displaystyle\lim_{r\to0}\oint_{|z|=r}\frac{z^3}{1+z^2}\,\mathrm dz=0$。

    **证。** 不妨设 $r<1$。在 $|z|=r$ 上，$|1+z^2|\ge1-|z|^2=1-r^2$，所以

    $$\left|\oint_{|z|=r}\frac{z^3}{1+z^2}\,\mathrm dz\right|\le\frac{r^3}{1-r^2}\cdot2\pi r\to0\quad(r\to0).$$

    注意 $|1+z^2|\ge\big|1-|z|^2\big|$ 用的是三角不等式 $|a+b|\ge\big||a|-|b|\big|$（第 1 章）。

---

## 3.2 柯西积分定理

### 3.2.1 柯西-古萨基本定理

例题 3.1 中 $z$ 的积分与路径无关、$\bar z$ 的积分与路径有关。等价地说，$\oint z\,\mathrm dz=0$，而 $\oint\bar z\,\mathrm dz\ne0$。什么时候闭路积分为零？

!!! theorem "定理 3.2（柯西-古萨基本定理）"

    若函数 $f(z)$ 在单连通区域 $D$ 内处处解析，则对 $D$ 内**任意一条**闭曲线 $C$（不必是简单曲线），

    $$\oint_C f(z)\,\mathrm dz=0 .$$

??? proof "证明思路（附加条件 $f'$ 连续）"

    柯西最初的证明假设 $f'(z)$ 连续，这时 $u,v$ 的偏导数连续，可用格林公式。设 $C$ 为简单闭曲线，围成区域 $G\subset D$：

    $$\oint_C u\,\mathrm dx-v\,\mathrm dy=\iint_G\left(-v_x-u_y\right)\mathrm dx\,\mathrm dy,\qquad
    \oint_C v\,\mathrm dx+u\,\mathrm dy=\iint_G\left(u_x-v_y\right)\mathrm dx\,\mathrm dy .$$

    由 C-R 方程 $u_x=v_y,\ u_y=-v_x$，两个二重积分的被积函数都为零，故 $\oint_Cf\,\mathrm dz=0$。古萨（Goursat）证明了“$f'$ 连续”这一条件可以去掉，证明见教材，这里从略。$\blacksquare$

!!! theorem "推论"

    设 $C$ 是区域 $D$ 的边界，$f(z)$ 在 $D$ 内解析、在闭区域 $\bar D=D\cup C$ 上连续，则 $\oint_C f(z)\,\mathrm dz=0$。

!!! warning "易错点：三个条件缺一不可"

    - **解析**：$\oint_{|z|=1}\bar z\,\mathrm dz=\int_0^{2\pi}e^{-i\theta}\cdot ie^{i\theta}\mathrm d\theta=2\pi i\ne0$；
    - **单连通**（或者说，$C$ 内部没有奇点）：$\frac1z$ 在 $0<|z|<2$ 内解析，但 $\oint_{|z|=1}\frac{\mathrm dz}{z}=2\pi i$，因为 $C$ 内部有“洞” $z=0$；
    - **闭曲线**：沿非闭曲线的积分一般不为零。

    反过来，$\oint_Cf\,\mathrm dz=0$ **不能**推出 $f$ 解析。例如 $\oint_{|z|=1}\frac{\mathrm dz}{z^2}=0$，但 $\frac{1}{z^2}$ 在 $z=0$ 不解析。

!!! example "例题 3.6"

    计算 (1) $\displaystyle\oint_{|z|=1}\frac{e^z}{z-2}\,\mathrm dz$；(2) $\displaystyle\oint_{|z|=1}\frac{\mathrm dz}{\cos z}$。

    **解。** (1) 奇点 $z=2$ 在 $|z|=1$ 外，被积函数在 $|z|\le1$ 上解析，积分为 $0$。

    (2) $\cos z=0\iff z=\frac\pi2+k\pi$（第 2 章），离原点最近的零点 $\pm\frac\pi2$ 也满足 $\frac\pi2>1$，都在单位圆外，所以积分为 $0$。

### 3.2.2 闭路变形原理

!!! theorem "定理 3.3（闭路变形原理）"

    设 $C$ 与 $C_1$ 是两条简单闭曲线，$C_1$ 在 $C$ 的内部，$D$ 是它们之间的二连通区域。若 $f(z)$ 在 $D$ 内及边界 $C,C_1$ 上解析，则

    $$\oint_C f(z)\,\mathrm dz=\oint_{C_1}f(z)\,\mathrm dz\qquad(C,\ C_1\ 均取逆时针方向).$$

??? proof "证明"

    作两条割线 $AA'$、$BB'$ 连接 $C$ 与 $C_1$，把 $D$ 分成两个单连通区域。分别对两块的边界用柯西定理，两式相加：割线各被来回走了一次，积分抵消，剩下

    $$\oint_C f\,\mathrm dz+\oint_{C_1^-}f\,\mathrm dz=0,$$

    即 $\oint_Cf\,\mathrm dz=\oint_{C_1}f\,\mathrm dz$。$\blacksquare$

**含义。** 只要变形过程中**不越过奇点**，闭路积分的值不变。所以计算时总可以把任意形状的回路换成绕奇点的小圆。结合例题 3.2：

$$\oint_C\frac{\mathrm dz}{z-z_0}=2\pi i\quad(\text{$C$ 为任意包围 $z_0$ 的正向简单闭曲线}).$$

### 3.2.3 复合闭路定理

!!! theorem "定理 3.4（复合闭路定理）"

    设 $C$ 是简单闭曲线，$C_1,C_2,\dots,C_n$ 是在 $C$ 内部的简单闭曲线，它们互不包含也互不相交，以 $C,C_1,\dots,C_n$ 为边界的多连通区域记为 $D$。若 $f(z)$ 在 $\bar D$ 上解析，则

    1. $\displaystyle\oint_\Gamma f(z)\,\mathrm dz=0$，其中 $\Gamma=C+C_1^-+\cdots+C_n^-$ 是 $D$ 的**正向边界**（外边界逆时针，内边界顺时针）；
    2. $\displaystyle\oint_C f(z)\,\mathrm dz=\sum_{k=1}^n\oint_{C_k}f(z)\,\mathrm dz$（全部取逆时针方向）。

<figure class="video" markdown>
<video controls preload="metadata" width="100%" poster="assets/posters/C3_Deformation.jpg">
<source src="assets/videos/C3_Deformation.mp4" type="video/mp4">
</video>
<figcaption>动画 3-3　闭路变形：任意回路收缩成绕奇点的小圆，积分始终是 $2\pi i$；移开到不包围奇点，积分立刻变成 $0$。复合闭路：用割缝把多连通域变成单连通域。</figcaption>
</figure>

!!! example "例题 3.7（教材例 3.7 改编）"

    计算 $\displaystyle I=\oint_C\frac{2z-1}{z^2-z}\,\mathrm dz$，其中 (1) $C:|z|=3$；(2) $C:\left|z-1\right|=\frac12$；(3) $C:|z+2|=1$。

    **解。** 奇点为 $z=0$ 和 $z=1$。部分分式：$\dfrac{2z-1}{z(z-1)}=\dfrac1z+\dfrac1{z-1}$。

    (1) 两个奇点都在 $|z|=3$ 内。由复合闭路定理及例题 3.2，

    $$I=\oint_C\frac{\mathrm dz}{z}+\oint_C\frac{\mathrm dz}{z-1}=2\pi i+2\pi i=4\pi i .$$

    (2) 只有 $z=1$ 在内。$\frac1z$ 在 $|z-1|\le\frac12$ 上解析，积分为 $0$，所以 $I=2\pi i$。

    (3) 两个奇点都不在 $|z+2|\le1$ 内，$I=0$。

### 3.2.4 路径无关性与原函数

!!! theorem "定理 3.5（路径无关）"

    若 $f(z)$ 在单连通区域 $D$ 内解析，则 $\int_C f(z)\,\mathrm dz$ 只与 $C$ 的起点 $z_0$ 和终点 $z_1$ 有关，记作 $\int_{z_0}^{z_1}f(z)\,\mathrm dz$。

（两条从 $z_0$ 到 $z_1$ 的路径 $C_1,C_2$ 拼成闭路 $C_1+C_2^-$，由柯西定理得 $\int_{C_1}=\int_{C_2}$。）

!!! definition "定义 3.2（原函数）"

    设在区域 $D$ 内 $F'(z)=f(z)$，则称 $F(z)$ 为 $f(z)$ 在 $D$ 内的一个**原函数**。$f$ 的全体原函数 $F(z)+c$ 称为 $f$ 的**不定积分**。

!!! theorem "定理 3.6（变上限积分与牛顿-莱布尼茨公式）"

    设 $f(z)$ 在单连通区域 $D$ 内解析，$z_0\in D$。

    1. $F(z)=\displaystyle\int_{z_0}^{z}f(\zeta)\,\mathrm d\zeta$ 在 $D$ 内解析，且 $F'(z)=f(z)$；
    2. $f$ 的任意两个原函数相差一个常数；
    3. 若 $G(z)$ 是 $f$ 的一个原函数，则

        $$\int_{z_0}^{z_1}f(z)\,\mathrm dz=G(z_1)-G(z_0)\qquad(z_0,z_1\in D).$$

??? proof "1 的证明"

    取 $z+\Delta z\in D$，由路径无关可沿线段从 $z$ 积到 $z+\Delta z$：

    $$\frac{F(z+\Delta z)-F(z)}{\Delta z}-f(z)=\frac{1}{\Delta z}\int_z^{z+\Delta z}\big[f(\zeta)-f(z)\big]\,\mathrm d\zeta .$$

    由 $f$ 连续，对任意 $\varepsilon>0$，当 $|\Delta z|$ 足够小时线段上 $|f(\zeta)-f(z)|<\varepsilon$，由长大不等式，上式的模 $\le\frac{1}{|\Delta z|}\cdot\varepsilon|\Delta z|=\varepsilon$。故 $F'(z)=f(z)$。$\blacksquare$

!!! example "例题 3.8"

    计算 (1) $\displaystyle\int_0^{1+i}z^2\,\mathrm dz$；(2) $\displaystyle\int_0^{i}z\cos z\,\mathrm dz$；(3) $\displaystyle\int_C\sin z\,\mathrm dz$，其中 $C$ 是上半圆周 $|z-1|=1$，从 $0$ 到 $2$（教材例 3.6）。

    **解。** 被积函数都在全平面解析，可以用牛顿-莱布尼茨公式。

    (1) $\displaystyle\int_0^{1+i}z^2\,\mathrm dz=\frac{z^3}{3}\Big|_0^{1+i}=\frac{(1+i)^3}{3}=\frac{-2+2i}{3}$。

    (2) 分部积分同样适用：

    $$\int_0^iz\cos z\,\mathrm dz=\big[z\sin z+\cos z\big]_0^i=i\sin i+\cos i-1 .$$

    由 $\sin i=i\sinh1$，$\cos i=\cosh1$，结果为 $-\sinh1+\cosh1-1=e^{-1}-1$。

    (3) $\displaystyle\int_C\sin z\,\mathrm dz=-\cos z\Big|_0^2=1-\cos2$。（也可以把路径换成实轴上的线段 $[0,2]$。）

!!! warning "易错点：牛顿-莱布尼茨公式的前提"

    需要 $f$ 在一个**包含积分路径的单连通区域**内解析。

    $\frac1z$ 在 $\mathbb C\setminus\{0\}$ 内解析，但这个区域不是单连通的。$\ln z$ 只是 $\frac1z$ 在“割去负实轴”的区域内的原函数，不能拿 $\ln z$ 去算绕原点一周的积分（否则会得到 $0$，正确答案是 $2\pi i$）。

    $\bar z$、$\operatorname{Re}z$、$|z|$ 等不解析的函数没有原函数，只能用参数法。

---

## 3.3 柯西积分公式

### 3.3.1 柯西积分公式

!!! theorem "定理 3.7（柯西积分公式）"

    设 $f(z)$ 在区域 $D$ 内解析，$C$ 为 $D$ 内的一条正向简单闭曲线，其内部完全属于 $D$，$z_0$ 为 $C$ 内任意一点，则

    $$f(z_0)=\frac{1}{2\pi i}\oint_C\frac{f(z)}{z-z_0}\,\mathrm dz .$$

??? proof "证明"

    $\frac{f(z)}{z-z_0}$ 在 $C$ 内只有奇点 $z_0$。由闭路变形原理，可把 $C$ 换成小圆 $K:|z-z_0|=\rho$：

    $$\oint_C\frac{f(z)}{z-z_0}\,\mathrm dz=\oint_K\frac{f(z)}{z-z_0}\,\mathrm dz=\underbrace{\oint_K\frac{f(z_0)}{z-z_0}\,\mathrm dz}_{=2\pi i f(z_0)}+\oint_K\frac{f(z)-f(z_0)}{z-z_0}\,\mathrm dz .$$

    由 $f$ 在 $z_0$ 连续，任给 $\varepsilon>0$，当 $\rho$ 足够小时 $K$ 上 $|f(z)-f(z_0)|<\varepsilon$，所以最后一项的模 $\le\frac{\varepsilon}{\rho}\cdot2\pi\rho=2\pi\varepsilon$。但左边与 $\rho$ 无关，故最后一项只能为 $0$。$\blacksquare$

**含义。**

1. 解析函数在 $C$ **内部**的值完全由它在**边界** $C$ 上的值决定。这和静电学“边界上的电势决定内部电势”是同一类事实（第 2 章：解析函数的实部是调和函数）。
2. 反过来读就是一个**计算积分**的公式：

    $$\oint_C\frac{f(z)}{z-z_0}\,\mathrm dz=2\pi i\,f(z_0)\qquad(f\ \text{在 $C$ 内及 $C$ 上解析，$z_0$ 在 $C$ 内}).$$

!!! success "用柯西积分公式计算闭路积分的步骤"

    1. 找出被积函数在 $C$ **内部**的奇点；
    2. 若只有一个奇点 $z_0$，把被积函数写成 $\dfrac{f(z)}{z-z_0}$，其中 $f$ 在 $C$ 内解析（把“外面”的因子都并进 $f$）；
    3. 若有多个奇点，用复合闭路定理在每个奇点周围取小圆分别计算，或者先做部分分式；
    4. 答案为 $2\pi i\,f(z_0)$ 之和。

!!! example "例题 3.9"

    计算 (1) $\displaystyle\oint_{|z|=2}\frac{e^z}{z-1}\,\mathrm dz$；(2) $\displaystyle\oint_{|z|=2}\frac{z}{(9-z^2)(z+i)}\,\mathrm dz$（教材例 3.10）。

    **解。** (1) $e^z$ 解析，$z=1$ 在圆内：$I=2\pi i\,e^1=2\pi e\,i$。

    (2) 奇点 $\pm3$ 在圆外，$-i$ 在圆内。令 $f(z)=\dfrac{z}{9-z^2}$，它在 $|z|\le2$ 上解析：

    $$I=\oint_{|z|=2}\frac{f(z)}{z-(-i)}\,\mathrm dz=2\pi i\cdot f(-i)=2\pi i\cdot\frac{-i}{9-(-i)^2}=2\pi i\cdot\frac{-i}{10}=\frac{\pi}{5}.$$

!!! example "例题 3.10"

    计算 $\displaystyle\oint_{|z|=2}\frac{\mathrm dz}{z^2+1}$。

    **解。** 奇点 $\pm i$ 都在圆内。以 $\pm i$ 为心作互不相交的小圆 $C_1,C_2$，由复合闭路定理：

    $$I=\oint_{C_1}\frac{\frac{1}{z+i}}{z-i}\,\mathrm dz+\oint_{C_2}\frac{\frac{1}{z-i}}{z+i}\,\mathrm dz=2\pi i\cdot\frac{1}{2i}+2\pi i\cdot\frac{1}{-2i}=0 .$$

    或用部分分式：$\dfrac1{z^2+1}=\dfrac1{2i}\left(\dfrac1{z-i}-\dfrac1{z+i}\right)$，两项各给出 $\frac{2\pi i}{2i}$，相减为 $0$。

!!! warning "易错点"

    - $z_0$ 必须在 $C$ **内部**；在外部时积分为 $0$（柯西定理），在 $C$ 上时积分本身没有定义；
    - 公式中的 $f$ 必须在 $C$ 内**解析**。例题 3.9(2) 中取 $f=\frac{z}{9-z^2}$ 之所以可行，是因为 $\pm3$ 在 $|z|=2$ 外；若换成 $|z|=4$，$f$ 在圆内不再解析，必须用复合闭路定理（练习 6）；
    - 结果是 $2\pi i\,f(z_0)$，不要漏掉 $2\pi i$，也不要把 $f$ 和被积函数 $\frac{f}{z-z_0}$ 搞混。

### 3.3.2 平均值公式

!!! theorem "定理 3.8（平均值公式）"

    若 $f(z)$ 在 $|z-z_0|<R$ 内解析，在 $|z-z_0|\le R$ 上连续，则

    $$f(z_0)=\frac{1}{2\pi}\int_0^{2\pi}f\big(z_0+Re^{i\theta}\big)\,\mathrm d\theta .$$

（在柯西积分公式中令 $z=z_0+Re^{i\theta}$ 即得。）解析函数在圆心的值等于它在圆周上的**平均值**。取实部：调和函数也有同样的平均值性质。

### 3.3.3 最大模原理

!!! theorem "定理 3.9（最大模原理）"

    若 $f(z)$ 在区域 $D$ 内解析且不为常数，则 $|f(z)|$ 在 $D$ 内**没有最大值**。

    **推论 1** 在区域 $D$ 内解析的函数，如果其模在 $D$ 内某点达到最大值，则它恒为常数。

    **推论 2** 若 $f$ 在有界区域 $D$ 内解析、在 $\bar D$ 上连续，则 $|f(z)|$ 的最大值在 $D$ 的**边界**上取得。

（直观理解：由平均值公式，$|f(z_0)|\le$ 周围圆周上 $|f|$ 的平均值，所以 $|f|$ 不可能在内部“冒尖”。）

!!! example "例题 3.11（教材例 3.11）"

    设 $f(z)$ 在 $|z|<R$ 内解析，$0<r<R$，令 $M(r)=\max_{|z|=r}|f(z)|$。证明 $M(r)$ 是 $r$ 的单调递增函数。

    **证。** 由推论 2，$|f|$ 在 $|z|\le r$ 上的最大值在边界 $|z|=r$ 上取得，即 $M(r)=\max_{|z|\le r}|f(z)|$。当 $r_1<r_2$ 时，$\{|z|\le r_1\}\subset\{|z|\le r_2\}$，故 $M(r_1)\le M(r_2)$。

---

## 3.4 解析函数的高阶导数

### 3.4.1 高阶导数公式

对柯西积分公式 $f(z)=\frac{1}{2\pi i}\oint_C\frac{f(\zeta)}{\zeta-z}\,\mathrm d\zeta$ 关于 $z$ 求导，形式上只作用在 $\frac{1}{\zeta-z}$ 上：$\frac{\mathrm d^n}{\mathrm dz^n}\frac1{\zeta-z}=\frac{n!}{(\zeta-z)^{n+1}}$。这可以严格证明。

!!! theorem "定理 3.10（高阶导数公式）"

    解析函数 $f(z)$ 的导数仍为解析函数，因而 $f$ 具有**任意阶导数**，且

    $$f^{(n)}(z_0)=\frac{n!}{2\pi i}\oint_C\frac{f(z)}{(z-z_0)^{n+1}}\,\mathrm dz\qquad(n=1,2,\dots),$$

    其中 $C$ 为 $f$ 的解析区域 $D$ 内围绕 $z_0$ 的任意正向简单闭曲线，其内部完全属于 $D$。

??? proof "证明（$n=1$）"

    由柯西积分公式，

    $$\frac{f(z_0+h)-f(z_0)}{h}=\frac{1}{2\pi i}\oint_C\frac{f(z)}{(z-z_0-h)(z-z_0)}\,\mathrm dz .$$

    它与 $\frac{1}{2\pi i}\oint_C\frac{f(z)}{(z-z_0)^2}\,\mathrm dz$ 的差为

    $$\frac{h}{2\pi i}\oint_C\frac{f(z)}{(z-z_0)^2(z-z_0-h)}\,\mathrm dz .$$

    设 $|f|\le M$（在 $C$ 上），$d$ 为 $z_0$ 到 $C$ 的距离，$L$ 为 $C$ 的长度。当 $|h|<\frac d2$ 时 $|z-z_0-h|\ge\frac d2$，差的模 $\le\frac{|h|}{2\pi}\cdot\frac{ML}{d^2\cdot d/2}\to0$。一般的 $n$ 用归纳法。$\blacksquare$

**这一点和实函数完全不同。** 实函数可导一次，导数可能不连续，更谈不上二阶导。而复函数只要在区域内**可导一次**（即解析），就自动无穷次可导。原因仍是那个方向条件：它比实可导苛刻得多。

!!! success "计算公式（把高阶导数公式反过来用）"

    $$\oint_C\frac{f(z)}{(z-z_0)^{n+1}}\,\mathrm dz=\frac{2\pi i}{n!}\,f^{(n)}(z_0)$$

    （$f$ 在 $C$ 内及 $C$ 上解析，$z_0$ 在 $C$ 内）。注意分母指数是 $n+1$，对应的是 $n$ 阶导数，阶乘也是 $n!$。

!!! example "例题 3.12（教材例 3.12 部分）"

    计算 (1) $\displaystyle\oint_{|z-i|=1}\frac{\cos z}{(z-i)^3}\,\mathrm dz$；(2) $\displaystyle\oint_{|z|=1}\frac{e^z}{z^{100}}\,\mathrm dz$。

    **解。** (1) $n+1=3$，$n=2$，$f=\cos z$：

    $$I=\frac{2\pi i}{2!}(\cos z)''\Big|_{z=i}=\pi i\,(-\cos i)=-\pi i\cosh1 .$$

    (2) $n=99$，$(e^z)^{(99)}=e^z$：$I=\dfrac{2\pi i}{99!}e^0=\dfrac{2\pi i}{99!}$。

!!! example "例题 3.13"

    计算 $\displaystyle I=\oint_{|z|=2}\frac{e^z}{(z^2+1)^2}\,\mathrm dz$。

    **解。** $(z^2+1)^2=(z-i)^2(z+i)^2$，二阶奇点 $\pm i$ 都在圆内。分别以 $\pm i$ 为心作小圆 $C_1,C_2$。

    在 $C_1$ 上取 $f_1(z)=\dfrac{e^z}{(z+i)^2}$，$f_1'(z)=\dfrac{e^z(z+i-2)}{(z+i)^3}$：

    $$\oint_{C_1}=2\pi i\,f_1'(i)=2\pi i\cdot\frac{e^i(2i-2)}{(2i)^3}=\frac{\pi(1-i)e^{i}}{2}.$$

    在 $C_2$ 上取 $f_2(z)=\dfrac{e^z}{(z-i)^2}$，$f_2'(z)=\dfrac{e^z(z-i-2)}{(z-i)^3}$：

    $$\oint_{C_2}=2\pi i\,f_2'(-i)=2\pi i\cdot\frac{e^{-i}(-2i-2)}{(-2i)^3}=-\frac{\pi(1+i)e^{-i}}{2}.$$

    相加，并用 $e^{\pm i}=\cos1\pm i\sin1$ 化简：

    $$I=\frac\pi2\Big[(1-i)e^{i}-(1+i)e^{-i}\Big]=\frac\pi2\cdot2i(\sin1-\cos1)=\pi i(\sin1-\cos1)=\sqrt2\,\pi i\,\sin\!\Big(1-\frac\pi4\Big).$$

### 3.4.2 柯西不等式与刘维尔定理

!!! theorem "定理 3.11（柯西不等式）"

    设 $f(z)$ 在 $|z-z_0|<R$ 内解析，且 $|f(z)|\le M$，则

    $$\big|f^{(n)}(z_0)\big|\le\frac{n!\,M}{R^n}\qquad(n=1,2,\dots).$$

（对 $|z-z_0|=r<R$ 用高阶导数公式和长大不等式：$|f^{(n)}(z_0)|\le\frac{n!}{2\pi}\cdot\frac{M}{r^{n+1}}\cdot2\pi r$，再令 $r\to R$。）

!!! theorem "定理 3.12（刘维尔定理）"

    在整个复平面上解析且有界的函数必为常数。

??? proof "证明"

    设 $|f(z)|\le M$。对任意 $z_0$ 和任意 $R>0$，由柯西不等式 $|f'(z_0)|\le\frac MR$。令 $R\to\infty$ 得 $f'(z_0)=0$。所以 $f'\equiv0$，$f$ 为常数（第 2 章练习 8）。$\blacksquare$

**推论（代数基本定理）。** 次数 $n\ge1$ 的复系数多项式 $P(z)$ 至少有一个零点。

否则 $\frac{1}{P(z)}$ 在全平面解析。又 $|P(z)|\to\infty$（$z\to\infty$），所以 $\frac1P$ 有界，由刘维尔定理 $\frac1P$ 为常数，矛盾。

!!! warning "易错点：$\sin z$ 与刘维尔定理"

    $\sin z$ 在全平面解析，但不是常数，所以它**无界**（第 2 章：$|\sin(iy)|=|\sinh y|\to\infty$）。“全平面解析 + 有界 ⇒ 常数”这一结论在实数中不成立（$\sin x$ 是反例）。

---

## 知识点归纳总结

!!! abstract "知识点归纳总结"

    ```mermaid
    flowchart TD
        DEF["复积分 = 黎曼和极限<br/>参数法 ∫ f(z(t)) z'(t) dt"] --> EST["长大不等式 |∫| ≤ ML"]
        DEF --> KEY["∮ dz/(z−z₀)ⁿ = 2πi (n=1), 0 (n≠1)"]
        CR["第2章：解析 / C-R"] --> CT["柯西定理<br/>单连通 + 解析 ⇒ ∮ = 0"]
        CT --> DEFORM["闭路变形 / 复合闭路<br/>只看围住哪些奇点"]
        CT --> PRIM["路径无关 ⇒ 原函数<br/>牛顿-莱布尼茨"]
        DEFORM --> CIF["柯西积分公式<br/>∮ f/(z−z₀) = 2πi f(z₀)"]
        KEY --> CIF
        CIF --> MEAN["平均值公式 → 最大模原理"]
        CIF --> HD["高阶导数公式<br/>∮ f/(z−z₀)ⁿ⁺¹ = 2πi f⁽ⁿ⁾(z₀)/n!"]
        HD --> INF["解析 ⇒ 无穷次可导"]
        HD --> LIOU["柯西不等式 → 刘维尔定理<br/>→ 代数基本定理"]
        INF -.->|"第4章"| TAY["泰勒级数"]
        HD -.->|"第5章"| RES["留数定理"]
    ```

    **核心结论**

    | 内容 | 结论 |
    |---|---|
    | 参数法 | $\int_Cf\,\mathrm dz=\int_\alpha^\beta f(z(t))z'(t)\,\mathrm dt$ |
    | 估值 | $\left\lvert\int_Cf\,\mathrm dz\right\rvert\le\int_C\lvert f\rvert\,\mathrm ds\le ML$ |
    | 关键积分 | $\oint_{\lvert z-z_0\rvert=r}\frac{\mathrm dz}{(z-z_0)^n}=2\pi i\ (n=1)$，其余整数 $n$ 为 $0$ |
    | 柯西定理 | 单连通域内解析 ⇒ $\oint_Cf\,\mathrm dz=0$ |
    | 复合闭路 | $\oint_Cf\,\mathrm dz=\sum_k\oint_{C_k}f\,\mathrm dz$ |
    | 原函数 | $\int_{z_0}^{z_1}f\,\mathrm dz=G(z_1)-G(z_0)$（单连通域内解析） |
    | 柯西积分公式 | $\oint_C\frac{f(z)}{z-z_0}\,\mathrm dz=2\pi if(z_0)$ |
    | 高阶导数 | $\oint_C\frac{f(z)}{(z-z_0)^{n+1}}\,\mathrm dz=\frac{2\pi i}{n!}f^{(n)}(z_0)$ |
    | 柯西不等式 | $\lvert f^{(n)}(z_0)\rvert\le n!M/R^n$ |
    | 刘维尔 | 全平面解析且有界 ⇒ 常数 |

    **闭路积分 $\oint_C\frac{P(z)}{Q(z)}\,\mathrm dz$ 的决策流程**

    1. 求 $Q$ 的零点，判断哪些在 $C$ **内部**；
    2. 内部没有奇点 → $0$（柯西定理）；
    3. 内部只有一个奇点 $z_0$（$m$ 阶）→ 写成 $\frac{f(z)}{(z-z_0)^m}$，用柯西积分公式（$m=1$）或高阶导数公式（$m\ge2$）；
    4. 内部有多个奇点 → 复合闭路定理拆成小圆，或部分分式；
    5. 被积函数不解析（含 $\bar z,\ |z|,\ \operatorname{Re}z$）→ 参数法；在圆周 $|z|=r$ 上可用 $\bar z=\frac{r^2}{z}$ 化成解析函数。

    **易错点汇总**

    - 柯西定理需要解析和单连通（$C$ 内无奇点）；$\oint=0$ 推不出解析；
    - 牛顿-莱布尼茨公式要求在包含路径的单连通区域内解析；
    - 柯西积分公式中 $z_0$ 必须在 $C$ 内部，$f$ 必须在 $C$ 内解析；
    - 高阶导数公式：分母 $(z-z_0)^{n+1}$ 对应 $\frac{2\pi i}{n!}f^{(n)}$，指数与阶数差 $1$；
    - $\mathrm dz$ 与 $\mathrm ds=|\mathrm dz|$ 不要混淆；
    - 闭路积分默认逆时针，顺时针要变号。

---

## 练习题

!!! example "练习 1"

    计算 $\displaystyle\int_C(x-y+ix^2)\,\mathrm dz$，其中 $C$ 为从 $0$ 到 $1+i$ 的直线段。

??? success "解答"

    $C:z=(1+i)t,\ t:0\to1$，即 $x=y=t$，$\mathrm dz=(1+i)\,\mathrm dt$。被积函数为 $0+it^2$：

    $$\int_0^1it^2(1+i)\,\mathrm dt=\frac{i(1+i)}{3}=\frac{-1+i}{3}.$$

    被积函数不解析，不能用原函数。

!!! example "练习 2"

    计算 $\displaystyle\oint_{|z|=2}\frac{\bar z}{|z|}\,\mathrm dz$。

??? success "解答"

    在 $|z|=2$ 上 $|z|=2$，$z\bar z=4\Rightarrow\bar z=\frac4z$，所以

    $$\oint_{|z|=2}\frac{\bar z}{|z|}\,\mathrm dz=\oint_{|z|=2}\frac{2}{z}\,\mathrm dz=2\cdot2\pi i=4\pi i .$$

    （也可以用参数法：$z=2e^{i\theta}$，被积式为 $e^{-i\theta}\cdot2ie^{i\theta}\mathrm d\theta=2i\,\mathrm d\theta$。）

!!! example "练习 3"

    证明 $\left|\displaystyle\oint_{|z|=2}\frac{\mathrm dz}{z^2+1}\right|\le\frac{4\pi}{3}$，并求该积分的准确值。

??? success "解答"

    在 $|z|=2$ 上 $|z^2+1|\ge|z|^2-1=3$，长度 $4\pi$，故积分的模 $\le\frac{4\pi}{3}$。准确值为 $0$（例题 3.10）。估值只给出上界，不一定接近真实值。

!!! example "练习 4"

    计算 $\displaystyle\oint_{|z|=2}\frac{\mathrm dz}{z^2(z-1)}$。

??? success "解答"

    部分分式：$\dfrac{1}{z^2(z-1)}=-\dfrac1z-\dfrac1{z^2}+\dfrac1{z-1}$（通分验证：$-z(z-1)-(z-1)+z^2=1$）。三个奇点项都在圆内，由例题 3.2：

    $$I=-2\pi i-0+2\pi i=0 .$$

    或者用高阶导数公式：在 $z=0$ 附近取 $f_1=\frac1{z-1}$，$f_1'(0)=-1$，贡献 $-2\pi i$；在 $z=1$ 附近取 $f_2=\frac1{z^2}$，$f_2(1)=1$，贡献 $2\pi i$。

!!! example "练习 5"

    计算 $\displaystyle\oint_{|z|=1}\frac{\sin z}{z^4}\,\mathrm dz$ 与 $\displaystyle\oint_{|z|=1}\frac{e^{2z}}{z^3}\,\mathrm dz$。

??? success "解答"

    (1) $n=3$：$\dfrac{2\pi i}{3!}(\sin z)'''\big|_0=\dfrac{2\pi i}{6}\cdot(-\cos0)=-\dfrac{\pi i}{3}$。

    (2) $n=2$：$\dfrac{2\pi i}{2!}(e^{2z})''\big|_0=\pi i\cdot4=4\pi i$。

!!! example "练习 6"

    计算 $\displaystyle\oint_{|z|=4}\frac{z}{(9-z^2)(z+i)}\,\mathrm dz$（例题 3.9(2) 换成 $|z|=4$）。

??? success "解答"

    此时 $-i,\ 3,\ -3$ 都在圆内。$\frac{z}{(9-z^2)(z+i)}=\frac{-z}{(z-3)(z+3)(z+i)}$，分别用柯西积分公式：

    - $z=-i$：$\frac{\pi}{5}$（例题 3.9）；
    - $z=3$：$2\pi i\cdot\dfrac{-3}{6(3+i)}=-\dfrac{\pi i}{3+i}=\dfrac{-\pi-3\pi i}{10}$；
    - $z=-3$：$2\pi i\cdot\dfrac{3}{(-6)(-3+i)}=\dfrac{\pi i}{3-i}=\dfrac{-\pi+3\pi i}{10}$。

    相加：$\frac\pi5-\frac{2\pi}{10}=0$。

    （原因：被积函数在 $\infty$ 附近像 $-\frac1{z^2}$，在大圆上积分 $\to0$，由闭路变形，所有大圆上的积分都相等，只能为 $0$。第 5 章将用“无穷远点的留数”解释这件事。）

!!! example "练习 7"

    计算 $\displaystyle\int_0^{\pi i}ze^{z}\,\mathrm dz$ 与 $\displaystyle\int_1^{i}\frac{\ln(z+1)}{z+1}\,\mathrm dz$（沿 $|z|=1$ 在第一象限的弧，$\ln$ 取主值）。

??? success "解答"

    (1) $\int ze^z\,\mathrm dz=(z-1)e^z$：$\big[(z-1)e^z\big]_0^{\pi i}=(\pi i-1)(-1)-(-1)=2-\pi i$。

    (2) $\ln(z+1)$ 在 $\mathbb C\setminus(-\infty,-1]$ 内解析，该区域是单连通的且包含积分路径。原函数为 $\frac12\ln^2(z+1)$：

    $$\frac12\Big[\ln^2(1+i)-\ln^22\Big],\qquad \ln(1+i)=\tfrac12\ln2+\tfrac{\pi}{4}i .$$

    展开得 $\frac12\left[\frac{\ln^22}{4}-\frac{\pi^2}{16}+\frac{\pi\ln2}{4}i-\ln^22\right]=-\frac{\pi^2}{32}-\frac{3\ln^22}{8}+\frac{\pi\ln2}{8}i$。

!!! example "练习 8"

    设 $f(z)$ 在全平面解析，且 $\operatorname{Re}f(z)\le M$ 对所有 $z$ 成立。证明 $f$ 为常数。

??? success "解答"

    令 $g(z)=e^{f(z)}$，它在全平面解析，且 $|g(z)|=e^{\operatorname{Re}f(z)}\le e^M$ 有界。由刘维尔定理 $g$ 为常数，所以 $g'=f'e^f\equiv0$。又 $e^f\ne0$，故 $f'\equiv0$，$f$ 为常数。

---

**下一章预告。** 高阶导数公式说解析函数无穷次可导，那么它能不能像 $e^x$ 那样展开成幂级数？第 4 章给出答案：能，而且收敛半径**恰好等于到最近奇点的距离**。在奇点附近则要用同时含正、负幂的洛朗级数，它的 $(z-z_0)^{-1}$ 项系数正是第 5 章的“留数”。→ [第 4 章 解析函数的级数表示](ch4.md)
