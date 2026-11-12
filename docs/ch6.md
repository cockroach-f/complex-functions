# 第 6 章 共形映射

!!! tip "本章目标"

    1. 理解导数的几何意义：$\arg f'(z_0)$ 是**旋转角**，$|f'(z_0)|$ 是**伸缩率**；理解**保角性**与**共形映射**的定义；
    2. 掌握求象区域的方法（保域性定理、边界对应原理），了解黎曼映射定理；
    3. 掌握**分式线性映射**的分解与四个性质（保形、保圆、保对称、三点唯一确定），熟记上半平面 $\to$ 单位圆、单位圆 $\to$ 单位圆、上半平面 $\to$ 上半平面的标准形式；
    4. 掌握幂函数 $w=z^n$（角形域）与指数函数 $w=e^z$（带形域）的映射性质，会**组合**这些映射，把月牙形、角形、带形、半圆等区域映射为上半平面或单位圆——期末“求象区域”“构造共形映射”两类题。

**本章主线。** 第 2 章的 C-R 方程说：可导的复函数在局部就是“乘以复数 $f'(z_0)$”，即旋转 + 伸缩。旋转和伸缩都不改变角度，所以解析函数（在 $f'\ne0$ 的点）**保持曲线的夹角**。本章研究这种映射能把区域变成什么样，以及怎样把复杂区域变成单位圆或上半平面。这是求解二维拉普拉斯方程（静电场、热传导、流体）的标准工具，也是 Smith 圆图的数学基础。

---

## 6.1 共形映射的概念

### 6.1.1 导数的几何意义

设 $w=f(z)$ 在区域 $D$ 内解析，$z_0\in D$，$f'(z_0)\ne0$。$C:z=z(t)$ 是过 $z_0=z(t_0)$ 的光滑曲线，$z'(t_0)\ne0$，它的切线方向角为 $\arg z'(t_0)$。像曲线 $\Gamma:w=f(z(t))$ 在 $w_0=f(z_0)$ 处的切向量为

$$w'(t_0)=f'(z_0)\,z'(t_0),\qquad\arg w'(t_0)=\arg f'(z_0)+\arg z'(t_0),\qquad|\mathrm dw|=|f'(z_0)|\,|\mathrm dz| .$$

!!! definition "定义 6.1（旋转角与伸缩率）"

    称 $\arg w'(t_0)-\arg z'(t_0)$ 为映射 $w=f(z)$ 沿曲线 $C$ 在 $z_0$ 处的**旋转角**；称

    $$\lim_{z\to z_0,\,z\in C}\frac{|f(z)-f(z_0)|}{|z-z_0|}$$

    为映射沿 $C$ 在 $z_0$ 处的**伸缩率**。

!!! theorem "定理 6.1（导数的几何意义）"

    设 $f$ 在 $z_0$ 处解析且 $f'(z_0)\ne0$，则映射 $w=f(z)$ 在 $z_0$ 处：

    1. **旋转角不变性**：沿过 $z_0$ 的任何光滑曲线，旋转角都等于 $\arg f'(z_0)$，与曲线无关；
    2. **伸缩率不变性**：沿过 $z_0$ 的任何光滑曲线，伸缩率都等于 $|f'(z_0)|$，与曲线无关；
    3. **保角性**：过 $z_0$ 的任意两条光滑曲线的夹角，在大小和方向上都保持不变。

（第 3 条：两条曲线的切线都转过同一个角 $\arg f'(z_0)$，所以夹角不变。）

<figure class="video" markdown>
<video controls preload="metadata" width="100%" poster="assets/posters/C6_Conformality.jpg">
<source src="assets/videos/C6_Conformality.mp4" type="video/mp4">
</video>
<figcaption>动画 6-1　保角性：在 $z_0=1+0.6i$ 处，$w=z^2$ 把两条切线都旋转 $\arg f'(z_0)$，夹角不变；在 $z_0=0$ 处 $f'(0)=0$，夹角被加倍，保角性失效。</figcaption>
</figure>

!!! example "例题 6.1"

    求 $w=z^2$ 在 $z_0=i$ 处的旋转角与伸缩率；$w=z^2$ 在哪些点保角？

    **解。** $f'(z)=2z$，$f'(i)=2i$，旋转角 $\arg(2i)=\frac\pi2$，伸缩率 $|2i|=2$。$f'(z)=2z\ne0\iff z\ne0$，所以除原点外处处保角。原点处，过原点夹角为 $\alpha$ 的两条射线被映成夹角为 $2\alpha$ 的射线，夹角不保持。

### 6.1.2 共形映射

!!! definition "定义 6.2（保角映射与共形映射）"

    1. 若 $w=f(z)$ 在区域 $D$ 内每一点都具有保角性与伸缩率不变性，则称 $f$ 为 $D$ 内的**第一类保角映射**；
    2. 若 $f$ 是 $D$ 内的第一类保角映射，且是**一一映射**（双方单值），则称 $f$ 为 $D$ 内的**共形映射**。

    （若保持夹角大小而方向相反，称为**第二类保角映射**，例如 $w=\bar z$。）

!!! theorem "定理 6.2"

    若 $f$ 在区域 $D$ 内解析且 $f'(z)\ne0$，则 $f$ 是 $D$ 内的第一类保角映射。若 $f$ 还在 $D$ 内是一一的，则 $f$ 是 $D$ 内的共形映射。

!!! example "例题 6.2"

    讨论 (1) $w=e^z$；(2) $w=z^2$ 的共形性。

    **解。** (1) $(e^z)'=e^z\ne0$，所以 $e^z$ 在整个复平面上是第一类保角映射。但 $e^{z+2\pi i}=e^z$，它在全平面不是一一的，不是共形映射；在任一宽度不超过 $2\pi$ 的水平带形内（如 $0<\operatorname{Im}z<2\pi$）是共形映射。

    (2) $z^2$ 在 $z\ne0$ 处保角，但 $z$ 与 $-z$ 的像相同。在上半平面 $\operatorname{Im}z>0$ 内它是一一的，是共形映射（把上半平面映成 $\mathbb C\setminus[0,+\infty)$）。

!!! warning "易错点：保角 ≠ 共形"

    “共形”还要求一一。$e^z$、$z^n$ 处处（或几乎处处）保角，但只在适当的区域内才是共形映射。题目给出区域时，要检查映射在这个区域内是否一一。

---

## 6.2 共形映射的基本问题

### 6.2.1 问题一：已知区域与映射，求象区域

!!! theorem "定理 6.3（保域性定理）"

    若 $w=f(z)$ 在区域 $D$ 内解析且不恒为常数，则 $D$ 的象 $G=f(D)$ 也是一个区域。

!!! theorem "定理 6.4（边界对应原理）"

    设区域 $D$ 的边界为简单闭曲线 $C$，$f$ 在闭域 $\bar D=D\cup C$ 上解析，且把 $C$ 一一地映射为简单闭曲线 $\Gamma$。则 $\Gamma$ 所围成的某个区域 $G$ 就是 $D$ 的象，且 $f$ 把 $D$ 共形映射为 $G$；当 $z$ 沿 $C$ 正向绕行时，$w$ 沿 $\Gamma$ 也按正向绕行（$G$ 在 $\Gamma$ 的左手边）。

!!! success "求象区域的一般步骤"

    1. 写出边界 $C$ 的方程（参数方程或 $x,y$ 的关系式）；
    2. 求出边界的像 $\Gamma$：把 $z=x+iy$ 代入得 $u(x,y),v(x,y)$，消去参数；
    3. 确定象区域在 $\Gamma$ 的哪一侧：取 $D$ 内一个**测试点**，看它的像落在哪一侧（或用“沿边界正向走，区域在左边”）。

!!! example "例题 6.3"

    求区域 $D=\{z:\operatorname{Re}z>1\}$ 在 $w=\dfrac1z$ 下的像。

    **解。** 边界 $x=1$：$z=1+iy$，

    $$w=\frac{1}{1+iy}=\frac{1-iy}{1+y^2},\qquad u=\frac{1}{1+y^2},\quad v=\frac{-y}{1+y^2} .$$

    $u^2+v^2=\frac{1}{1+y^2}=u$，即 $\left(u-\frac12\right)^2+v^2=\frac14$，是圆 $\left|w-\frac12\right|=\frac12$（$y\to\pm\infty$ 对应 $w=0$）。测试点 $z=2\in D$ 的像 $w=\frac12$ 在圆内，所以象区域为 $\left|w-\frac12\right|<\frac12$。

### 6.2.2 问题二：已知两个区域，求共形映射

!!! theorem "定理 6.5（黎曼映射定理）"

    设 $D$ 是单连通区域，其边界至少包含两个点，则存在解析函数 $w=f(z)$，把 $D$ 共形映射为单位圆盘 $|w|<1$。若再给定 $z_0\in D$，要求 $f(z_0)=0$，$\arg f'(z_0)=\theta_0$，则这样的映射是**唯一**的。

**说明。**

1. 由此，任意两个这样的单连通区域之间都存在共形映射（经单位圆盘中转）。
2. “边界至少两个点”不能去掉：若 $D=\mathbb C$（边界为空集或只有 $\infty$ 一点），共形映射 $f:D\to\{|w|<1\}$ 将是全平面上解析且有界的函数，由刘维尔定理（定理 3.12）$f$ 为常数，矛盾。
3. 黎曼定理只断言存在，不给构造方法。下面两节给出构造这类映射的基本工具：**分式线性映射**与**初等函数**。

---

## 6.3 分式线性映射

### 6.3.1 定义与分解

!!! definition "定义 6.3（分式线性映射）"

    称

    $$w=\frac{az+b}{cz+d}\qquad(a,b,c,d\in\mathbb C,\ ad-bc\ne0)$$

    为**分式线性映射**（默比乌斯变换）。在扩充复平面上补充定义：$c\ne0$ 时，$-\frac dc\mapsto\infty$，$\infty\mapsto\frac ac$；$c=0$ 时 $\infty\mapsto\infty$。

$ad-bc\ne0$ 保证它不是常数：$w'=\frac{ad-bc}{(cz+d)^2}\ne0$。它的逆映射 $z=\frac{-dw+b}{cw-a}$ 也是分式线性映射，所以分式线性映射是扩充复平面到自身的一一映射。

!!! theorem "定理 6.6（分解）"

    任一分式线性映射都可以分解为下列简单映射的复合：

    1. **平移** $w=z+b$；
    2. **旋转与伸缩** $w=az$（$a=re^{i\alpha}$：旋转 $\alpha$，伸缩 $r$ 倍）；
    3. **反演（倒数映射）** $w=\dfrac1z$。

    事实上，$c=0$ 时 $w=\frac adz+\frac bd$；$c\ne0$ 时

    $$w=\frac ac+\frac{bc-ad}{c}\cdot\frac{1}{cz+d}.$$

**反演的几何意义。** $w=\frac1z=\frac{1}{r}e^{-i\theta}$ 可分两步完成：先作关于单位圆的对称 $w_1=\frac{1}{\bar z}$（模变为倒数，辐角不变），再作关于实轴的对称 $w=\bar w_1$。

!!! definition "定义 6.4（关于圆周的对称点）"

    设圆周 $C:|z-z_0|=R$。若点 $z_1,z_2$ 在从 $z_0$ 出发的同一条射线上，且

    $$|z_1-z_0|\cdot|z_2-z_0|=R^2,$$

    则称 $z_1,z_2$ 关于圆周 $C$ **对称**（约定 $z_0$ 与 $\infty$ 关于 $C$ 对称）。等价地，$z_2=z_0+\dfrac{R^2}{\overline{z_1-z_0}}$；又等价于：过 $z_1,z_2$ 的任何圆周都与 $C$ 正交。

### 6.3.2 分式线性映射的性质

!!! theorem "定理 6.7（保形性）"

    分式线性映射在扩充复平面上是共形映射。

（在 $z\ne-\frac dc,\infty$ 处 $w'\ne0$；在 $-\frac dc$ 与 $\infty$ 处，通过 $\zeta=\frac1z$ 定义角度后也保角。）

!!! theorem "定理 6.8（保圆性）"

    在扩充复平面上，分式线性映射把圆周映射为圆周（直线看作过 $\infty$ 的圆周）。

??? proof "证明"

    平移、旋转、伸缩显然把圆周映为圆周，只需证明 $w=\frac1z$ 的情形。由第 1 章，圆周或直线的方程为

    $$Az\bar z+\bar\beta z+\beta\bar z+C=0\qquad(A,C\in\mathbb R,\ |\beta|^2>AC).$$

    令 $z=\frac1w$，乘以 $w\bar w$ 得 $Cw\bar w+\beta w+\bar\beta\bar w+A=0$，仍是同一类方程。$\blacksquare$

!!! success "像是圆还是直线？"

    设 $w=\frac{az+b}{cz+d}$，$c\ne0$。$z$ 平面上的圆（或直线）$C$：

    - 若 $C$ **经过** $-\frac dc$（映到 $\infty$ 的点），则像是**直线**；
    - 否则像是**有限的圆**。

<figure class="video" markdown>
<video controls preload="metadata" width="100%" poster="assets/posters/C6_Inversion.jpg">
<source src="assets/videos/C6_Inversion.mp4" type="video/mp4">
</video>
<figcaption>动画 6-2　反演 $w=\frac1z$：竖直线变成过原点的圆，过原点的圆变成直线，不过原点的圆变成圆；单位圆内部的点被映到外部。</figcaption>
</figure>

!!! theorem "定理 6.9（保对称性）"

    若 $z_1,z_2$ 关于圆周 $C$ 对称，则在分式线性映射下，它们的像 $w_1,w_2$ 关于像圆周 $\Gamma$ 对称。

（过 $z_1,z_2$ 的圆都与 $C$ 正交；由保圆性和保角性，它们的像是过 $w_1,w_2$ 且与 $\Gamma$ 正交的圆，所以 $w_1,w_2$ 关于 $\Gamma$ 对称。）

!!! theorem "定理 6.10（三点唯一确定）"

    给定扩充 $z$ 平面上三个互异点 $z_1,z_2,z_3$ 与扩充 $w$ 平面上三个互异点 $w_1,w_2,w_3$，存在唯一的分式线性映射把 $z_k$ 映为 $w_k$（$k=1,2,3$），它由

    $$\frac{w-w_1}{w-w_2}\cdot\frac{w_3-w_2}{w_3-w_1}=\frac{z-z_1}{z-z_2}\cdot\frac{z_3-z_2}{z_3-z_1}$$

    确定（**交比不变性**）。若某点为 $\infty$，把含它的因子换成 $1$（例如 $w_3=\infty$ 时，$\frac{w_3-w_2}{w_3-w_1}$ 换成 $1$）。

!!! example "例题 6.4"

    求把 $z_1=2,\ z_2=i,\ z_3=-2$ 分别映为 $w_1=-1,\ w_2=i,\ w_3=1$ 的分式线性映射。

    **解。** 由交比公式，

    $$\frac{w+1}{w-i}\cdot\frac{1-i}{2}=\frac{z-2}{z-i}\cdot\frac{-2-i}{-4}.$$

    整理（把 $\frac{w+1}{w-i}$ 解出来再解 $w$）得

    $$w=\frac{z-6i}{3iz-2}.$$

    **验算：** $z=2$：$\frac{2-6i}{6i-2}=-1$ ✓；$z=i$：$\frac{-5i}{-5}=i$ ✓；$z=-2$：$\frac{-2-6i}{-6i-2}=1$ ✓。

    （几何上：$2,i,-2$ 在圆 $|z|=2$ 上，像 $-1,i,1$ 在单位圆上；按保圆性，该映射把圆 $|z|=2$ 映为单位圆。）

### 6.3.3 三个标准映射

!!! theorem "定理 6.11（上半平面 → 单位圆盘）"

    把上半平面 $\operatorname{Im}z>0$ 共形映射为单位圆盘 $|w|<1$、且把上半平面内的点 $\lambda$ 映为 $w=0$ 的分式线性映射为

    $$w=e^{i\theta}\,\frac{z-\lambda}{z-\bar\lambda}\qquad(\operatorname{Im}\lambda>0,\ \theta\in\mathbb R).$$

??? proof "推导"

    $\lambda$ 与 $\bar\lambda$ 关于实轴对称。实轴要映为单位圆，$\lambda\mapsto0$，由保对称性 $\bar\lambda$ 必须映为 $0$ 关于单位圆的对称点 $\infty$。所以 $w=k\frac{z-\lambda}{z-\bar\lambda}$。实轴上 $|x-\lambda|=|x-\bar\lambda|$，要使 $|w|=1$，需 $|k|=1$，即 $k=e^{i\theta}$。$\blacksquare$

!!! theorem "定理 6.12（单位圆盘 → 单位圆盘）"

    把 $|z|<1$ 共形映射为 $|w|<1$、且把 $\alpha$（$|\alpha|<1$）映为 $0$ 的分式线性映射为

    $$w=e^{i\theta}\,\frac{z-\alpha}{1-\bar\alpha z}\qquad(\theta\in\mathbb R).$$

    （$\alpha$ 关于单位圆的对称点为 $\frac{1}{\bar\alpha}$，推导同上。）此时 $f'(\alpha)=\dfrac{e^{i\theta}}{1-|\alpha|^2}$，所以 $\arg f'(\alpha)=\theta$。

!!! theorem "定理 6.13（上半平面 → 上半平面）"

    $w=\dfrac{az+b}{cz+d}$，其中 $a,b,c,d$ 为**实数**且 $ad-bc>0$。

<figure class="video" markdown>
<video controls preload="metadata" width="100%" poster="assets/posters/C6_HalfPlaneToDisk.jpg">
<source src="assets/videos/C6_HalfPlaneToDisk.mp4" type="video/mp4">
</video>
<figcaption>动画 6-3　$w=\frac{z-i}{z+i}$：上半平面的水平线、竖直线变成单位圆内的两族正交圆弧，实轴变成单位圆周，$i\mapsto0$，$0\mapsto-1$，$1\mapsto-i$。</figcaption>
</figure>

!!! example "例题 6.5"

    (1) 求把上半平面映为单位圆盘、且满足 $f(i)=0$、$\arg f'(i)=-\frac\pi2$ 的映射；(2) 求把单位圆盘映为单位圆盘、且满足 $f\!\left(\frac12\right)=0$、$\arg f'\!\left(\frac12\right)=0$ 的映射。

    **解。** (1) 由定理 6.11，$w=e^{i\theta}\frac{z-i}{z+i}$。$f'(z)=e^{i\theta}\frac{2i}{(z+i)^2}$，$f'(i)=e^{i\theta}\frac{2i}{(2i)^2}=\frac{e^{i\theta}}{2i}$，$\arg f'(i)=\theta-\frac\pi2=-\frac\pi2$，得 $\theta=0$：

    $$w=\frac{z-i}{z+i}.$$

    (2) 由定理 6.12，$\arg f'(\alpha)=\theta=0$：

    $$w=\frac{z-\frac12}{1-\frac12z}=\frac{2z-1}{2-z}.$$

!!! example "例题 6.6（教材例题）"

    求圆弧 $C_1:|z|=1,\ \operatorname{Im}z\ge0$ 与线段 $C_2:[-1,1]$ 在 $w=\dfrac{z+1}{1-z}$ 下的像，并求上半单位圆盘 $D=\{|z|<1,\operatorname{Im}z>0\}$ 的像。

    **解。** $z=-1\mapsto0$，$z=1\mapsto\infty$，所以 $C_1,C_2$ 都映成从原点出发的射线（过 $\infty$ 的“圆”）。

    - $C_2$ 上取 $z=0\mapsto1$：$C_2$ 映为正实轴；
    - $C_1$ 上取 $z=i\mapsto\frac{1+i}{1-i}=i$：$C_1$ 映为正虚轴；
    - 测试点 $z=\frac i2\mapsto\frac{1+i/2}{1-i/2}=\frac{(2+i)^2}{5}=\frac{3+4i}{5}$，在第一象限。

    所以 $D$ 的像为第一象限 $0<\arg w<\frac\pi2$。（也可用保角性：$C_1,C_2$ 在 $-1$ 处夹角为 $\frac\pi2$，像射线夹角也是 $\frac\pi2$。）

    **启示。** 两段圆弧围成的区域，把一个交点映为 $0$、另一个交点映为 $\infty$，就变成**角形域**。这是 6.4 节组合映射的第一步。

!!! warning "易错点"

    - 定理 6.11 中分母是 $z-\bar\lambda$，不是 $z+\lambda$（只有 $\lambda=i$ 时才有 $z-\bar\lambda=z+i$）；
    - 定理 6.12 中分母是 $1-\bar\alpha z$，共轭不能漏；
    - 用交比公式时，$z_k$ 与 $w_k$ 的对应顺序必须一致；
    - 判断象区域时一定要用测试点或方向确认是圆的“内部”还是“外部”。

---

## 6.4 几个初等函数构成的共形映射

### 6.4.1 幂函数 $w=z^n$

$z=re^{i\theta}\mapsto w=r^ne^{in\theta}$：模取 $n$ 次方，辐角乘以 $n$。$(z^n)'=nz^{n-1}\ne0$（$z\ne0$），所以：

!!! success "幂函数的映射性质"

    $w=z^n$（$n\ge2$ 为整数）把角形域 $0<\arg z<\alpha$（$0<\alpha\le\frac{2\pi}{n}$）共形映射为角形域 $0<\arg w<n\alpha$；张角**放大 $n$ 倍**。

    特别地，$w=z^n$ 把 $0<\arg z<\frac{\pi}{n}$ 映为上半平面。反过来，$w=\sqrt[n]{z}$（取适当分支）把角形域的张角**缩小**为原来的 $\frac1n$。更一般地，$w=z^{\pi/\alpha}$ 把 $0<\arg z<\alpha$ 映为上半平面。

!!! example "例题 6.7"

    求把角形域 $0<\arg z<\frac\pi4$ 映为单位圆盘 $|w|<1$ 的一个共形映射。

    **解。** $\zeta=z^4$ 把它映为上半平面 $0<\arg\zeta<\pi$；再用定理 6.11（取 $\lambda=i$）：

    $$w=\frac{\zeta-i}{\zeta+i}=\frac{z^4-i}{z^4+i}.$$

### 6.4.2 指数函数 $w=e^z$

$z=x+iy\mapsto w=e^xe^{iy}$：$x$ 控制模，$y$ 控制辐角（第 2 章动画 2-3）。

!!! success "指数函数的映射性质"

    1. 竖直线 $x=c$ $\mapsto$ 圆周 $|w|=e^c$；水平线 $y=c$ $\mapsto$ 射线 $\arg w=c$；
    2. 带形域 $0<\operatorname{Im}z<a$（$0<a\le2\pi$）$\mapsto$ 角形域 $0<\arg w<a$；特别地，$0<\operatorname{Im}z<\pi$ $\mapsto$ 上半平面；
    3. 半带形域 $0<\operatorname{Im}z<\pi,\ \operatorname{Re}z<0$ $\mapsto$ 上半单位圆盘 $|w|<1,\ \operatorname{Im}w>0$；
    4. 反过来，对数函数 $w=\ln z$ 把角形域映为带形域。

!!! example "例题 6.8"

    求把带形域 $a<\operatorname{Re}z<b$ 映为上半平面的共形映射。

    **解。** 先把竖直带形转成标准的水平带形 $0<\operatorname{Im}\zeta<\pi$：平移 $z-a$，旋转 $90^\circ$（乘以 $i$），伸缩 $\frac{\pi}{b-a}$，即 $\zeta=\frac{\pi i}{b-a}(z-a)$。（验证：$\operatorname{Re}z=a\mapsto\operatorname{Im}\zeta=0$，$\operatorname{Re}z=b\mapsto\operatorname{Im}\zeta=\pi$。）再取指数：

    $$w=e^{\zeta}=\exp\!\left(\frac{\pi i(z-a)}{b-a}\right).$$

### 6.4.3 组合映射：构造共形映射的套路

!!! success "把区域映为上半平面 / 单位圆的四步法"

    1. **（预处理）** 平移、旋转、伸缩，使区域处于便于处理的位置；
    2. **化为角形域或带形域**：区域由两段相交圆弧（或直线）围成时，用分式线性映射把一个交点映为 $0$、另一个映为 $\infty$，得到**角形域**；两段相切的圆弧，把切点映为 $\infty$，得到**带形域**；
    3. **化为上半平面**：角形域用幂函数 $\zeta^{\pi/\alpha}$；带形域用指数函数；
    4. **化为单位圆盘**：用定理 6.11 $w=e^{i\theta}\frac{\eta-\lambda}{\eta-\bar\lambda}$。

!!! example "例题 6.9（上半单位圆盘 → 单位圆盘）"

    求把 $D=\{|z|<1,\ \operatorname{Im}z>0\}$ 映为单位圆盘 $|w|<1$ 的共形映射。

    **解。**

    1. 由例题 6.6，$\zeta=\dfrac{1+z}{1-z}$ 把 $D$ 映为第一象限 $0<\arg\zeta<\frac\pi2$；
    2. $\eta=\zeta^2$ 把第一象限映为上半平面；
    3. $w=\dfrac{\eta-i}{\eta+i}$ 把上半平面映为单位圆盘。

    复合得

    $$w=\frac{\left(\frac{1+z}{1-z}\right)^2-i}{\left(\frac{1+z}{1-z}\right)^2+i}=\frac{(1+z)^2-i(1-z)^2}{(1+z)^2+i(1-z)^2}.$$

    特别地，**上半单位圆盘 → 上半平面**的映射为 $w=\left(\dfrac{1+z}{1-z}\right)^2$，值得记住。

<figure class="video" markdown>
<video controls preload="metadata" width="100%" poster="assets/posters/C6_Pipeline.jpg">
<source src="assets/videos/C6_Pipeline.mp4" type="video/mp4">
</video>
<figcaption>动画 6-4　组合映射流水线：上半单位圆盘 $\xrightarrow{\frac{1+z}{1-z}}$ 第一象限 $\xrightarrow{\zeta^2}$ 上半平面 $\xrightarrow{\frac{\eta-i}{\eta+i}}$ 单位圆盘。红色直径、蓝色半圆以及内部的同心半圆和射线一路被追踪，始终保持正交。</figcaption>
</figure>

!!! example "例题 6.10（半圆外部）"

    求把 $D=\{\operatorname{Im}z>0,\ |z|>1\}$（上半平面去掉上半单位圆盘）映为上半平面的共形映射。

    **解。** 边界由实轴的两段 $(-\infty,-1]$、$[1,+\infty)$ 与上半单位圆周组成，它们交于 $\pm1$，夹角为 $\frac\pi2$。令 $\zeta=\frac{z-1}{z+1}$（$1\mapsto0$，$-1\mapsto\infty$）：

    - $[1,+\infty)\mapsto[0,1)$（正实轴的一部分），$(-\infty,-1)\mapsto(1,+\infty)$，合起来是正实轴；
    - 上半单位圆周：$z=e^{i\theta}\mapsto\frac{e^{i\theta}-1}{e^{i\theta}+1}=i\tan\frac\theta2$，是正虚轴；
    - 测试点 $z=2i\mapsto\frac{2i-1}{2i+1}=\frac{3+4i}{5}$，在第一象限。

    所以 $\zeta$ 把 $D$ 映为第一象限，再令 $w=\zeta^2$：

    $$w=\left(\frac{z-1}{z+1}\right)^2 .$$

    （另一个常用答案是茹科夫斯基函数 $w=z+\frac1z$：它把 $D$ 映为上半平面，在机翼理论中很重要。）

!!! example "例题 6.11（半带形域）"

    求把半带形域 $D=\{0<\operatorname{Im}z<\pi,\ \operatorname{Re}z<0\}$ 映为上半平面的共形映射。

    **解。** $\zeta=e^z$ 把 $D$ 映为上半单位圆盘（$\operatorname{Re}z<0\iff|\zeta|<1$，$0<\operatorname{Im}z<\pi\iff0<\arg\zeta<\pi$）。再用例题 6.9 的结论：

    $$w=\left(\frac{1+e^z}{1-e^z}\right)^2 .$$

!!! warning "易错点：上半单位圆盘不能直接平方"

    $w=z^2$ 把上半单位圆盘映为**去掉 $[0,1)$ 的单位圆盘**，不是上半平面，也不是单位圆盘。原因：上半单位圆盘不是角形域，它的边界在 $\pm1$ 处有两个“角”，在 $0$ 处没有角。正确做法是先把角点 $\pm1$ 映到 $0$ 和 $\infty$。

---

## 知识点归纳总结

!!! abstract "知识点归纳总结"

    ```mermaid
    flowchart TD
        CR["第2章：解析 / C-R<br/>局部 = 乘以 f'(z₀)"] --> GEO["导数的几何意义<br/>旋转角 arg f'，伸缩率 |f'|"]
        GEO --> CONF["保角性（f' ≠ 0）<br/>+ 一一 ⇒ 共形映射"]
        CONF --> P1["问题一：求象区域<br/>保域性 + 边界对应 + 测试点"]
        CONF --> P2["问题二：求映射<br/>黎曼映射定理（存在唯一）"]
        P2 --> MOB["分式线性映射 (az+b)/(cz+d)<br/>= 平移 ∘ 旋转伸缩 ∘ 反演"]
        MOB --> PROP["保形 · 保圆 · 保对称 · 三点确定"]
        PROP --> STD["标准映射<br/>上半平面→圆：e^{iθ}(z−λ)/(z−λ̄)<br/>圆→圆：e^{iθ}(z−α)/(1−ᾱz)"]
        P2 --> ELEM["初等函数<br/>zⁿ：角形张角 ×n<br/>e^z：带形 → 角形"]
        STD --> COMB["组合：交点→0,∞ → 角形<br/>→ 幂函数 → 上半平面 → 单位圆"]
        ELEM --> COMB
        COMB -.应用.-> APP["二维静电场 / Smith 圆图"]
    ```

    **核心结论**

    | 内容 | 结论 |
    |---|---|
    | 几何意义 | 旋转角 $\arg f'(z_0)$，伸缩率 $\lvert f'(z_0)\rvert$ |
    | 共形映射 | 解析 + $f'\ne0$ + 一一 |
    | 分式线性映射 | $w=\frac{az+b}{cz+d}$，$ad-bc\ne0$；保形、保圆、保对称 |
    | 三点确定 | $\frac{w-w_1}{w-w_2}\cdot\frac{w_3-w_2}{w_3-w_1}=\frac{z-z_1}{z-z_2}\cdot\frac{z_3-z_2}{z_3-z_1}$ |
    | 上半平面 → 单位圆 | $w=e^{i\theta}\frac{z-\lambda}{z-\bar\lambda}$，$\operatorname{Im}\lambda>0$ |
    | 单位圆 → 单位圆 | $w=e^{i\theta}\frac{z-\alpha}{1-\bar\alpha z}$，$\lvert\alpha\rvert<1$ |
    | 上半平面 → 上半平面 | 实系数，$ad-bc>0$ |
    | 幂函数 | $0<\arg z<\alpha\ \xrightarrow{z^n}\ 0<\arg w<n\alpha$ |
    | 指数函数 | $0<\operatorname{Im}z<a\ \xrightarrow{e^z}\ 0<\arg w<a$ |
    | 上半单位圆盘 → 上半平面 | $w=\left(\frac{1+z}{1-z}\right)^2$ |

    **方法清单**

    1. 求象区域：边界参数化 → 求像曲线 → 测试点定侧；分式线性映射可用“三点定圆”；
    2. 构造映射：找出边界的交点（角点），用分式线性映射送到 $0$ 和 $\infty$ → 角形域 → 幂函数 → 上半平面 → 标准映射到单位圆；
    3. 带附加条件（$f(z_0)=w_0$、$\arg f'(z_0)=\theta_0$）：先写出含参数的通式，再代入条件定参数。

    **易错点汇总**

    - 保角不等于共形（还要一一）；$f'(z_0)=0$ 处保角性失效；
    - 标准映射的分母 $z-\bar\lambda$、$1-\bar\alpha z$ 漏写共轭；
    - 判断像在圆内还是圆外，不用测试点就下结论；
    - 直接对非角形域用幂函数；
    - 复合映射的顺序写反。

---

## 练习题

!!! example "练习 1"

    求 $w=z^2$ 在 $z_0=1+i$ 处的旋转角与伸缩率；求 $w=e^z$ 在 $z_0=\ln2+\frac{\pi}{2}i$ 处的旋转角与伸缩率。

??? success "解答"

    (1) $f'(1+i)=2+2i=2\sqrt2e^{i\pi/4}$：旋转角 $\frac\pi4$，伸缩率 $2\sqrt2$。

    (2) $f'(z_0)=e^{z_0}=2e^{i\pi/2}=2i$：旋转角 $\frac\pi2$，伸缩率 $2$。

!!! example "练习 2"

    求单位圆盘 $|z|<1$ 在 $w=\dfrac{z-1}{z+1}$ 下的像。

??? success "解答"

    $z=-1\mapsto\infty$ 在单位圆上，所以单位圆的像是直线。$z=1\mapsto0$，$z=i\mapsto\frac{i-1}{i+1}=i$，所以像直线是虚轴。测试点 $z=0\mapsto-1$，所以像为左半平面 $\operatorname{Re}w<0$。

!!! example "练习 3"

    求带形域 $0<\operatorname{Re}z<1$ 在 $w=\dfrac1z$ 下的像。

??? success "解答"

    $\operatorname{Re}z=0$（虚轴）$\mapsto$ 虚轴；$\operatorname{Re}z=1\mapsto$ 圆 $\left|w-\frac12\right|=\frac12$（例题 6.3）。测试点 $z=\frac12\mapsto2$：在右半平面且在该圆外。所以像为 $\left\{\operatorname{Re}w>0,\ \left|w-\frac12\right|>\frac12\right\}$。

!!! example "练习 4"

    求把 $-1,0,1$ 分别映为 $0,i,\infty$ 的分式线性映射。

??? success "解答"

    $w_3=\infty$ 时交比公式中 $\frac{w_3-w_2}{w_3-w_1}$ 换成 $1$：

    $$\frac{w-0}{w-i}=\frac{z+1}{z-0}\cdot\frac{1-0}{1+1}=\frac{z+1}{2z}.$$

    解得 $2zw=(z+1)(w-i)$，$w=\dfrac{i(z+1)}{1-z}$。验算：$-1\mapsto0$，$0\mapsto i$，$1\mapsto\infty$ ✓。

!!! example "练习 5"

    求把上半平面映为单位圆盘的映射 $f$，使 $f(2i)=0$，$f'(2i)>0$。

??? success "解答"

    $w=e^{i\theta}\frac{z-2i}{z+2i}$，$f'(z)=e^{i\theta}\frac{4i}{(z+2i)^2}$，$f'(2i)=e^{i\theta}\frac{4i}{(4i)^2}=\frac{e^{i\theta}}{4i}=-\frac{i}{4}e^{i\theta}$。要使它为正实数，取 $e^{i\theta}=i$：

    $$w=i\,\frac{z-2i}{z+2i}.$$

!!! example "练习 6"

    求把单位圆盘映为单位圆盘的映射 $f$，使 $f\!\left(\frac i2\right)=0$，$\arg f'\!\left(\frac i2\right)=\frac\pi2$。

??? success "解答"

    $\alpha=\frac i2$，$\bar\alpha=-\frac i2$，$\theta=\arg f'(\alpha)=\frac\pi2$：

    $$w=i\,\frac{z-\frac i2}{1+\frac i2z}=i\,\frac{2z-i}{2+iz}.$$

!!! example "练习 7"

    求把扇形 $\left\{|z|<1,\ 0<\arg z<\frac\pi4\right\}$ 映为上半平面的共形映射。

??? success "解答"

    $\zeta=z^4$ 把扇形映为上半单位圆盘（张角 $\frac\pi4\times4=\pi$，模仍小于 $1$）。再用例题 6.9：

    $$w=\left(\frac{1+z^4}{1-z^4}\right)^2.$$

!!! example "练习 8"

    求把带形域 $0<\operatorname{Im}z<\pi$ 映为单位圆盘、且把 $z=\frac{\pi}{2}i$ 映为 $w=0$ 的共形映射。

??? success "解答"

    $\zeta=e^z$ 把带形映为上半平面，$\frac{\pi}{2}i\mapsto e^{i\pi/2}=i$。再用定理 6.11（$\lambda=i$）：

    $$w=e^{i\theta}\,\frac{e^z-i}{e^z+i}\qquad(\theta\ \text{任意，取}\ \theta=0\ \text{即得一个答案}).$$

---

**全课程到此结束。** 想看六章是怎样一环扣一环的，请读[全课程总结](summary.md)。
