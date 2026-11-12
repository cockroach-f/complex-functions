# 第 4 章 解析函数的级数表示

!!! tip "本章目标"

    1. 理解复数序列、复数项级数收敛与绝对收敛的概念，会把问题化为实部、虚部（或模）的实级数；
    2. 掌握**阿贝尔定理**和幂级数的**收敛圆**，会用比值法、根值法求**收敛半径**；
    3. 掌握**泰勒展开定理**：解析函数在一点附近能展开成幂级数，收敛半径等于**到最近奇点的距离**；牢记 $e^z,\ \sin z,\ \cos z,\ \frac1{1-z},\ \ln(1+z),\ (1+z)^\alpha$ 的展开式；
    4. 掌握**洛朗展开**：先按奇点把平面划分成圆环，再用部分分式与几何级数在每个圆环内展开——期末必考；
    5. 理解洛朗系数 $c_{-1}$ 与闭路积分 $\oint_C f\,\mathrm dz=2\pi i\,c_{-1}$ 的关系，为第 5 章的留数做准备。

**本章主线。** 第 3 章的高阶导数公式说解析函数无穷次可导，于是可以写出泰勒级数。复变函数的一个优点在这里显现出来：泰勒级数的收敛范围是一个**圆盘**，半径由**奇点**的位置决定，这也解释了实分析里一些“莫名其妙”的发散现象。在奇点附近，函数无法展开成泰勒级数，只好让负幂项也参与进来，于是得到洛朗级数。

---

## 4.1 复数项级数

### 4.1.1 复数序列

!!! definition "定义 4.1（复数序列的极限）"

    设 $\{z_n\}$ 是复数序列，$z_0$ 为复数。若对任意 $\varepsilon>0$，存在正整数 $N$，当 $n>N$ 时 $|z_n-z_0|<\varepsilon$，则称 $\{z_n\}$ **收敛**于 $z_0$，记作 $\lim\limits_{n\to\infty}z_n=z_0$。

!!! theorem "定理 4.1"

    设 $z_n=x_n+iy_n$，$z_0=x_0+iy_0$。则 $\lim\limits_{n\to\infty}z_n=z_0\iff\lim\limits_{n\to\infty}x_n=x_0$ 且 $\lim\limits_{n\to\infty}y_n=y_0$。

（由 $|x_n-x_0|,|y_n-y_0|\le|z_n-z_0|\le|x_n-x_0|+|y_n-y_0|$ 即得，和第 1 章复变函数极限的定理 1.1 一样。）

!!! example "例题 4.1"

    讨论下列序列的收敛性：(1) $z_n=\left(1+\dfrac1n\right)e^{i\pi/n}$；(2) $z_n=\dfrac{10000\,i^n}{n!}$；(3) $z_n=i^n$。

    **解。** (1) $x_n=\left(1+\frac1n\right)\cos\frac\pi n\to1$，$y_n=\left(1+\frac1n\right)\sin\frac\pi n\to0$，所以 $z_n\to1$。

    (2) $|z_n|=\dfrac{10000}{n!}\to0$，所以 $z_n\to0$。（$|z_n|\to0\iff z_n\to0$。）

    (3) $z_n$ 在 $i,-1,-i,1$ 之间循环，不收敛；但 $|z_n|=1$ 收敛。**模收敛推不出序列收敛**（极限不是 $0$ 时）。

### 4.1.2 复数项级数

!!! definition "定义 4.2（级数收敛）"

    对复数项级数 $\sum_{n=1}^\infty z_n$，记部分和 $S_n=z_1+z_2+\cdots+z_n$。若 $\{S_n\}$ 收敛于 $S$，则称级数**收敛**，$S$ 称为它的**和**；否则称级数**发散**。

!!! theorem "定理 4.2"

    设 $z_n=x_n+iy_n$，则 $\sum z_n$ 收敛 $\iff$ $\sum x_n$ 与 $\sum y_n$ 都收敛。

    **必要条件：** 若 $\sum z_n$ 收敛，则 $\lim\limits_{n\to\infty}z_n=0$。

!!! definition "定义 4.3（绝对收敛与条件收敛）"

    若 $\sum|z_n|$ 收敛，称 $\sum z_n$ **绝对收敛**；若 $\sum z_n$ 收敛而 $\sum|z_n|$ 发散，称 $\sum z_n$ **条件收敛**。

!!! theorem "定理 4.3"

    1. 若 $\sum|z_n|$ 收敛，则 $\sum z_n$ 收敛（绝对收敛 ⇒ 收敛）；
    2. $\sum z_n$ 绝对收敛 $\iff$ $\sum x_n$ 与 $\sum y_n$ 都绝对收敛。

$\sum|z_n|$ 是**正项**实级数，可以直接用高等数学中的比较法、比值法、根值法。这是判断复级数收敛最常用的手段。

!!! example "例题 4.2（教材例 4.2 改编）"

    判断收敛性：(1) $\displaystyle\sum_{n=1}^\infty\left(\frac1n+\frac{i}{2^n}\right)$；(2) $\displaystyle\sum_{n=1}^\infty\frac{(3+4i)^n}{6^n}$；(3) $\displaystyle\sum_{n=1}^\infty\frac{i^n}{n}$。

    **解。** (1) 实部 $\sum\frac1n$ 发散，故原级数发散（虽然虚部收敛）。

    (2) $\left|\dfrac{(3+4i)^n}{6^n}\right|=\left(\dfrac56\right)^n$，几何级数收敛，所以原级数**绝对收敛**。

    (3) $\left|\frac{i^n}{n}\right|=\frac1n$，$\sum\frac1n$ 发散，所以不绝对收敛。但实部 $-\frac12+\frac14-\frac16+\cdots=\sum_{k\ge1}\frac{(-1)^k}{2k}$ 与虚部 $1-\frac13+\frac15-\cdots$ 都是交错级数，都收敛（莱布尼茨判别法），所以原级数**条件收敛**。

---

## 4.2 幂级数

### 4.2.1 函数项级数与幂级数

!!! definition "定义 4.4（函数项级数、幂级数）"

    设 $\{f_n(z)\}$ 在区域 $D$ 内有定义，$\sum_{n=1}^\infty f_n(z)$ 称为**复变函数项级数**。若对 $D$ 内每一点 $z$ 级数都收敛，其和 $S(z)$ 称为**和函数**。形如

    $$\sum_{n=0}^\infty c_n(z-z_0)^n=c_0+c_1(z-z_0)+c_2(z-z_0)^2+\cdots$$

    的级数称为**幂级数**，$c_n$ 为复常数。令 $\zeta=z-z_0$，总可以化为 $\sum c_n\zeta^n$ 的形式来研究。

最基本的例子是几何级数：

$$\sum_{n=0}^\infty z^n=1+z+z^2+\cdots,\qquad S_N=\frac{1-z^{N+1}}{1-z}\ (z\ne1).$$

当 $|z|<1$ 时 $z^{N+1}\to0$，$S_N\to\frac1{1-z}$；当 $|z|\ge1$ 时 $|z^n|\ge1$，通项不趋于 $0$，发散。

<figure class="video" markdown>
<video controls preload="metadata" width="100%" poster="assets/posters/C4_GeometricPartialSums.jpg">
<source src="assets/videos/C4_GeometricPartialSums.mp4" type="video/mp4">
</video>
<figcaption>动画 4-1　几何级数的部分和：每一步加上 $z^n$（把上一步的向量转 $\arg z$、缩放 $|z|$）。$|z|<1$ 时盘旋收敛到 $\frac1{1-z}$，$|z|>1$ 时向外发散。</figcaption>
</figure>

### 4.2.2 阿贝尔定理与收敛圆

!!! theorem "定理 4.4（阿贝尔定理）"

    对幂级数 $\sum_{n=0}^\infty c_nz^n$：

    1. 若它在 $z_1\ne0$ 处收敛，则在 $|z|<|z_1|$ 内**绝对收敛**；
    2. 若它在 $z_2$ 处发散，则在 $|z|>|z_2|$ 内发散。

??? proof "证明"

    (1) $\sum c_nz_1^n$ 收敛 ⇒ $c_nz_1^n\to0$ ⇒ 存在 $M$ 使 $|c_nz_1^n|\le M$。当 $|z|<|z_1|$ 时，令 $q=\frac{|z|}{|z_1|}<1$，

    $$|c_nz^n|=|c_nz_1^n|\left|\frac{z}{z_1}\right|^n\le Mq^n,$$

    由比较法 $\sum|c_nz^n|$ 收敛。

    (2) 反证：若存在 $|z_3|>|z_2|$ 使级数在 $z_3$ 收敛，由 (1) 级数在 $z_2$ 处绝对收敛，矛盾。$\blacksquare$

由阿贝尔定理，收敛点与发散点之间存在一个分界圆：

!!! definition "定义 4.5（收敛圆与收敛半径）"

    存在 $R\in[0,+\infty]$，使幂级数 $\sum c_n(z-z_0)^n$ 在 $|z-z_0|<R$ 内绝对收敛，在 $|z-z_0|>R$ 内发散。$R$ 称为**收敛半径**，$|z-z_0|<R$ 称为**收敛圆**。$R=0$ 表示只在 $z_0$ 收敛；$R=+\infty$ 表示在全平面收敛。

!!! warning "易错点：收敛圆周上的情况"

    在圆周 $|z-z_0|=R$ 上，级数可能处处收敛、处处发散，或部分收敛部分发散，需要逐点讨论：

    - $\sum z^n$（$R=1$）：在 $|z|=1$ 上通项模为 $1$，处处发散；
    - $\sum\frac{z^n}{n^2}$（$R=1$）：在 $|z|=1$ 上 $\sum\frac1{n^2}$ 收敛，处处绝对收敛；
    - $\sum\frac{z^n}{n}$（$R=1$）：$z=1$ 时为调和级数，发散；$z=-1$ 时为交错级数，收敛。

### 4.2.3 收敛半径的求法

!!! theorem "定理 4.5（比值法与根值法）"

    对幂级数 $\sum c_n(z-z_0)^n$，若

    $$\lim_{n\to\infty}\left|\frac{c_{n+1}}{c_n}\right|=\lambda\qquad\text{或}\qquad\lim_{n\to\infty}\sqrt[n]{|c_n|}=\lambda,$$

    则收敛半径 $R=\dfrac1\lambda$（约定 $\lambda=0$ 时 $R=+\infty$，$\lambda=+\infty$ 时 $R=0$）。

!!! example "例题 4.3"

    求下列幂级数的收敛半径：

    (1) $\displaystyle\sum_{n=1}^\infty\frac{z^n}{n^2}$；(2) $\displaystyle\sum_{n=0}^\infty\frac{z^n}{n!}$；(3) $\displaystyle\sum_{n=1}^\infty(nz)^n$；(4) $\displaystyle\sum_{n=1}^\infty\left(1+\frac1n\right)^{n^2}(z-1)^n$；(5) $\displaystyle\sum_{n=0}^\infty\cos(in)\,z^n$。

    **解。**

    (1) $\left|\frac{c_{n+1}}{c_n}\right|=\frac{n^2}{(n+1)^2}\to1$，$R=1$，收敛圆为 $|z|<1$。

    (2) $\frac{n!}{(n+1)!}=\frac1{n+1}\to0$，$R=+\infty$。

    (3) $\sqrt[n]{n^n}=n\to\infty$，$R=0$，只在 $z=0$ 收敛。

    (4) $\sqrt[n]{|c_n|}=\left(1+\frac1n\right)^n\to e$，$R=\frac1e$，收敛圆为 $|z-1|<\frac1e$。

    (5) $\cos(in)=\cosh n=\frac{e^n+e^{-n}}{2}$，$\frac{\cosh(n+1)}{\cosh n}\to e$，$R=\frac1e$。

!!! warning "易错点：缺项幂级数"

    对 $\sum_{n=0}^\infty2^nz^{2n}$ 这类缺项级数，奇数次项的系数为 $0$，$\left|\frac{c_{n+1}}{c_n}\right|$ 没有意义。应当把它看成 $\zeta=z^2$ 的幂级数（$\sum2^n\zeta^n$ 的收敛半径为 $\frac12$，所以 $|z^2|<\frac12$，即 $R=\frac1{\sqrt2}$），或者直接对通项 $|2^nz^{2n}|$ 用比值法。

### 4.2.4 幂级数的性质

!!! theorem "定理 4.6（幂级数的运算与分析性质）"

    设 $\sum a_nz^n$ 的收敛半径为 $r_1$，$\sum b_nz^n$ 的收敛半径为 $r_2$，$R=\min(r_1,r_2)$。

    1. **四则运算**：在 $|z|<R$ 内，$\sum a_nz^n\pm\sum b_nz^n=\sum(a_n\pm b_n)z^n$，

        $$\Big(\sum a_nz^n\Big)\Big(\sum b_nz^n\Big)=\sum_{n=0}^\infty\big(a_0b_n+a_1b_{n-1}+\cdots+a_nb_0\big)z^n ;$$

    2. **代换（复合）**：若 $|z|<r$ 时 $\sum a_n z^n$ 收敛于 $f(z)$，$g(z)$ 在 $|z|<R$ 内解析且 $|g(z)|<r$，则 $f(g(z))=\sum a_n[g(z)]^n$；
    3. **和函数的分析性质**：设 $\sum c_n(z-z_0)^n$ 的收敛半径为 $R>0$，和函数为 $f(z)$。则 $f$ 在收敛圆内**解析**，并且可以**逐项求导、逐项积分**，所得级数的收敛半径仍为 $R$：

        $$f'(z)=\sum_{n=1}^\infty nc_n(z-z_0)^{n-1},\qquad \int_{z_0}^{z}f(\zeta)\,\mathrm d\zeta=\sum_{n=0}^\infty\frac{c_n}{n+1}(z-z_0)^{n+1}.$$

!!! example "例题 4.4（代换法的典型用法）"

    把 $\dfrac1{z-b}$ 展开成 $z-a$ 的幂级数（$a\ne b$）。

    **解。** 凑出“$\frac{1}{1-q}$”的形式：

    $$\frac1{z-b}=\frac{1}{(z-a)-(b-a)}=-\frac{1}{b-a}\cdot\frac{1}{1-\frac{z-a}{b-a}}=-\sum_{n=0}^\infty\frac{(z-a)^n}{(b-a)^{n+1}},$$

    在 $\left|\frac{z-a}{b-a}\right|<1$，即 $|z-a|<|b-a|$ 时成立。收敛半径正好是展开中心 $a$ 到奇点 $b$ 的距离。这是下一节泰勒定理的一个预演。

---

## 4.3 泰勒级数

### 4.3.1 泰勒展开定理

!!! theorem "定理 4.7（泰勒展开定理）"

    设 $f(z)$ 在区域 $D$ 内解析，$z_0\in D$，$R$ 为 $z_0$ 到 $D$ 的边界的距离（$D$ 为全平面时 $R=+\infty$）。则在 $|z-z_0|<R$ 内

    $$f(z)=\sum_{n=0}^\infty c_n(z-z_0)^n,\qquad c_n=\frac{f^{(n)}(z_0)}{n!}=\frac{1}{2\pi i}\oint_C\frac{f(\zeta)}{(\zeta-z_0)^{n+1}}\,\mathrm d\zeta,$$

    其中 $C$ 为圆周 $|\zeta-z_0|=\rho<R$。并且展开式是**唯一**的。

??? proof "证明思路"

    任取 $z$ 满足 $|z-z_0|<\rho<R$，由柯西积分公式 $f(z)=\frac{1}{2\pi i}\oint_C\frac{f(\zeta)}{\zeta-z}\,\mathrm d\zeta$。关键是把柯西核展开成几何级数：在 $C$ 上 $\left|\frac{z-z_0}{\zeta-z_0}\right|<1$，

    $$\frac{1}{\zeta-z}=\frac{1}{(\zeta-z_0)-(z-z_0)}=\frac{1}{\zeta-z_0}\cdot\frac{1}{1-\frac{z-z_0}{\zeta-z_0}}=\sum_{n=0}^\infty\frac{(z-z_0)^n}{(\zeta-z_0)^{n+1}} .$$

    代入后逐项积分（级数在 $C$ 上一致收敛），再用高阶导数公式（定理 3.10）得 $c_n=\frac{f^{(n)}(z_0)}{n!}$。唯一性：若 $f=\sum a_n(z-z_0)^n$，逐项求导后令 $z=z_0$，得 $a_n=\frac{f^{(n)}(z_0)}{n!}$。$\blacksquare$

!!! theorem "推论"

    1. $f$ 在 $z_0$ **解析** $\iff$ $f$ 在 $z_0$ 的某邻域内可以展开成幂级数（第 2 章的“可导”与这里的“可展开”是同一回事）；
    2. 若 $f$ 在 $z_0$ 解析，$z_1$ 是 $f$ 的**离 $z_0$ 最近的奇点**，则 $f$ 在 $z_0$ 的泰勒级数的收敛半径 $R=|z_1-z_0|$。

推论 2 的实用价值很大：**不用计算任何系数**，就能知道收敛半径。

<figure class="video" markdown>
<video controls preload="metadata" width="100%" poster="assets/posters/C4_RadiusOfConvergence.jpg">
<source src="assets/videos/C4_RadiusOfConvergence.mp4" type="video/mp4">
</video>
<figcaption>动画 4-2　$\frac1{1+x^2}$ 在实轴上处处光滑，但泰勒多项式只在 $|x|<1$ 内逼近它。原因在复平面上：奇点 $\pm i$ 限定了收敛圆的半径 $R=1$。换到 $z_0=1$ 展开，$R=|1-i|=\sqrt2$。</figcaption>
</figure>

!!! warning "易错点：实分析的直觉在这里会失效"

    $\frac{1}{1+x^2}$ 和 $e^x$ 在实轴上都处处光滑，但前者的麦克劳林级数在 $|x|\ge1$ 发散，后者处处收敛。原因只有在复平面上才看得到：$\frac{1}{1+z^2}$ 有奇点 $\pm i$，而 $e^z$ 没有奇点。

### 4.3.2 常用展开式

!!! success "必背的泰勒展开式（在 $z_0=0$ 处）"

    | 函数 | 展开式 | 收敛范围 |
    |---|---|---|
    | $\dfrac1{1-z}$ | $\displaystyle\sum_{n=0}^\infty z^n=1+z+z^2+\cdots$ | $\lvert z\rvert<1$ |
    | $\dfrac1{1+z}$ | $\displaystyle\sum_{n=0}^\infty(-1)^nz^n=1-z+z^2-\cdots$ | $\lvert z\rvert<1$ |
    | $e^z$ | $\displaystyle\sum_{n=0}^\infty\frac{z^n}{n!}=1+z+\frac{z^2}{2!}+\cdots$ | $\lvert z\rvert<\infty$ |
    | $\sin z$ | $\displaystyle\sum_{n=0}^\infty\frac{(-1)^nz^{2n+1}}{(2n+1)!}=z-\frac{z^3}{3!}+\cdots$ | $\lvert z\rvert<\infty$ |
    | $\cos z$ | $\displaystyle\sum_{n=0}^\infty\frac{(-1)^nz^{2n}}{(2n)!}=1-\frac{z^2}{2!}+\cdots$ | $\lvert z\rvert<\infty$ |
    | $\ln(1+z)$（主值） | $\displaystyle\sum_{n=0}^\infty\frac{(-1)^nz^{n+1}}{n+1}=z-\frac{z^2}{2}+\frac{z^3}{3}-\cdots$ | $\lvert z\rvert<1$ |
    | $(1+z)^\alpha$（主值） | $\displaystyle1+\alpha z+\frac{\alpha(\alpha-1)}{2!}z^2+\cdots$ | $\lvert z\rvert<1$ |

    收敛范围都可以由推论 2 直接读出：$\frac1{1\mp z}$ 的奇点为 $\pm1$；$\ln(1+z)$、$(1+z)^\alpha$ 的支点为 $-1$；$e^z,\sin z,\cos z$ 没有奇点。

### 4.3.3 展开方法

!!! success "两种方法"

    1. **直接法**：计算 $c_n=\frac{f^{(n)}(z_0)}{n!}$。只适合导数有规律的函数（$e^z,\sin z,\cos z$）；
    2. **间接法**（主要方法）：由唯一性，用任何方法得到的幂级数就是泰勒级数。利用已知展开式，通过**代换、四则运算、逐项求导、逐项积分**得到。最后一定要写出收敛范围（用推论 2）。

!!! example "例题 4.5（直接法）"

    把 $e^z$ 在 $z=0$ 展开。

    **解。** $(e^z)^{(n)}=e^z$，$c_n=\frac{e^0}{n!}=\frac1{n!}$，所以 $e^z=\sum_{n=0}^\infty\frac{z^n}{n!}$。$e^z$ 在全平面解析，$R=+\infty$。

!!! example "例题 4.6（逐项积分）"

    把 $f(z)=\ln(1+z)$（主值）分别在 $z=0$ 和 $z=1$ 处展开成幂级数。

    **解。** $\ln(1+z)$ 在 $\mathbb C\setminus(-\infty,-1]$ 内解析，且 $f'(z)=\frac1{1+z}$。

    (1) 在 $z=0$ 处：$f'(z)=\sum_{n=0}^\infty(-1)^nz^n$，$|z|<1$。从 $0$ 到 $z$ 逐项积分，并注意 $f(0)=0$：

    $$\ln(1+z)=\sum_{n=0}^\infty\frac{(-1)^n}{n+1}z^{n+1},\qquad|z|<1 .$$

    (2) 在 $z=1$ 处：$\dfrac1{1+z}=\dfrac{1}{2+(z-1)}=\dfrac12\cdot\dfrac{1}{1+\frac{z-1}{2}}=\displaystyle\sum_{n=0}^\infty\frac{(-1)^n(z-1)^n}{2^{n+1}}$，$|z-1|<2$。从 $1$ 到 $z$ 逐项积分，$f(1)=\ln2$：

    $$\ln(1+z)=\ln2+\sum_{n=0}^\infty\frac{(-1)^n}{(n+1)2^{n+1}}(z-1)^{n+1},\qquad|z-1|<2 .$$

    收敛半径 $2$ 正好是 $z=1$ 到支点 $-1$ 的距离。

!!! example "例题 4.7（三角恒等式 + 代换）"

    把 $\sin^2z$ 在 $z=0$ 处展开。

    **解。** $\sin^2z=\dfrac{1-\cos2z}{2}$，把 $\cos w$ 的展开式中的 $w$ 换成 $2z$：

    $$\sin^2z=\frac12\left[1-\sum_{n=0}^\infty\frac{(-1)^n(2z)^{2n}}{(2n)!}\right]=\sum_{n=1}^\infty\frac{(-1)^{n+1}2^{2n-1}}{(2n)!}z^{2n}=z^2-\frac{z^4}{3}+\cdots,\quad|z|<\infty .$$

!!! example "例题 4.8（换展开中心）"

    把 $\sin z$ 在 $z=1$ 处展开。

    **解。** 令 $\zeta=z-1$，$\sin z=\sin(1+\zeta)=\sin1\cos\zeta+\cos1\sin\zeta$：

    $$\sin z=\sin1\sum_{n=0}^\infty\frac{(-1)^n(z-1)^{2n}}{(2n)!}+\cos1\sum_{n=0}^\infty\frac{(-1)^n(z-1)^{2n+1}}{(2n+1)!},\qquad|z-1|<\infty .$$

!!! example "例题 4.9（逐项求导）"

    把 $\dfrac{1}{(1+z)^2}$ 在 $z=0$ 处展开。

    **解。** $\dfrac1{(1+z)^2}=-\left(\dfrac1{1+z}\right)'=-\displaystyle\sum_{n=1}^\infty(-1)^nnz^{n-1}=\sum_{n=0}^\infty(-1)^n(n+1)z^n$，$|z|<1$。

!!! example "例题 4.10（有理函数：部分分式 + 几何级数）"

    把 $f(z)=\dfrac{1}{(z-1)(z-2)}$ 在 $z=0$ 处展开。

    **解。** $f(z)=\dfrac1{1-z}-\dfrac1{2-z}$（通分验证：$\frac{(2-z)-(1-z)}{(1-z)(2-z)}=\frac{1}{(z-1)(z-2)}$）。

    $$f(z)=\sum_{n=0}^\infty z^n-\frac12\sum_{n=0}^\infty\frac{z^n}{2^n}=\sum_{n=0}^\infty\left(1-\frac{1}{2^{n+1}}\right)z^n,\qquad|z|<1 .$$

    最近的奇点是 $z=1$，$R=1$。

---

## 4.4 洛朗级数

### 4.4.1 含负幂项的级数

$\frac1{1-z}$ 只有一个奇点 $z=1$，但它的泰勒级数在 $|z|>1$ 上就发散了。如果非要在 $|z|>1$ 上把它写成级数，可以换一个方向凑几何级数：

$$\frac{1}{1-z}=-\frac1z\cdot\frac{1}{1-\frac1z}=-\sum_{n=1}^\infty\frac{1}{z^n},\qquad\left|\frac1z\right|<1\iff|z|>1 .$$

这里出现的是 $z$ 的**负幂**。一般地，考虑**双边幂级数**

$$\sum_{n=-\infty}^{+\infty}c_n(z-z_0)^n=\underbrace{\sum_{n=1}^{\infty}c_{-n}(z-z_0)^{-n}}_{\text{负幂部分}}+\underbrace{\sum_{n=0}^{\infty}c_n(z-z_0)^n}_{\text{正幂部分}} .$$

正幂部分在某个圆 $|z-z_0|<R_2$ 内收敛；负幂部分是 $\zeta=\frac1{z-z_0}$ 的幂级数，在 $|z-z_0|>R_1$ 内收敛。当 $R_1<R_2$ 时，双边幂级数在**圆环** $R_1<|z-z_0|<R_2$ 内收敛，和函数在圆环内解析。

### 4.4.2 洛朗定理

!!! theorem "定理 4.8（洛朗定理）"

    设 $f(z)$ 在圆环 $D:R_1<|z-z_0|<R_2$（$0\le R_1<R_2\le+\infty$）内解析，则在 $D$ 内

    $$f(z)=\sum_{n=-\infty}^{+\infty}c_n(z-z_0)^n,\qquad c_n=\frac{1}{2\pi i}\oint_C\frac{f(\zeta)}{(\zeta-z_0)^{n+1}}\,\mathrm d\zeta\quad(n\in\mathbb Z),$$

    其中 $C$ 为圆环内绕 $z_0$ 的任意一条正向简单闭曲线。该展开式（**洛朗级数**）是**唯一**的。

??? proof "证明思路"

    取 $R_1<r<|z-z_0|<\rho<R_2$，作圆周 $\Gamma_1:|\zeta-z_0|=r$，$\Gamma_2:|\zeta-z_0|=\rho$。在二连通区域上用柯西积分公式（复合闭路）：

    $$f(z)=\frac{1}{2\pi i}\oint_{\Gamma_2}\frac{f(\zeta)}{\zeta-z}\,\mathrm d\zeta-\frac{1}{2\pi i}\oint_{\Gamma_1}\frac{f(\zeta)}{\zeta-z}\,\mathrm d\zeta .$$

    在 $\Gamma_2$ 上 $\left|\frac{z-z_0}{\zeta-z_0}\right|<1$，按泰勒定理那样展开得到**正幂**部分。在 $\Gamma_1$ 上 $\left|\frac{\zeta-z_0}{z-z_0}\right|<1$，于是

    $$-\frac{1}{\zeta-z}=\frac{1}{(z-z_0)-(\zeta-z_0)}=\sum_{n=1}^\infty\frac{(\zeta-z_0)^{n-1}}{(z-z_0)^{n}},$$

    逐项积分得到**负幂**部分。最后由闭路变形把 $\Gamma_1,\Gamma_2$ 都换成 $C$。$\blacksquare$

!!! warning "易错点：洛朗系数不是导数"

    当 $f$ 在 $z_0$ 不解析时，$c_n=\frac{1}{2\pi i}\oint_C\frac{f(\zeta)}{(\zeta-z_0)^{n+1}}\mathrm d\zeta$ **不等于** $\frac{f^{(n)}(z_0)}{n!}$（$f^{(n)}(z_0)$ 根本不存在）。只有当 $f$ 在整个圆盘 $|z-z_0|<R_2$ 内解析时，洛朗级数才退化为泰勒级数。

!!! theorem "最重要的一个系数：$c_{-1}$"

    在洛朗系数公式中取 $n=-1$：

    $$\oint_C f(z)\,\mathrm dz=2\pi i\,c_{-1} .$$

    圆环内的闭路积分只由 $(z-z_0)^{-1}$ 项的系数决定。第 5 章把 $c_{-1}$ 称为**留数**。这其实是例题 3.2 的推广：逐项积分时，只有 $n=-1$ 那一项给出 $2\pi i$，其余各项都给出 $0$。

### 4.4.3 洛朗展开的方法

!!! success "标准步骤（间接展开法）"

    1. **找奇点、分圆环**：以展开中心 $z_0$ 为圆心，过各奇点作圆，把平面分成若干个解析圆环（若题目已指定圆环，就只做那一个）；
    2. **部分分式**：把有理函数分解成 $\frac{A}{z-a}$、$\frac{B}{(z-a)^2}$ 等简单项；
    3. **逐项凑几何级数**：对每个 $\frac{1}{z-a}$，比较 $|z-z_0|$ 与 $|a-z_0|$ 的大小：
        - 在圆环内若 $|z-z_0|<|a-z_0|$，提出 $(a-z_0)$，得到**正幂**：$\dfrac{1}{z-a}=-\dfrac{1}{a-z_0}\cdot\dfrac{1}{1-\frac{z-z_0}{a-z_0}}$；
        - 若 $|z-z_0|>|a-z_0|$，提出 $(z-z_0)$，得到**负幂**：$\dfrac{1}{z-a}=\dfrac{1}{z-z_0}\cdot\dfrac{1}{1-\frac{a-z_0}{z-z_0}}$；
    4. 非有理函数（如 $e^{1/z}$、$\sin\frac1z$）直接把 $\frac1{z-z_0}$ 代入已知的泰勒展开式；
    5. 合并，写出成立的圆环。

    一句话：**保证公比的模小于 1**。

!!! example "例题 4.11（代入已知展开式）"

    把 (1) $z^3e^{1/z}$，(2) $\dfrac{e^z}{z^2}$ 在 $0<|z|<\infty$ 内展开成洛朗级数。

    **解。** (1) $e^{1/z}=\sum_{n=0}^\infty\frac{1}{n!\,z^n}$（$z\ne0$），所以

    $$z^3e^{1/z}=\sum_{n=0}^\infty\frac{z^{3-n}}{n!}=z^3+z^2+\frac{z}{2!}+\frac{1}{3!}+\frac{1}{4!\,z}+\frac{1}{5!\,z^2}+\cdots .$$

    负幂项有无穷多个。$c_{-1}=\frac1{4!}=\frac1{24}$，所以 $\oint_{|z|=1}z^3e^{1/z}\,\mathrm dz=\frac{2\pi i}{24}=\frac{\pi i}{12}$。

    (2) $\dfrac{e^z}{z^2}=\dfrac{1}{z^2}\left(1+z+\dfrac{z^2}{2!}+\dfrac{z^3}{3!}+\cdots\right)=\dfrac1{z^2}+\dfrac1z+\dfrac1{2!}+\dfrac{z}{3!}+\cdots$。负幂项只有两项，$c_{-1}=1$。

!!! example "例题 4.12（同一函数，不同圆环）"

    把 $f(z)=\dfrac{1}{(z-1)(z-2)}$ 分别在 (1) $|z|<1$；(2) $1<|z|<2$；(3) $2<|z|<\infty$；(4) $0<|z-1|<1$；(5) $1<|z-1|<\infty$ 内展开成洛朗级数。

    **解。** 奇点为 $1,2$。部分分式 $f(z)=\dfrac{1}{1-z}-\dfrac{1}{2-z}$。

    (1) $|z|<1$：$|z|<1$ 且 $\left|\frac z2\right|<1$，两项都给正幂，得到例题 4.10 的泰勒级数：

    $$f(z)=\sum_{n=0}^\infty\left(1-\frac1{2^{n+1}}\right)z^n .$$

    (2) $1<|z|<2$：$\left|\frac1z\right|<1$，$\left|\frac z2\right|<1$：

    $$\frac1{1-z}=-\frac1z\cdot\frac1{1-\frac1z}=-\sum_{n=1}^\infty\frac{1}{z^n},\qquad\frac1{2-z}=\frac12\cdot\frac{1}{1-\frac z2}=\sum_{n=0}^\infty\frac{z^n}{2^{n+1}},$$

    $$f(z)=-\sum_{n=1}^\infty\frac{1}{z^n}-\sum_{n=0}^\infty\frac{z^n}{2^{n+1}} .$$

    (3) $|z|>2$：$\left|\frac1z\right|<1$，$\left|\frac2z\right|<1$：

    $$\frac1{2-z}=-\frac1z\cdot\frac1{1-\frac2z}=-\sum_{n=1}^\infty\frac{2^{n-1}}{z^n},\qquad f(z)=\sum_{n=1}^\infty\frac{2^{n-1}-1}{z^n}=\sum_{n=2}^\infty\frac{2^{n-1}-1}{z^n} .$$

    (4) 以 $z=1$ 为中心，令 $w=z-1$，$f=\dfrac{1}{w(w-1)}$。在 $0<|w|<1$ 内：

    $$f(z)=-\frac1w\cdot\frac{1}{1-w}=-\sum_{n=0}^\infty w^{n-1}=-\frac{1}{z-1}-1-(z-1)-(z-1)^2-\cdots .$$

    (5) 在 $|w|>1$ 内：$\dfrac{1}{w-1}=\dfrac1w\cdot\dfrac{1}{1-\frac1w}=\displaystyle\sum_{n=1}^\infty\frac{1}{w^n}$，所以 $f(z)=\displaystyle\sum_{n=2}^\infty\frac{1}{(z-1)^n}$。

    **观察。** (3) 和 (5) 的展开式只有负幂，但这并不说明展开中心是奇点：(3) 的中心 $z=0$ 是 $f$ 的解析点。圆环内展开出现负幂，只说明圆环**内边界以内**有奇点。只有 (4) 这种“去心邻域”型圆环中的负幂部分，才刻画奇点 $z=1$ 本身（第 5 章的奇点分类）。

<figure class="video" markdown>
<video controls preload="metadata" width="100%" poster="assets/posters/C4_LaurentAnnuli.jpg">
<source src="assets/videos/C4_LaurentAnnuli.mp4" type="video/mp4">
</video>
<figcaption>动画 4-3　奇点 $1,2$ 把平面切成三个以原点为心的圆环，每个圆环内的展开式各不相同；以奇点 $z=1$ 为心的去心邻域中，$(z-1)^{-1}$ 的系数 $c_{-1}=-1$ 就是留数。</figcaption>
</figure>

!!! example "例题 4.13（以奇点为中心）"

    把 $f(z)=\dfrac{1}{z^2+1}$ 在以 $z=i$ 为中心的圆环内展开：(1) $0<|z-i|<2$；(2) $|z-i|>2$。

    **解。** 令 $w=z-i$，$f=\dfrac{1}{(z-i)(z+i)}=\dfrac{1}{w(w+2i)}$。不必做部分分式，$\frac1w$ 本身已经是 $w$ 的幂。

    (1) $|w|<2=|2i|$：$\dfrac1{w+2i}=\dfrac{1}{2i}\cdot\dfrac{1}{1+\frac{w}{2i}}=\displaystyle\sum_{n=0}^\infty\frac{(-1)^nw^n}{(2i)^{n+1}}$，所以

    $$f(z)=\sum_{n=0}^\infty\frac{(-1)^n}{(2i)^{n+1}}(z-i)^{n-1},\qquad c_{-1}=\frac1{2i}.$$

    (2) $|w|>2$：$\dfrac1{w+2i}=\dfrac1w\cdot\dfrac{1}{1+\frac{2i}{w}}=\displaystyle\sum_{n=0}^\infty\frac{(-1)^n(2i)^n}{w^{n+1}}$，所以

    $$f(z)=\sum_{n=0}^\infty\frac{(-1)^n(2i)^n}{(z-i)^{n+2}} .$$

!!! example "例题 4.14（含二阶因子）"

    把 $f(z)=\dfrac{z}{(z-1)^2(2-z)}$ 在 (1) $0<|z-1|<1$；(2) $1<|z-1|<\infty$ 内展开。

    **解。** 令 $w=z-1$，则 $z=1+w$，$2-z=1-w$：

    $$f(z)=\frac{1+w}{w^2}\cdot\frac{1}{1-w}=\left(\frac{1}{w^2}+\frac1w\right)\frac{1}{1-w}.$$

    (1) $|w|<1$：$\frac1{1-w}=\sum_{n\ge0}w^n$，

    $$f(z)=\sum_{n=0}^\infty w^{n-2}+\sum_{n=0}^\infty w^{n-1}=\frac{1}{(z-1)^2}+2\sum_{n=-1}^\infty(z-1)^n .$$

    (2) $|w|>1$：$\frac1{1-w}=-\sum_{n\ge1}w^{-n}$，

    $$f(z)=-\sum_{n=1}^\infty w^{-n-2}-\sum_{n=1}^\infty w^{-n-1}=-\frac{1}{(z-1)^2}-2\sum_{n=3}^\infty\frac{1}{(z-1)^n} .$$

!!! example "例题 4.15（用 $c_{-1}$ 计算积分）"

    计算 $\displaystyle\oint_{|z|=3/2}\frac{\mathrm dz}{(z-1)(z-2)}$ 与 $\displaystyle\oint_{|z|=3}\frac{\mathrm dz}{(z-1)(z-2)}$。

    **解。** $|z|=\frac32$ 位于圆环 $1<|z|<2$ 内。由例题 4.12(2)，$z^{-1}$ 的系数 $c_{-1}=-1$，故积分为 $2\pi i\cdot(-1)=-2\pi i$。

    $|z|=3$ 位于圆环 $|z|>2$ 内。由例题 4.12(3)，展开式从 $z^{-2}$ 开始，$c_{-1}=0$，积分为 $0$。

    （用第 3 章的方法验证：第一个积分只包含奇点 $z=1$，$2\pi i\cdot\frac{1}{1-2}=-2\pi i$；第二个积分为 $2\pi i(-1+1)=0$。两种方法结果一致。）

!!! warning "易错点"

    - 展开前**必须先确定圆环**，同一个函数在不同圆环内的展开式不同，答案一定要注明成立范围；
    - 凑几何级数时“公比”的模必须小于 $1$。在 $1<|z|<2$ 内把 $\frac1{1-z}$ 写成 $\sum z^n$ 是最常见的错误；
    - 以 $z_0$ 为中心展开，就必须写成 $(z-z_0)$ 的幂，不能写成 $z$ 的幂；
    - 用 $\oint_Cf\,\mathrm dz=2\pi ic_{-1}$ 时，所用的必须是**包含 $C$ 的那个圆环**内的展开式。

---

## 知识点归纳总结

!!! abstract "知识点归纳总结"

    ```mermaid
    flowchart TD
        SEQ["复数序列 / 级数<br/>⇔ 实部、虚部分别收敛"] --> ABS["绝对收敛 Σ|zₙ|<br/>正项级数判别法"]
        ABS --> PS["幂级数 Σ cₙ(z−z₀)ⁿ"]
        PS --> ABEL["阿贝尔定理 → 收敛圆 |z−z₀| < R"]
        ABEL --> RAD["R = 1/lim|cₙ₊₁/cₙ| = 1/lim ⁿ√|cₙ|"]
        PS --> PROP["和函数解析<br/>逐项求导 / 积分"]
        HD["第3章：高阶导数公式"] --> TAY["泰勒定理<br/>cₙ = f⁽ⁿ⁾(z₀)/n!"]
        TAY --> NEAR["R = 到最近奇点的距离"]
        TAY --> STD["e^z, sin z, cos z,<br/>1/(1−z), ln(1+z), (1+z)^α"]
        PROP --> TAY
        TAY --> LAU["洛朗定理（圆环内）<br/>正幂 + 负幂"]
        LAU --> METH["分圆环 → 部分分式<br/>→ 公比模 < 1"]
        LAU --> CM1["∮ f dz = 2πi c₋₁"]
        CM1 -.第5章.-> RES["留数定理"]
        LAU -.第5章.-> SING["孤立奇点分类"]
    ```

    **核心结论**

    | 内容 | 结论 |
    |---|---|
    | 序列 / 级数 | $z_n\to z_0\iff x_n\to x_0,\ y_n\to y_0$；$\sum z_n$ 收敛 ⇒ $z_n\to0$ |
    | 绝对收敛 | $\sum\lvert z_n\rvert$ 收敛 ⇒ $\sum z_n$ 收敛 |
    | 阿贝尔定理 | 在 $z_1$ 收敛 ⇒ $\lvert z\rvert<\lvert z_1\rvert$ 内绝对收敛 |
    | 收敛半径 | $R=1/\lambda$，$\lambda=\lim\lvert c_{n+1}/c_n\rvert$ 或 $\lim\sqrt[n]{\lvert c_n\rvert}$ |
    | 泰勒 | $c_n=\frac{f^{(n)}(z_0)}{n!}$，$R=$ 到最近奇点的距离 |
    | 洛朗 | $c_n=\frac{1}{2\pi i}\oint_C\frac{f(\zeta)}{(\zeta-z_0)^{n+1}}\mathrm d\zeta$，$n\in\mathbb Z$，圆环内唯一 |
    | 关键系数 | $\oint_Cf\,\mathrm dz=2\pi i\,c_{-1}$ |

    **方法清单**

    1. 判断复级数收敛：先看通项是否趋于 $0$ → 看 $\sum|z_n|$ → 不行再拆实部、虚部；
    2. 求收敛半径：比值法 / 根值法；缺项级数换元；已知和函数时直接找最近奇点；
    3. 泰勒展开：间接法（代换、四则、逐项求导或积分），并写出收敛圆；
    4. 洛朗展开：分圆环 → 部分分式 → 每项凑 $\frac1{1-q}$ 且 $|q|<1$ → 合并，并写出圆环；
    5. 用 $c_{-1}$ 求闭路积分：选包含积分路径的那个圆环内的展开式。

    **易错点汇总**

    - $|z_n|$ 收敛 ⇏ $z_n$ 收敛；
    - 收敛圆周上要逐点讨论；
    - 缺项级数不能直接用 $|c_{n+1}/c_n|$；
    - 泰勒级数要写收敛范围；实轴上光滑不代表收敛半径为无穷；
    - 洛朗展开中公比的模必须小于 $1$，展开式必须与圆环对应；
    - 洛朗系数 $c_n\ne f^{(n)}(z_0)/n!$。

---

## 练习题

!!! example "练习 1"

    判断收敛性，收敛时求和：(1) $\displaystyle\sum_{n=1}^\infty\left(\frac{1+i}{2}\right)^n$；(2) $\displaystyle\sum_{n=1}^\infty\left(\frac1n+\frac{i}{n^2}\right)$。

??? success "解答"

    (1) 公比 $q=\frac{1+i}{2}$，$|q|=\frac{\sqrt2}{2}<1$，绝对收敛，

    $$\sum_{n=1}^\infty q^n=\frac{q}{1-q}=\frac{(1+i)/2}{(1-i)/2}=\frac{1+i}{1-i}=i .$$

    (2) 实部 $\sum\frac1n$ 发散，原级数发散。

!!! example "练习 2"

    求收敛半径：(1) $\displaystyle\sum_{n=1}^\infty\frac{n!}{n^n}z^n$；(2) $\displaystyle\sum_{n=0}^\infty2^nz^{2n}$；(3) $\displaystyle\sum_{n=1}^\infty\frac{(z-i)^n}{n\,3^n}$。

??? success "解答"

    (1) $\left|\dfrac{c_{n+1}}{c_n}\right|=\dfrac{(n+1)!\,n^n}{(n+1)^{n+1}\,n!}=\left(\dfrac{n}{n+1}\right)^n\to\dfrac1e$，$R=e$。

    (2) 看成 $\zeta=z^2$ 的级数 $\sum2^n\zeta^n$，收敛半径 $\frac12$，即 $|z|^2<\frac12$，$R=\frac{1}{\sqrt2}$。

    (3) $\sqrt[n]{\frac{1}{n3^n}}\to\frac13$，$R=3$，收敛圆为 $|z-i|<3$。

!!! example "练习 3"

    把 $f(z)=\dfrac{z}{(z+1)(z+2)}$ 在 $z=2$ 处展开成泰勒级数，并指出收敛半径。

??? success "解答"

    部分分式：$f(z)=-\dfrac{1}{z+1}+\dfrac{2}{z+2}$。令 $w=z-2$：

    $$\frac1{z+1}=\frac{1}{3+w}=\frac13\sum_{n=0}^\infty\left(-\frac w3\right)^n\ (|w|<3),\qquad \frac1{z+2}=\frac1{4+w}=\frac14\sum_{n=0}^\infty\left(-\frac w4\right)^n\ (|w|<4).$$

    $$f(z)=\sum_{n=0}^\infty(-1)^n\left(\frac{2}{4^{n+1}}-\frac{1}{3^{n+1}}\right)(z-2)^n,\qquad |z-2|<3 .$$

    $R=3$ 是 $z=2$ 到最近奇点 $-1$ 的距离。

!!! example "练习 4"

    求 $\dfrac{e^z}{1-z}$ 在 $z=0$ 处泰勒展开的前四项，并指出收敛半径。

??? success "解答"

    柯西乘积：$\left(1+z+\frac{z^2}{2}+\frac{z^3}{6}+\cdots\right)(1+z+z^2+z^3+\cdots)$，$z^n$ 的系数为 $\sum_{k=0}^n\frac1{k!}$：

    $$\frac{e^z}{1-z}=1+2z+\frac52z^2+\frac83z^3+\cdots,\qquad R=1 .$$

!!! example "练习 5"

    把 $\arctan z=\displaystyle\int_0^z\frac{\mathrm d\zeta}{1+\zeta^2}$ 在 $z=0$ 处展开成幂级数。

??? success "解答"

    $\dfrac1{1+\zeta^2}=\displaystyle\sum_{n=0}^\infty(-1)^n\zeta^{2n}$（$|\zeta|<1$），逐项积分：

    $$\arctan z=\sum_{n=0}^\infty\frac{(-1)^n}{2n+1}z^{2n+1}=z-\frac{z^3}{3}+\frac{z^5}{5}-\cdots,\qquad|z|<1 .$$

    收敛半径为 $1$，因为被积函数有奇点 $\pm i$。

!!! example "练习 6"

    把 $f(z)=\dfrac{1}{z(z-1)}$ 分别在 (1) $0<|z|<1$；(2) $|z|>1$ 内展开成洛朗级数，并计算 $\displaystyle\oint_{|z|=2}f(z)\,\mathrm dz$。

??? success "解答"

    (1) $f=-\dfrac1z\cdot\dfrac{1}{1-z}=-\displaystyle\sum_{n=0}^\infty z^{n-1}=-\frac1z-1-z-z^2-\cdots$。

    (2) $f=\dfrac{1}{z^2}\cdot\dfrac{1}{1-\frac1z}=\displaystyle\sum_{n=0}^\infty\frac{1}{z^{n+2}}=\frac1{z^2}+\frac1{z^3}+\cdots$。

    $|z|=2$ 在圆环 $|z|>1$ 内，$c_{-1}=0$，积分为 $0$。（验证：$\frac{1}{z(z-1)}=\frac{1}{z-1}-\frac1z$，$2\pi i-2\pi i=0$。）

!!! example "练习 7"

    把 $f(z)=\dfrac{1}{z^2(z-i)}$ 在 $0<|z-i|<1$ 内展开成洛朗级数，写出 $c_{-1}$。

??? success "解答"

    令 $w=z-i$，$f=\dfrac{1}{w}\cdot\dfrac{1}{(w+i)^2}$。先展开

    $$\frac{1}{w+i}=\frac1i\cdot\frac{1}{1+\frac wi}=\sum_{n=0}^\infty\frac{(-1)^nw^n}{i^{n+1}}\quad(|w|<1),$$

    求导再取负号：$\dfrac{1}{(w+i)^2}=-\dfrac{\mathrm d}{\mathrm dw}\dfrac1{w+i}=\displaystyle\sum_{n=1}^\infty\frac{(-1)^{n+1}n\,w^{n-1}}{i^{n+1}}$。所以

    $$f(z)=\sum_{n=1}^\infty\frac{(-1)^{n+1}n}{i^{n+1}}(z-i)^{n-2}=-\frac{1}{z-i}+\cdots .$$

    $n=1$ 项给出 $c_{-1}=\frac{1}{i^2}=-1$。（验证：用柯西积分公式，在 $z=i$ 附近 $\frac{1}{z^2}\big|_{z=i}=-1$。）

!!! example "练习 8"

    计算 $\displaystyle\oint_{|z|=1}z^2\sin\frac1z\,\mathrm dz$。

??? success "解答"

    $z=0$ 是 $\sin\frac1z$ 的奇点，被积函数不能写成“$\frac{\text{解析函数}}{(z-z_0)^{n+1}}$”的形式，第 3 章的公式用不上。在 $0<|z|<\infty$ 内展开：

    $$z^2\sin\frac1z=z^2\left(\frac1z-\frac{1}{3!\,z^3}+\frac{1}{5!\,z^5}-\cdots\right)=z-\frac{1}{6z}+\frac{1}{120z^3}-\cdots .$$

    $c_{-1}=-\frac16$，积分为 $2\pi i\cdot\left(-\frac16\right)=-\frac{\pi i}{3}$。

---

**下一章预告。** 洛朗级数的负幂部分刻画了函数在孤立奇点附近的行为：负幂项没有、只有有限项、有无穷多项，分别对应可去奇点、极点、本性奇点。$c_{-1}$（留数）则把“算闭路积分”变成“算一个系数”，还能用来计算实积分，比如 $\int_{-\infty}^{\infty}\frac{\mathrm dx}{1+x^2}$。→ [第 5 章 留数及其应用](ch5.md)
