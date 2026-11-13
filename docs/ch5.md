# 第 5 章 留数及其应用

!!! tip "本章目标"

    1. 会求解析函数**零点的阶数**；会用洛朗级数、极限、零点阶数三种方法判别**孤立奇点的类型**（可去奇点、$m$ 阶极点、本性奇点），包括 $z=\infty$；
    2. 理解**留数** $\operatorname{Res}[f,z_0]=c_{-1}$，熟练掌握留数的计算规则（一阶极点、$\frac PQ$ 型、$m$ 阶极点、洛朗展开法、无穷远点的留数）；
    3. 熟练运用**留数定理**计算闭路积分——期末必考；
    4. 会用留数计算三类实积分：$\displaystyle\int_0^{2\pi}R(\cos\theta,\sin\theta)\,\mathrm d\theta$、$\displaystyle\int_{-\infty}^{+\infty}R(x)\,\mathrm dx$、$\displaystyle\int_{-\infty}^{+\infty}R(x)e^{iax}\,\mathrm dx$——期末必考。

**本章主线。** 第 3 章说闭路积分只取决于回路围住的奇点，第 4 章说在奇点附近要用洛朗级数，并且 $\oint_Cf\,\mathrm dz=2\pi i\,c_{-1}$。本章把这两件事合在一起：先按洛朗级数的负幂部分给奇点**分类**，再给出不展开级数也能求 $c_{-1}$（**留数**）的公式，最后把闭路积分归结为“求几个留数再相加”。这个方法还能用来计算许多用实分析很难算的定积分。

---

## 5.1 孤立奇点

### 5.1.1 零点

!!! definition "定义 5.1（零点及其阶数）"

    设 $f(z)$ 在 $z_0$ 处解析。若 $f(z_0)=0$，则称 $z_0$ 为 $f$ 的**零点**。若 $f$ 在 $z_0$ 的某邻域内可表示为

    $$f(z)=(z-z_0)^m\varphi(z),$$

    其中 $\varphi(z)$ 在 $z_0$ 处解析且 $\varphi(z_0)\ne0$，$m$ 为正整数，则称 $z_0$ 为 $f$ 的 **$m$ 阶零点**。

!!! theorem "定理 5.1（零点阶数的判别）"

    设 $f(z)$ 在 $z_0$ 处解析，则下列条件等价：

    1. $z_0$ 为 $f$ 的 $m$ 阶零点；
    2. $f(z_0)=f'(z_0)=\cdots=f^{(m-1)}(z_0)=0$，且 $f^{(m)}(z_0)\ne0$；
    3. $f$ 在 $z_0$ 处的泰勒级数从 $(z-z_0)^m$ 项开始，即 $c_0=\cdots=c_{m-1}=0$，$c_m\ne0$。

（由泰勒系数 $c_n=\frac{f^{(n)}(z_0)}{n!}$ 立即可得。）

!!! example "例题 5.1"

    求零点及其阶数：(1) $f(z)=z^3-1$；(2) $f(z)=z\sin z$ 在 $z=0$ 处；(3) $f(z)=z-\sin z$ 在 $z=0$ 处。

    **解。** (1) 零点为 $z_k=e^{2k\pi i/3}$（$k=0,1,2$）。$f'(z_k)=3z_k^2\ne0$，都是一阶零点。

    (2) $f(0)=0$，$f'(z)=\sin z+z\cos z$，$f'(0)=0$；$f''(z)=2\cos z-z\sin z$，$f''(0)=2\ne0$。所以 $z=0$ 为二阶零点。（也可以看级数：$z\sin z=z^2-\frac{z^4}{3!}+\cdots$。）

    (3) $z-\sin z=\frac{z^3}{3!}-\frac{z^5}{5!}+\cdots$，从 $z^3$ 开始，所以 $z=0$ 为三阶零点。

!!! theorem "定理 5.2（零点的孤立性）"

    若 $f$ 在 $z_0$ 处解析、$f(z_0)=0$，且 $f$ 在 $z_0$ 的邻域内不恒为零，则存在 $z_0$ 的一个去心邻域，在其中 $f(z)\ne0$。

??? proof "证明"

    由定理 5.1，$f(z)=(z-z_0)^m\varphi(z)$，$\varphi(z_0)\ne0$。$\varphi$ 在 $z_0$ 连续，取 $\varepsilon=\frac{|\varphi(z_0)|}{2}$，存在 $\delta>0$，当 $|z-z_0|<\delta$ 时 $|\varphi(z)-\varphi(z_0)|<\varepsilon$，从而 $|\varphi(z)|>\frac{|\varphi(z_0)|}{2}>0$。于是当 $0<|z-z_0|<\delta$ 时 $f(z)\ne0$。$\blacksquare$

### 5.1.2 孤立奇点及其分类

!!! definition "定义 5.2（孤立奇点）"

    若 $z_0$ 为 $f(z)$ 的奇点，且存在 $\delta>0$，使 $f$ 在去心邻域 $0<|z-z_0|<\delta$ 内处处解析，则称 $z_0$ 为 $f$ 的**孤立奇点**。

**不是所有奇点都是孤立的。** $\ln z$ 的奇点是原点和整条负实轴，它们都不孤立。$\frac{1}{\sin(1/z)}$ 的奇点为 $z=0$ 和 $z_k=\frac1{k\pi}$（$k=\pm1,\pm2,\dots$），$z_k\to0$，所以 $z=0$ 不是孤立奇点（它的任意去心邻域里都有别的奇点）。

在孤立奇点 $z_0$ 的去心邻域内，$f$ 可以展开为洛朗级数（定理 4.8）。按负幂部分的项数分类：

!!! definition "定义 5.3（孤立奇点的分类）"

    设 $z_0$ 为 $f$ 的孤立奇点，$f$ 在 $0<|z-z_0|<\delta$ 内的洛朗级数为 $\sum_{n=-\infty}^{+\infty}c_n(z-z_0)^n$。

    1. 若没有负幂项，称 $z_0$ 为**可去奇点**；
    2. 若只有有限个负幂项，且最低次项为 $c_{-m}(z-z_0)^{-m}$（$c_{-m}\ne0$），称 $z_0$ 为 **$m$ 阶极点**（一阶极点也称**简单极点**）；
    3. 若有无穷多个负幂项，称 $z_0$ 为**本性奇点**。

!!! theorem "定理 5.3（用极限判别）"

    设 $z_0$ 为 $f$ 的孤立奇点，则

    1. $z_0$ 为可去奇点 $\iff\lim\limits_{z\to z_0}f(z)$ 存在且有限；
    2. $z_0$ 为极点 $\iff\lim\limits_{z\to z_0}f(z)=\infty$；
    3. $z_0$ 为本性奇点 $\iff\lim\limits_{z\to z_0}f(z)$ 既不存在也不为 $\infty$。

    并且，$z_0$ 为 $m$ 阶极点 $\iff$ $f(z)=\dfrac{g(z)}{(z-z_0)^m}$，其中 $g$ 在 $z_0$ 处解析且 $g(z_0)\ne0$。

**可去奇点**的名字来源：令 $f(z_0)=c_0=\lim\limits_{z\to z_0}f(z)$，$f$ 就在 $z_0$ 处解析了，奇点被“去掉”。例如 $\frac{\sin z}{z}$，令其在 $0$ 处的值为 $1$。

**本性奇点**附近的情况最复杂。魏尔斯特拉斯定理说，$f$ 在本性奇点的任意小去心邻域内的值在复平面上稠密；皮卡定理更进一步：除了至多一个例外值以外，$f$ 在那里取遍所有复数无穷多次。例如 $e^{1/z}$ 在 $z=0$ 附近取遍除 $0$ 以外的所有值。

<figure class="video" markdown>
<video controls preload="metadata" width="100%" poster="assets/posters/C5_DomainColoring.jpg">
<source src="assets/videos/C5_DomainColoring.mp4" type="video/mp4">
</video>
<figcaption>动画 5-1　域着色（颜色 = 辐角，亮度 = 模）：$\frac{\sin z}{z}$ 在原点颜色平滑（可去奇点）；$\frac1{z^2}$ 在原点变白、颜色绕两圈（二阶极点）；$e^{1/z}$ 在原点附近越放大越混乱（本性奇点）。</figcaption>
</figure>

### 5.1.3 用零点判别极点的阶数

!!! theorem "定理 5.4（零点与极点的关系）"

    1. $z_0$ 为 $f$ 的 $m$ 阶极点 $\iff$ $z_0$ 为 $\dfrac1{f}$ 的 $m$ 阶零点（约定 $\frac1f$ 在 $z_0$ 处的值为 $0$）；
    2. 设 $f(z)=\dfrac{P(z)}{Q(z)}$，$P,Q$ 在 $z_0$ 处解析，$z_0$ 为 $P$ 的 $n$ 阶零点（$n=0$ 表示 $P(z_0)\ne0$）、$Q$ 的 $m$ 阶零点。当 $m>n$ 时，$z_0$ 为 $f$ 的 $m-n$ 阶极点；当 $m\le n$ 时，$z_0$ 为 $f$ 的可去奇点。

（第 2 条：$f=\frac{(z-z_0)^n\varphi(z)}{(z-z_0)^m\psi(z)}=\frac{1}{(z-z_0)^{m-n}}\cdot\frac{\varphi(z)}{\psi(z)}$，再用定理 5.3。）

!!! example "例题 5.2（教材例题改编）"

    判断下列函数的孤立奇点类型（极点要写出阶数）：

    (1) $\dfrac{1}{z(z^2+1)^2}$；(2) $\dfrac{1}{\cos z}$；(3) $\dfrac{\cos z}{\sin^2z}$；(4) $\dfrac{e^z-1-z}{z^4}$；(5) $\dfrac{\sin z}{z(e^{z^2}-1)}$；(6) $e^{1/z}$。

    **解。**

    (1) $\frac{1}{z(z-i)^2(z+i)^2}$：$z=0$ 为一阶极点，$z=\pm i$ 为二阶极点。

    (2) $\cos z$ 的零点 $z_k=k\pi+\frac\pi2$ 都是一阶零点（$(\cos z)'=-\sin z_k=\mp1\ne0$），所以 $z_k$ 都是 $\frac1{\cos z}$ 的一阶极点。

    (3) $z_k=k\pi$ 是 $\sin z$ 的一阶零点，因而是 $\sin^2z$ 的二阶零点，而 $\cos z_k=\pm1\ne0$，所以 $z_k$ 都是二阶极点。

    (4) $z=0$ 是 $z^4$ 的四阶零点。$e^z-1-z=\frac{z^2}{2!}+\frac{z^3}{3!}+\cdots$，$z=0$ 是分子的二阶零点，所以 $z=0$ 为 $4-2=2$ 阶极点。直接展开也可以：

    $$\frac{e^z-1-z}{z^4}=\frac{1}{2!\,z^2}+\frac{1}{3!\,z}+\frac{1}{4!}+\frac{z}{5!}+\cdots .$$

    (5) 奇点为 $z=0$ 及 $z^2=2k\pi i$（$k\ne0$）的根。在 $z=0$ 处：$z(e^{z^2}-1)=z\left(z^2+\frac{z^4}{2}+\cdots\right)$ 为三阶零点，$\sin z$ 为一阶零点，所以 $z=0$ 为二阶极点。其余各点处分母是一阶零点、分子不为零，都是一阶极点。

    (6) $e^{1/z}=\sum_{n=0}^\infty\frac{1}{n!\,z^n}$ 有无穷多负幂项，$z=0$ 为本性奇点。（也可以用极限：沿正实轴 $z\to0^+$ 时 $e^{1/z}\to+\infty$，沿负实轴 $z\to0^-$ 时 $e^{1/z}\to0$，极限既不存在也不是 $\infty$。）

!!! warning "易错点：不要只看分母"

    $\frac{\sin z}{z}$ 的分母在 $0$ 处有零点，但 $z=0$ 是可去奇点，不是极点。判断极点阶数时必须同时看**分子在该点的零点阶数**。例题 5.2(4) 若只看分母会误判为四阶极点。

### 5.1.4 无穷远点

!!! definition "定义 5.4（无穷远点作为孤立奇点）"

    若 $f(z)$ 在 $R<|z|<+\infty$ 内解析，则称 $z=\infty$ 为 $f$ 的**孤立奇点**。令 $t=\frac1z$，$g(t)=f\!\left(\frac1t\right)$，按 $t=0$ 是 $g$ 的可去奇点、$m$ 阶极点、本性奇点，相应地称 $z=\infty$ 是 $f$ 的**可去奇点、$m$ 阶极点、本性奇点**。

    等价地，看 $f$ 在 $R<|z|<+\infty$ 内的洛朗级数的**正幂**部分：没有正幂项、只有有限个正幂项（最高次为 $z^m$）、有无穷多个正幂项。

!!! example "例题 5.3"

    判断 $z=\infty$ 的类型：(1) $P(z)=a_nz^n+\cdots+a_0$（$a_n\ne0$）；(2) $e^z$；(3) $\dfrac{z}{z+1}$；(4) $\dfrac{1}{\sin z}$。

    **解。** (1) 正幂最高为 $z^n$：$n$ 阶极点。(2) $e^z=\sum\frac{z^n}{n!}$ 有无穷多正幂项：本性奇点。(3) $\frac{z}{z+1}=\frac{1}{1+\frac1z}=1-\frac1z+\frac1{z^2}-\cdots$（$|z|>1$），没有正幂项：可去奇点。(4) 奇点 $k\pi\to\infty$，$\frac1{\sin z}$ 在任何 $R<|z|<\infty$ 内都不解析，$\infty$ **不是孤立奇点**。

---

## 5.2 留数

### 5.2.1 留数的定义与留数定理

!!! definition "定义 5.5（留数）"

    设 $z_0$ 为 $f(z)$ 的孤立奇点，$f$ 在 $0<|z-z_0|<\delta$ 内的洛朗级数为 $\sum_{n=-\infty}^{+\infty}c_n(z-z_0)^n$。称 $c_{-1}$ 为 $f$ 在 $z_0$ 处的**留数**（residue），记作

    $$\operatorname{Res}[f(z),z_0]=c_{-1}=\frac{1}{2\pi i}\oint_Cf(z)\,\mathrm dz,$$

    其中 $C$ 为该去心邻域内绕 $z_0$ 的任意一条正向简单闭曲线。

名字的含义：把洛朗级数沿 $C$ 逐项积分，其余各项积分都是 $0$，只有 $c_{-1}$ 这一项“留下”了 $2\pi i\,c_{-1}$（例题 3.2）。

!!! theorem "定理 5.5（留数定理）"

    设 $f(z)$ 在区域 $D$ 内除有限个孤立奇点 $z_1,z_2,\dots,z_n$ 外处处解析，$C$ 是 $D$ 内包围诸奇点的一条正向简单闭曲线，则

    $$\oint_Cf(z)\,\mathrm dz=2\pi i\sum_{k=1}^n\operatorname{Res}[f(z),z_k].$$

??? proof "证明"

    以 $z_k$ 为心作互不相交、互不包含且都在 $C$ 内部的小圆 $C_k$。由复合闭路定理（定理 3.4），

    $$\oint_Cf(z)\,\mathrm dz=\sum_{k=1}^n\oint_{C_k}f(z)\,\mathrm dz=\sum_{k=1}^n2\pi i\operatorname{Res}[f(z),z_k].\qquad\blacksquare$$

<figure class="video" markdown>
<video controls preload="metadata" width="100%" poster="assets/posters/C5_ResidueTheorem.jpg">
<source src="assets/videos/C5_ResidueTheorem.mp4" type="video/mp4">
</video>
<figcaption>动画 5-2　留数定理 = 复合闭路定理 + 洛朗级数逐项积分。大回路换成绕各奇点的小圆，每个小圆只留下 $2\pi i\,c_{-1}$；回路外的奇点不参与。</figcaption>
</figure>

留数定理把闭路积分变成了**求留数**。下面给出各种奇点处留数的算法。

### 5.2.2 留数的计算规则

!!! theorem "定理 5.6（留数的计算规则）"

    1. **可去奇点**：$\operatorname{Res}[f,z_0]=0$。
    2. **本性奇点**：只能把 $f$ 在去心邻域内作洛朗展开，读出 $c_{-1}$。
    3. **$m$ 阶极点**（规则 III）：

        $$\operatorname{Res}[f,z_0]=\frac{1}{(m-1)!}\lim_{z\to z_0}\frac{\mathrm d^{m-1}}{\mathrm dz^{m-1}}\Big[(z-z_0)^mf(z)\Big].$$

    4. **一阶极点**（规则 I）：$\operatorname{Res}[f,z_0]=\lim\limits_{z\to z_0}(z-z_0)f(z)$。
    5. **$\frac PQ$ 型一阶极点**（规则 II）：若 $f=\frac{P}{Q}$，$P,Q$ 在 $z_0$ 处解析，$P(z_0)\ne0$，$Q(z_0)=0$，$Q'(z_0)\ne0$，则

        $$\operatorname{Res}[f,z_0]=\frac{P(z_0)}{Q'(z_0)}.$$

??? proof "规则 III 与规则 II 的推导"

    规则 III：$f=\dfrac{c_{-m}}{(z-z_0)^m}+\cdots+\dfrac{c_{-1}}{z-z_0}+c_0+\cdots$，两边乘以 $(z-z_0)^m$：

    $$(z-z_0)^mf(z)=c_{-m}+\cdots+c_{-1}(z-z_0)^{m-1}+c_0(z-z_0)^m+\cdots .$$

    右边是泰勒级数，$(z-z_0)^{m-1}$ 的系数等于 $\frac{1}{(m-1)!}\frac{\mathrm d^{m-1}}{\mathrm dz^{m-1}}[(z-z_0)^mf]\big|_{z_0}$，即为 $c_{-1}$。

    规则 II：$z_0$ 是一阶极点，由规则 I，

    $$\operatorname{Res}[f,z_0]=\lim_{z\to z_0}\frac{P(z)}{\frac{Q(z)-Q(z_0)}{z-z_0}}=\frac{P(z_0)}{Q'(z_0)}.\qquad\blacksquare$$

!!! tip "规则 III 中的 $m$ 取大了也没关系"

    若 $z_0$ 实际上是 $k$ 阶极点（$k<m$），用 $m$ 代入规则 III 仍得正确结果，因为 $(z-z_0)^mf(z)$ 仍是解析函数，它的 $(z-z_0)^{m-1}$ 项系数仍是 $c_{-1}$。但取**小**了就错了。拿不准阶数时宁可取大，只是求导会更麻烦。

!!! example "例题 5.4（一阶极点）"

    求 (1) $f(z)=\dfrac{1}{z(z-1)}$ 在 $0,1$ 处的留数；(2) $f(z)=\dfrac{1}{1+z^4}$ 在各奇点处的留数。

    **解。** (1) 规则 I：

    $$\operatorname{Res}[f,0]=\lim_{z\to0}\frac{1}{z-1}=-1,\qquad\operatorname{Res}[f,1]=\lim_{z\to1}\frac1z=1 .$$

    (2) 奇点 $z_k=e^{i\pi(2k+1)/4}$（$k=0,1,2,3$），都是一阶极点。用规则 I 需要计算 $\prod_{j\ne k}(z_k-z_j)$，很麻烦；用规则 II 则很简单：

    $$\operatorname{Res}[f,z_k]=\frac{1}{4z_k^3}=\frac{z_k}{4z_k^4}=-\frac{z_k}{4}\qquad(z_k^4=-1).$$

    即 $\operatorname{Res}[f,e^{i\pi/4}]=-\frac14e^{i\pi/4}$，$\operatorname{Res}[f,e^{3i\pi/4}]=-\frac14e^{3i\pi/4}$，等等。

!!! example "例题 5.5（高阶极点：规则 III 与洛朗展开的比较）"

    求 $f(z)=\dfrac{z-\sin z}{z^6}$ 在 $z=0$ 处的留数。

    **解法一（洛朗展开）。** $z-\sin z=\frac{z^3}{3!}-\frac{z^5}{5!}+\frac{z^7}{7!}-\cdots$，所以

    $$f(z)=\frac{1}{3!\,z^3}-\frac{1}{5!\,z}+\frac{z}{7!}-\cdots,\qquad\operatorname{Res}[f,0]=-\frac{1}{5!}=-\frac{1}{120}.$$

    **解法二（规则 III）。** 分子是三阶零点，$z=0$ 为三阶极点：

    $$\operatorname{Res}[f,0]=\frac{1}{2!}\lim_{z\to0}\left(\frac{z-\sin z}{z^3}\right)'' ,$$

    需要求二阶导数后再连续用洛必达法则，非常繁琐。若把 $m$ 取成 $6$，则

    $$\operatorname{Res}[f,0]=\frac{1}{5!}\lim_{z\to0}(z-\sin z)^{(5)}=\frac{1}{5!}\lim_{z\to0}(-\cos z)=-\frac{1}{120},$$

    反而简单。这里 $m=6$ 能用，正是上面“$m$ 取大了也没关系”的原因。

    **经验。** 分子、分母都有现成的泰勒展开时，优先用洛朗展开直接读 $c_{-1}$。

!!! example "例题 5.6（本性奇点）"

    求留数：(1) $z^2\cos\dfrac1z$ 在 $z=0$；(2) $(1+z)e^{1/z}$ 在 $z=0$；(3) $z^2\cos\dfrac{1}{z-1}$ 在 $z=1$；(4) $\dfrac{e^{1/z}}{z(z-1)}$ 在 $z=0,1$。

    **解。** 本性奇点只能用洛朗展开。

    (1) $z^2\cos\frac1z=z^2\left(1-\frac{1}{2!\,z^2}+\frac{1}{4!\,z^4}-\cdots\right)=z^2-\frac12+\frac{1}{4!\,z^2}-\cdots$，只有偶次幂，$\operatorname{Res}=0$。

    (2) $(1+z)\left(1+\frac1z+\frac{1}{2!\,z^2}+\cdots\right)$ 中 $z^{-1}$ 的系数为 $1+\frac1{2!}=\frac32$。

    (3) 令 $w=z-1$，$z^2=w^2+2w+1$，$\cos\frac1w=1-\frac{1}{2!\,w^2}+\frac{1}{4!\,w^4}-\cdots$。乘积中 $w^{-1}$ 的系数只来自 $2w\cdot\left(-\frac{1}{2!\,w^2}\right)$，等于 $-1$。所以 $\operatorname{Res}[f,1]=-1$。

    (4) $z=1$ 是一阶极点：$\operatorname{Res}[f,1]=\lim\limits_{z\to1}\dfrac{e^{1/z}}{z}=e$。

    $z=0$ 是本性奇点。在 $0<|z|<1$ 内，

    $$f(z)=-\frac1z\cdot\frac{1}{1-z}\cdot e^{1/z}=-\frac1z\left(\sum_{k=0}^\infty z^k\right)\left(\sum_{n=0}^\infty\frac{1}{n!\,z^n}\right).$$

    $z^{-1}$ 的系数等于两个括号乘积中常数项的相反数。常数项来自 $k=n$ 的各项，为 $\sum_{n\ge0}\frac{1}{n!}=e$，所以 $\operatorname{Res}[f,0]=-e$。

### 5.2.3 用留数定理计算闭路积分

!!! example "例题 5.7（教材例题）"

    计算：

    (1) $\displaystyle\oint_{|z|=2}\frac{\sin^2z}{z^2(z-1)}\,\mathrm dz$；(2) $\displaystyle\oint_{|z|=2}\frac{e^z}{z(z-1)^2}\,\mathrm dz$；(3) $\displaystyle\oint_{|z|=2}\sin\frac{z}{z-1}\,\mathrm dz$；(4) $\displaystyle\oint_{|z|=\frac12}\frac{\mathrm dz}{z^{101}(1-z^2)}$。

    **解。**

    (1) $z=0$：$\frac{\sin^2z}{z^2}\to1$，是可去奇点，留数为 $0$。$z=1$：一阶极点，留数 $\lim\limits_{z\to1}\frac{\sin^2z}{z^2}=\sin^21$。故 $I=2\pi i\sin^21$。

    (2) $z=0$ 为一阶极点，$\operatorname{Res}=\lim\limits_{z\to0}\frac{e^z}{(z-1)^2}=1$；$z=1$ 为二阶极点，

    $$\operatorname{Res}[f,1]=\lim_{z\to1}\left(\frac{e^z}{z}\right)'=\lim_{z\to1}\frac{e^z(z-1)}{z^2}=0 .$$

    故 $I=2\pi i(1+0)=2\pi i$。

    (3) $z=1$ 为本性奇点。$\sin\frac{z}{z-1}=\sin\left(1+\frac{1}{z-1}\right)=\sin1\cos\frac{1}{z-1}+\cos1\sin\frac1{z-1}$。$\cos\frac{1}{z-1}$ 只含偶次负幂，$\sin\frac{1}{z-1}=\frac{1}{z-1}-\cdots$。所以 $\operatorname{Res}=\cos1$，$I=2\pi i\cos1$。

    (4) $z=0$ 为 $101$ 阶极点，用规则 III 要求 $100$ 阶导数。改用洛朗展开：在 $0<|z|<1$ 内

    $$f(z)=\frac{1}{z^{101}}\sum_{n=0}^\infty z^{2n},$$

    $z^{-1}$ 项对应 $2n=100$，系数为 $1$，故 $I=2\pi i$。

!!! warning "易错点"

    - 只算**回路内部**奇点的留数；
    - 先判断奇点类型再选规则：例题 5.7(1) 中 $z=0$ 看起来像二阶极点，实际是可去奇点；
    - 规则 II 要求 $P(z_0)\ne0$、$Q'(z_0)\ne0$（一阶极点）；
    - 本性奇点不能用规则 I、II、III。

### 5.2.4 无穷远点的留数

!!! definition "定义 5.6（无穷远点的留数）"

    设 $f$ 在 $R<|z|<+\infty$ 内解析，$C$ 为该圆环内绕原点的任意一条正向简单闭曲线，定义

    $$\operatorname{Res}[f,\infty]=\frac{1}{2\pi i}\oint_{C^-}f(z)\,\mathrm dz=-c_{-1},$$

    其中 $c_{-1}$ 是 $f$ 在 $R<|z|<+\infty$ 内洛朗级数中 $z^{-1}$ 的系数。$C^-$ 为顺时针方向，即绕 $\infty$ 的“正向”（$\infty$ 在左手边）。

!!! theorem "定理 5.7（留数和定理）"

    若 $f$ 在扩充复平面上只有有限个孤立奇点 $z_1,\dots,z_n,\infty$，则

    $$\sum_{k=1}^n\operatorname{Res}[f,z_k]+\operatorname{Res}[f,\infty]=0 .$$

    计算公式：

    $$\operatorname{Res}[f,\infty]=-\operatorname{Res}\left[f\!\left(\frac1z\right)\frac{1}{z^2},\,0\right].$$

（取包围所有有限奇点的大圆 $C$，$\oint_C f\,\mathrm dz=2\pi i\sum_k\operatorname{Res}[f,z_k]=-2\pi i\operatorname{Res}[f,\infty]$。计算公式由 $z=\frac1t$ 换元得到。）

**用途。** 当回路内部的奇点很多或阶数很高、而外部只有少数简单奇点时，改为计算外部的留数：

$$\oint_Cf\,\mathrm dz=-2\pi i\Big(\sum_{\text{外部}}\operatorname{Res}+\operatorname{Res}[f,\infty]\Big).$$

**实用判别。** 若 $f=\frac{P}{Q}$ 是有理函数且 $\deg Q\ge\deg P+2$，则 $\operatorname{Res}[f,\infty]=0$。若 $\deg Q=\deg P+1$，则 $\operatorname{Res}[f,\infty]=-\frac{a}{b}$，其中 $a,b$ 分别为 $P,Q$ 的首项系数（因为 $f\sim\frac{a}{bz}$）。

!!! example "例题 5.8"

    计算 (1) $\displaystyle\oint_{|z|=2}\frac{z^{15}}{(z^2+1)^2(z^4+2)^3}\,\mathrm dz$；(2) $\displaystyle\oint_{|z|=2}\frac{\mathrm dz}{(z+i)^{10}(z-1)(z-3)}$。

    **解。** (1) 奇点 $\pm i$（二阶）与 $z^4=-2$ 的四个根（三阶，模为 $2^{1/4}<2$），全在圆内，直接算非常麻烦。改算 $\infty$ 处的留数：分母次数 $4+12=16$，分子 $15$，$f\sim\frac1z$，所以 $c_{-1}=1$，$\operatorname{Res}[f,\infty]=-1$。由定理 5.7，

    $$I=2\pi i\sum_{\text{圆内}}\operatorname{Res}=-2\pi i\operatorname{Res}[f,\infty]=2\pi i .$$

    (2) 圆内奇点 $-i$（$10$ 阶）与 $1$，圆外奇点 $3$。分母次数 $12$，$\operatorname{Res}[f,\infty]=0$；$\operatorname{Res}[f,3]=\frac{1}{(3+i)^{10}\cdot2}$。所以

    $$I=-2\pi i\left(\frac{1}{2(3+i)^{10}}+0\right)=-\frac{\pi i}{(3+i)^{10}} .$$

!!! warning "易错点：$\operatorname{Res}[f,\infty]$ 的符号"

    $\operatorname{Res}[f,\infty]=-c_{-1}$，前面有负号；而且即使 $\infty$ 是可去奇点，$\operatorname{Res}[f,\infty]$ 也**可以不为零**（例如 $f=\frac1z$，$\operatorname{Res}[f,\infty]=-1$）。这与有限点处“可去奇点留数为零”不同。

---

## 5.3 留数在定积分计算中的应用

基本思路：把实积分看成某个闭路积分的一部分（或者通过换元变成闭路积分），再用留数定理计算。

### 5.3.1 $\displaystyle\int_0^{2\pi}R(\cos\theta,\sin\theta)\,\mathrm d\theta$ 型

$R(u,v)$ 为 $u,v$ 的有理函数，且在单位圆上分母不为零。令 $z=e^{i\theta}$，$\theta$ 从 $0$ 到 $2\pi$ 时 $z$ 沿单位圆正向绕一周：

$$\mathrm d\theta=\frac{\mathrm dz}{iz},\qquad\cos\theta=\frac{z+z^{-1}}{2},\qquad\sin\theta=\frac{z-z^{-1}}{2i},$$

!!! theorem "公式 I"

    $$\int_0^{2\pi}R(\cos\theta,\sin\theta)\,\mathrm d\theta=\oint_{|z|=1}R\!\left(\frac{z^2+1}{2z},\frac{z^2-1}{2iz}\right)\frac{\mathrm dz}{iz}=2\pi i\sum_{|z_k|<1}\operatorname{Res}[F,z_k],$$

    $F$ 为右边闭路积分的被积函数，求和只取**单位圆内**的奇点。

<figure class="video" markdown>
<video controls preload="metadata" width="100%" poster="assets/posters/C5_UnitCircleSub.jpg">
<source src="assets/videos/C5_UnitCircleSub.mp4" type="video/mp4">
</video>
<figcaption>动画 5-3　令 $z=e^{i\theta}$，区间 $[0,2\pi]$ 卷成单位圆，三角有理式变成 $z$ 的有理式，只取单位圆内的极点。</figcaption>
</figure>

!!! example "例题 5.9"

    计算 $\displaystyle I=\int_0^{2\pi}\frac{\mathrm d\theta}{5+3\sin\theta}$。

    **解。** 代入：

    $$I=\oint_{|z|=1}\frac{1}{5+\frac{3(z^2-1)}{2iz}}\cdot\frac{\mathrm dz}{iz}=\oint_{|z|=1}\frac{2\,\mathrm dz}{3z^2+10iz-3}.$$

    $3z^2+10iz-3=0$ 的根为 $z=\frac{-10i\pm8i}{6}$，即 $-\frac i3$（圆内）与 $-3i$（圆外）。由规则 II，

    $$\operatorname{Res}\left[F,-\tfrac i3\right]=\frac{2}{6z+10i}\bigg|_{z=-i/3}=\frac{2}{8i}=\frac{1}{4i},\qquad I=2\pi i\cdot\frac1{4i}=\frac\pi2 .$$

!!! example "例题 5.10"

    计算 $\displaystyle I=\int_0^{2\pi}\frac{\cos2\theta}{5-4\cos\theta}\,\mathrm d\theta$。

    **解。** $\cos2\theta=\frac{z^2+z^{-2}}{2}$，$5-4\cos\theta=5-2(z+z^{-1})=-\frac{(2z-1)(z-2)}{z}$，

    $$I=\oint_{|z|=1}\frac{z^4+1}{2z^2}\cdot\frac{-z}{(2z-1)(z-2)}\cdot\frac{\mathrm dz}{iz}=\oint_{|z|=1}\frac{i(z^4+1)}{2z^2(2z-1)(z-2)}\,\mathrm dz .$$

    圆内奇点：$z=0$（二阶极点）与 $z=\frac12$（一阶极点）。

    $$\operatorname{Res}[F,0]=\frac i2\left(\frac{z^4+1}{2z^2-5z+2}\right)'\bigg|_{z=0}=\frac i2\cdot\frac{0\cdot2-1\cdot(-5)}{4}=\frac{5i}{8},$$

    $$\operatorname{Res}\left[F,\tfrac12\right]=\frac{i(z^4+1)}{2z^2\cdot2(z-2)}\bigg|_{z=1/2}=\frac{i\cdot\frac{17}{16}}{2\cdot\frac14\cdot2\cdot\left(-\frac32\right)}=-\frac{17i}{24}.$$

    $$I=2\pi i\left(\frac{15i}{24}-\frac{17i}{24}\right)=2\pi i\cdot\left(-\frac{i}{12}\right)=\frac{\pi}{6}.$$

    **注意** $\cos2\theta$ 带来了 $z^{-2}$，使 $z=0$ 也成为奇点。这类题最常见的错误就是漏掉 $z=0$ 处的留数。

### 5.3.2 $\displaystyle\int_{-\infty}^{+\infty}R(x)\,\mathrm dx$ 型

!!! theorem "引理 5.1（大圆弧引理）"

    设 $C_R$ 为上半圆周 $z=Re^{i\theta}$（$0\le\theta\le\pi$）。若在 $C_R$ 上 $zf(z)\to0$（$R\to\infty$，对 $\theta$ 一致），则

    $$\lim_{R\to\infty}\int_{C_R}f(z)\,\mathrm dz=0 .$$

（由长大不等式：$\left|\int_{C_R}f\,\mathrm dz\right|\le\max_{C_R}|f|\cdot\pi R=\pi\max_{C_R}|zf(z)|\to0$。）

!!! theorem "公式 II"

    设 $R(x)=\frac{P(x)}{Q(x)}$ 为有理函数，$\deg Q\ge\deg P+2$，且 $Q$ 在实轴上没有零点，则

    $$\int_{-\infty}^{+\infty}R(x)\,\mathrm dx=2\pi i\sum_{\operatorname{Im}z_k>0}\operatorname{Res}[R(z),z_k].$$

**推导。** 取闭路 $[-R,R]+C_R$，$R$ 足够大，使上半平面的奇点全在内部：$\int_{-R}^RR(x)\,\mathrm dx+\int_{C_R}R(z)\,\mathrm dz=2\pi i\sum_{\operatorname{Im}z_k>0}\operatorname{Res}$。由于 $\deg Q\ge\deg P+2$，$zR(z)\to0$，令 $R\to\infty$，由引理 5.1 即得。

<figure class="video" markdown>
<video controls preload="metadata" width="100%" poster="assets/posters/C5_SemicircleContour.jpg">
<source src="assets/videos/C5_SemicircleContour.mp4" type="video/mp4">
</video>
<figcaption>动画 5-4　补上半圆周：$\int_{-R}^{R}\frac{\mathrm dx}{1+x^2}$ 与大圆弧上的积分之和始终是 $2\pi i\operatorname{Res}[\frac{1}{1+z^2},i]=\pi$；$R$ 增大时大圆弧上的积分趋于 $0$。</figcaption>
</figure>

!!! example "例题 5.11"

    计算 (1) $\displaystyle\int_{-\infty}^{+\infty}\frac{\mathrm dx}{1+x^4}$；(2) $\displaystyle\int_{-\infty}^{+\infty}\frac{x^2}{(x^2+1)(x^2+4)}\,\mathrm dx$；(3) $\displaystyle\int_{-\infty}^{+\infty}\frac{\mathrm dx}{(x^2+1)^2}$。

    **解。** 三个积分都满足 $\deg Q\ge\deg P+2$，实轴上无奇点。

    (1) 上半平面的奇点为 $e^{i\pi/4},e^{3i\pi/4}$。由例题 5.4，留数为 $-\frac{z_k}{4}$：

    $$I=2\pi i\cdot\left(-\frac14\right)\left(e^{i\pi/4}+e^{3i\pi/4}\right)=-\frac{\pi i}{2}\cdot i\sqrt2=\frac{\pi}{\sqrt2}.$$

    (2) 上半平面的奇点为 $i,2i$（一阶极点）：

    $$\operatorname{Res}[R,i]=\frac{i^2}{2i\,(i^2+4)}=\frac{-1}{6i},\qquad\operatorname{Res}[R,2i]=\frac{(2i)^2}{\big((2i)^2+1\big)\cdot4i}=\frac{1}{3i},$$

    $$I=2\pi i\left(-\frac{1}{6i}+\frac{1}{3i}\right)=\frac\pi3 .$$

    (3) $z=i$ 为二阶极点：

    $$\operatorname{Res}[R,i]=\lim_{z\to i}\left(\frac{1}{(z+i)^2}\right)'=\lim_{z\to i}\frac{-2}{(z+i)^3}=\frac{-2}{(2i)^3}=\frac{1}{4i},\qquad I=\frac\pi2 .$$

### 5.3.3 $\displaystyle\int_{-\infty}^{+\infty}R(x)e^{iax}\,\mathrm dx$（$a>0$）型

!!! theorem "引理 5.2（约当引理）"

    设 $a>0$，在上半圆周 $C_R$ 上 $g(z)\to0$（$R\to\infty$，对 $\theta$ 一致），则

    $$\lim_{R\to\infty}\int_{C_R}g(z)e^{iaz}\,\mathrm dz=0 .$$

（关键在于上半平面内 $|e^{iaz}|=e^{-ay}\le1$，且在 $\theta$ 接近 $\frac\pi2$ 时指数衰减。证明利用 $\sin\theta\ge\frac{2\theta}{\pi}$（$0\le\theta\le\frac\pi2$），这里从略。）

!!! theorem "公式 III"

    设 $R(x)=\frac{P(x)}{Q(x)}$，$\deg Q\ge\deg P+1$，$Q$ 在实轴上没有零点，$a>0$，则

    $$\int_{-\infty}^{+\infty}R(x)e^{iax}\,\mathrm dx=2\pi i\sum_{\operatorname{Im}z_k>0}\operatorname{Res}\big[R(z)e^{iaz},z_k\big].$$

    分别取实部、虚部，就得到 $\displaystyle\int_{-\infty}^{+\infty}R(x)\cos ax\,\mathrm dx$ 与 $\displaystyle\int_{-\infty}^{+\infty}R(x)\sin ax\,\mathrm dx$（$R$ 为实有理函数时）。

!!! example "例题 5.12"

    计算 (1) $\displaystyle\int_0^{+\infty}\frac{\cos x}{x^2+1}\,\mathrm dx$；(2) $\displaystyle\int_{-\infty}^{+\infty}\frac{x\sin x}{x^2+1}\,\mathrm dx$。

    **解。** (1) 被积函数是偶函数，原积分 $=\frac12\int_{-\infty}^{+\infty}\frac{\cos x}{x^2+1}\,\mathrm dx=\frac12\operatorname{Re}\int_{-\infty}^{+\infty}\frac{e^{ix}}{x^2+1}\,\mathrm dx$。上半平面只有 $z=i$：

    $$\int_{-\infty}^{+\infty}\frac{e^{ix}}{x^2+1}\,\mathrm dx=2\pi i\cdot\frac{e^{i\cdot i}}{2i}=\frac{\pi}{e},\qquad\text{原积分}=\frac{\pi}{2e}.$$

    (2) $\deg Q-\deg P=1$，公式 II 不能用，但公式 III 可以：

    $$\int_{-\infty}^{+\infty}\frac{xe^{ix}}{x^2+1}\,\mathrm dx=2\pi i\cdot\frac{i\,e^{-1}}{2i}=\frac{\pi i}{e},$$

    取虚部得 $\displaystyle\int_{-\infty}^{+\infty}\frac{x\sin x}{x^2+1}\,\mathrm dx=\frac{\pi}{e}$。

!!! warning "易错点：为什么不能直接用 $\cos az$"

    $\cos z=\frac{e^{iz}+e^{-iz}}{2}$ 中的 $e^{-iz}$ 在上半平面 $|e^{-iz}|=e^{y}\to\infty$，大圆弧上的积分不趋于零。所以必须先把 $\cos ax$ 换成 $e^{iax}$，算完再取实部。$a<0$ 时应改用下半平面（或对 $a$ 取共轭）。

### 5.3.4 被积函数在实轴上有一阶极点

!!! theorem "引理 5.3（小半圆引理）"

    设 $z_0$ 为 $f$ 的一阶极点，$C_r$ 为上半圆周 $z=z_0+re^{i\theta}$，$\theta$ 从 $\pi$ 变到 $0$（从 $z_0-r$ 经上方到 $z_0+r$，顺时针方向），则

    $$\lim_{r\to0}\int_{C_r}f(z)\,\mathrm dz=-\pi i\operatorname{Res}[f,z_0].$$

（$f=\frac{c_{-1}}{z-z_0}+g(z)$，$g$ 在 $z_0$ 解析。第一项在顺时针半圆上的积分为 $c_{-1}\int_\pi^0i\,\mathrm d\theta=-\pi ic_{-1}$，第二项 $\to0$。）

!!! example "例题 5.13（教材例 5.27）"

    计算狄利克雷积分 $\displaystyle\int_0^{+\infty}\frac{\sin x}{x}\,\mathrm dx$。

    **解。** 考虑 $f(z)=\frac{e^{iz}}{z}$，$z=0$ 在实轴上，是一阶极点。取闭路：$[-R,-r]$ $+$ 小半圆 $C_r$（顺时针，从 $-r$ 经上方到 $r$）$+$ $[r,R]$ $+$ 大半圆 $C_R$。闭路内没有奇点，积分为 $0$：

    $$\int_{-R}^{-r}\frac{e^{ix}}{x}\,\mathrm dx+\int_{r}^{R}\frac{e^{ix}}{x}\,\mathrm dx+\int_{C_r}f\,\mathrm dz+\int_{C_R}f\,\mathrm dz=0 .$$

    在前两个积分中，把第一个积分的 $x$ 换成 $-x$，合并为 $\int_r^R\frac{e^{ix}-e^{-ix}}{x}\,\mathrm dx=2i\int_r^R\frac{\sin x}{x}\,\mathrm dx$。

    令 $R\to\infty$：由约当引理，$\int_{C_R}\to0$。令 $r\to0$：由引理 5.3，$\int_{C_r}\to-\pi i\operatorname{Res}[f,0]=-\pi i$。于是

    $$2i\int_0^{+\infty}\frac{\sin x}{x}\,\mathrm dx-\pi i=0\quad\Longrightarrow\quad\int_0^{+\infty}\frac{\sin x}{x}\,\mathrm dx=\frac\pi2 .$$

!!! success "三类积分对照"

    | 类型 | 条件 | 闭路 | 结果 |
    |---|---|---|---|
    | $\int_0^{2\pi}R(\cos\theta,\sin\theta)\,\mathrm d\theta$ | 单位圆上分母不为零 | 单位圆 $z=e^{i\theta}$ | $2\pi i\sum_{\lvert z_k\rvert<1}\operatorname{Res}$ |
    | $\int_{-\infty}^{\infty}R(x)\,\mathrm dx$ | $\deg Q\ge\deg P+2$，实轴无奇点 | 线段 + 上半大圆 | $2\pi i\sum_{\operatorname{Im}z_k>0}\operatorname{Res}$ |
    | $\int_{-\infty}^{\infty}R(x)e^{iax}\,\mathrm dx$ | $a>0$，$\deg Q\ge\deg P+1$ | 线段 + 上半大圆（约当引理） | $2\pi i\sum_{\operatorname{Im}z_k>0}\operatorname{Res}[Re^{iaz}]$ |
    | 实轴上有一阶极点 | 同上 | 再挖去小半圆 | 实轴上的极点贡献 $\pi i\operatorname{Res}$ |

---

## 知识点归纳总结

!!! abstract "知识点归纳总结"

    ```mermaid
    flowchart TD
        LAU["第4章：洛朗级数<br/>∮ f dz = 2πi c₋₁"] --> CLS["孤立奇点分类<br/>负幂项：无 / 有限 / 无穷"]
        CLS --> RM["可去奇点<br/>lim f 有限"]
        CLS --> PO["m 阶极点<br/>lim f = ∞"]
        CLS --> ES["本性奇点<br/>无极限"]
        ZERO["零点阶数<br/>z₀ 处前 m−1 阶导数为 0"] --> PO
        LAU --> RES["留数 Res = c₋₁"]
        RES --> RULE["计算规则<br/>I: lim (z−z₀)f · II: P/Q' · III: m 阶公式<br/>本性奇点：洛朗展开"]
        RES --> INF["Res[f,∞] = −c₋₁<br/>所有留数之和 = 0"]
        RULE --> THM["留数定理<br/>∮ f dz = 2πi Σ Res"]
        INF --> THM
        THM --> T1["∫(0→2π) R(cos, sin) dθ<br/>z = exp(iθ)"]
        THM --> T2["∫ R(x) dx<br/>上半大圆"]
        THM --> T3["∫ R(x) exp(iax) dx<br/>约当引理"]
        THM -.->|"后续课程"| LAP["拉普拉斯 / Z 反变换"]
    ```

    **核心结论**

    | 内容 | 结论 |
    |---|---|
    | $m$ 阶零点 | $f(z_0)=\cdots=f^{(m-1)}(z_0)=0,\ f^{(m)}(z_0)\ne0$ |
    | 极点阶数 | $z_0$ 是 $Q$ 的 $m$ 阶、$P$ 的 $n$ 阶零点，$m>n$ ⇒ $m-n$ 阶极点 |
    | 留数定理 | $\oint_Cf\,\mathrm dz=2\pi i\sum\operatorname{Res}[f,z_k]$ |
    | 一阶极点 | $\operatorname{Res}=\lim(z-z_0)f=\frac{P(z_0)}{Q'(z_0)}$ |
    | $m$ 阶极点 | $\operatorname{Res}=\frac{1}{(m-1)!}\lim\frac{\mathrm d^{m-1}}{\mathrm dz^{m-1}}[(z-z_0)^mf]$ |
    | 无穷远点 | $\operatorname{Res}[f,\infty]=-c_{-1}=-\operatorname{Res}[f(\frac1z)\frac1{z^2},0]$ |
    | 留数和 | 有限奇点留数之和 $+\operatorname{Res}[f,\infty]=0$ |
    | 实积分 | 见上面“三类积分对照”表 |

    **方法清单**

    1. 判断奇点类型：有现成展开式 → 看洛朗级数；分式 → 比较分子分母零点阶数；其余 → 求极限；
    2. 求留数：一阶极点用规则 I / II；低阶极点用规则 III；高阶极点或本性奇点用洛朗展开；内部奇点太多时改算外部和 $\infty$ 处的留数；
    3. 闭路积分：列出内部奇点 → 求留数 → 乘 $2\pi i$；
    4. 三角有理式积分：$z=e^{i\theta}$，化简后只取单位圆内的奇点，特别注意 $z=0$；
    5. 无穷积分：先验证次数条件；含 $\cos ax,\sin ax$ 时换成 $e^{iax}$，最后取实部或虚部；偶函数先化成全实轴积分。

    **易错点汇总**

    - 可去奇点误判为极点（分子也有零点）；
    - 规则 III 中 $m$ 取小了；
    - 本性奇点用了极点公式；
    - 算了回路外部的奇点，或者漏了内部的奇点（尤其是 $z=0$）；
    - $\operatorname{Res}[f,\infty]=-c_{-1}$ 漏掉负号；
    - 无穷积分中取了下半平面的奇点，或直接对 $\cos az$ 用约当引理；
    - 忘了验证 $\deg Q\ge\deg P+2$（或 $+1$）以及实轴上有没有奇点。

---

## 练习题

!!! example "练习 1"

    指出下列函数的所有孤立奇点（含 $\infty$）及其类型：(1) $\dfrac{1}{z^3-z^5}$；(2) $\dfrac{z}{\sin z}$；(3) $e^{\frac{1}{z-1}}$。

??? success "解答"

    (1) $\frac{1}{z^3(1-z)(1+z)}$：$z=0$ 为三阶极点，$z=\pm1$ 为一阶极点；$z=\infty$：$f\sim-\frac{1}{z^5}\to0$，是可去奇点（并且是 $f$ 的五阶零点）。

    (2) $z=0$：$\frac{z}{\sin z}\to1$，可去奇点。$z=k\pi$（$k\ne0$）：分母一阶零点、分子不为零，一阶极点。$z=\infty$：极点 $k\pi\to\infty$，不是孤立奇点。

    (3) $z=1$：$e^{\frac1{z-1}}=\sum\frac{1}{n!(z-1)^n}$，本性奇点。$z=\infty$：$\lim\limits_{z\to\infty}e^{\frac1{z-1}}=1$，可去奇点。

!!! example "练习 2"

    $z=0$ 是 $f(z)=z^2\left(e^{z^2}-1\right)$ 的几阶零点？是 $\dfrac{1}{f(z)}$ 的什么奇点？

??? success "解答"

    $e^{z^2}-1=z^2+\frac{z^4}{2}+\cdots$，所以 $f=z^4+\frac{z^6}{2}+\cdots$，$z=0$ 为四阶零点。由定理 5.4，它是 $\frac1f$ 的四阶极点。

!!! example "练习 3"

    求留数：(1) $\dfrac{z}{(z-1)(z+1)^2}$ 在 $\pm1$；(2) $\dfrac{1-\cos z}{z^5}$ 在 $0$；(3) $\tan z$ 在 $\dfrac\pi2$。

??? success "解答"

    (1) $z=1$ 一阶：$\lim\limits_{z\to1}\frac{z}{(z+1)^2}=\frac14$。$z=-1$ 二阶：$\lim\limits_{z\to-1}\left(\frac{z}{z-1}\right)'=\lim\limits_{z\to-1}\frac{-1}{(z-1)^2}=-\frac14$。

    (2) $1-\cos z=\frac{z^2}{2!}-\frac{z^4}{4!}+\frac{z^6}{6!}-\cdots$，除以 $z^5$ 后 $z^{-1}$ 的系数为 $-\frac{1}{4!}=-\frac1{24}$。

    (3) $\tan z=\frac{\sin z}{\cos z}$，规则 II：$\frac{\sin(\pi/2)}{-\sin(\pi/2)}=-1$。

!!! example "练习 4"

    计算 $\displaystyle\oint_{|z|=2}\frac{z}{z^4-1}\,\mathrm dz$。

??? success "解答"

    四个一阶极点 $\pm1,\pm i$ 都在圆内。规则 II：$\operatorname{Res}=\frac{z}{4z^3}=\frac{1}{4z^2}$，在 $\pm1$ 处为 $\frac14$，在 $\pm i$ 处为 $-\frac14$，总和为 $0$，所以积分为 $0$。

    或者：$f\sim\frac1{z^3}$，$\operatorname{Res}[f,\infty]=0$，由留数和定理直接得 $0$。

!!! example "练习 5"

    计算 $\displaystyle\oint_{|z|=1}z^2e^{1/z}\,\mathrm dz$。

??? success "解答"

    $z^2e^{1/z}=z^2\sum_{n\ge0}\frac{1}{n!\,z^n}$，$z^{-1}$ 项对应 $n=3$，系数为 $\frac{1}{3!}=\frac16$。积分为 $2\pi i\cdot\frac16=\frac{\pi i}{3}$。

!!! example "练习 6"

    计算 $\displaystyle\oint_{|z|=3}\frac{\mathrm dz}{(z-1)^5(z-4)}$。

??? success "解答"

    圆内只有 $z=1$，但它是五阶极点，用规则 III 要求四阶导数。改用留数和定理：$\operatorname{Res}[f,\infty]=0$（分母次数比分子高 $6$），$\operatorname{Res}[f,4]=\frac{1}{3^5}=\frac1{243}$，所以 $\operatorname{Res}[f,1]=-\frac{1}{243}$，积分为 $-\dfrac{2\pi i}{243}$。

!!! example "练习 7"

    计算 $\displaystyle\int_0^{2\pi}\frac{\mathrm d\theta}{2+\cos\theta}$。

??? success "解答"

    $z=e^{i\theta}$：$2+\frac{z+z^{-1}}{2}=\frac{z^2+4z+1}{2z}$，

    $$I=\oint_{|z|=1}\frac{2z}{z^2+4z+1}\cdot\frac{\mathrm dz}{iz}=\oint_{|z|=1}\frac{2}{i(z^2+4z+1)}\,\mathrm dz .$$

    根为 $-2\pm\sqrt3$，圆内为 $-2+\sqrt3$。$\operatorname{Res}=\frac{2}{i(2z+4)}\Big|_{-2+\sqrt3}=\frac{1}{i\sqrt3}$，$I=2\pi i\cdot\frac{1}{i\sqrt3}=\dfrac{2\pi}{\sqrt3}$。

!!! example "练习 8"

    计算 $\displaystyle\int_{-\infty}^{+\infty}\frac{\cos x}{x^2+4x+5}\,\mathrm dx$。

??? success "解答"

    $x^2+4x+5=(x+2)^2+1$，零点 $-2\pm i$，上半平面为 $-2+i$。次数差为 $2\ge1$，用公式 III：

    $$\int_{-\infty}^{+\infty}\frac{e^{ix}}{x^2+4x+5}\,\mathrm dx=2\pi i\cdot\frac{e^{i(-2+i)}}{2i}=\pi e^{-1}e^{-2i}=\frac{\pi}{e}(\cos2-i\sin2).$$

    取实部：$\displaystyle\int_{-\infty}^{+\infty}\frac{\cos x}{x^2+4x+5}\,\mathrm dx=\frac{\pi\cos2}{e}$。（顺便得到 $\int_{-\infty}^{+\infty}\frac{\sin x}{x^2+4x+5}\,\mathrm dx=-\frac{\pi\sin2}{e}$。）

---

**下一章预告。** 前五章都在算：导数、积分、级数、留数。第 6 章回到第 1 章的几何观点：解析函数在 $f'(z_0)\ne0$ 处局部是“旋转 + 伸缩”，所以它**保持角度**。利用这个性质，可以把复杂区域变成单位圆或上半平面，这在二维电磁场、流体问题中很有用。→ [第 6 章 共形映射](ch6.md)
