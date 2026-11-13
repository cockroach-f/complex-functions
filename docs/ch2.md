# 第 2 章 解析函数

!!! tip "本章目标"

    1. 掌握复导数的定义，理解它比实二元函数的可微**苛刻得多**（$\Delta z\to0$ 的方向任意）；
    2. 熟练运用**柯西-黎曼（C-R）方程**判断函数在一点可导、在区域内解析，并求导数；
    3. 理解解析函数与调和函数的关系，会**已知实部（或虚部）构造解析函数**——期末必考题型；
    4. 掌握初等函数 $e^z,\ \operatorname{Ln}z,\ z^\alpha,\ \sin z,\ \cos z$ 的定义、多值性、主值与解析区域，知道它们与实函数的**不同之处**。

**本章主线。** 第 1 章 1.5.3 说过：复极限要求 $z$ 从**任何方向**趋近都得到同一个值。把这个要求用在差商 $\frac{\Delta w}{\Delta z}$ 上，就得到复导数；沿实轴、沿虚轴两个方向算出来必须一致，立刻逼出 C-R 方程；C-R 方程又推出实部虚部都是调和函数。整章就是这一个条件的连锁反应。

---

## 2.1 解析函数的概念

### 2.1.1 导数与微分

!!! definition "定义 2.1（导数）"

    设函数 $w=f(z)$ 在点 $z_0$ 的某邻域 $U(z_0)$ 内有定义，$z_0+\Delta z\in U(z_0)$。若极限

    $$\lim_{\Delta z\to0}\frac{\Delta w}{\Delta z}=\lim_{\Delta z\to0}\frac{f(z_0+\Delta z)-f(z_0)}{\Delta z}$$

    存在且有限，则称 $f(z)$ 在 $z_0$ **可导**，该极限称为 $f$ 在 $z_0$ 的**导数**，记作 $f'(z_0)$ 或 $\left.\dfrac{\mathrm dw}{\mathrm dz}\right|_{z=z_0}$。若 $f$ 在区域 $D$ 内处处可导，则称 $f$ 在 $D$ 内可导，得到导函数 $f'(z)$。

形式上和一元实函数完全相同，区别全在“$\Delta z\to0$ 的方式任意”这一点上。

!!! example "例题 2.1"

    (1) 求 $f(z)=z^2$ 与 $g(z)=\dfrac1z$ 的导数；(2) 讨论 $h(z)=\bar z$ 的可导性。

    **解。** (1)
    $$\lim_{\Delta z\to0}\frac{(z+\Delta z)^2-z^2}{\Delta z}=\lim_{\Delta z\to0}(2z+\Delta z)=2z,\qquad
    \lim_{\Delta z\to0}\frac{\frac{1}{z+\Delta z}-\frac1z}{\Delta z}=\lim_{\Delta z\to0}\frac{-1}{z(z+\Delta z)}=-\frac{1}{z^2}\ (z\ne0).$$

    同理 $(z^n)'=nz^{n-1}$（$n$ 为正整数），$(C)'=0$。

    (2) $\dfrac{\Delta h}{\Delta z}=\dfrac{\overline{z+\Delta z}-\bar z}{\Delta z}=\dfrac{\overline{\Delta z}}{\Delta z}$。令 $\Delta z=re^{i\theta}$，则它等于 $e^{-2i\theta}$，随方向变化，极限不存在（正是第 1 章例题 1.12 的 $\bar z/z$）。所以 $\bar z$ **处处不可导**——尽管它的实部 $x$、虚部 $-y$ 都是最光滑的函数。

!!! definition "定义（微分）"

    若存在复常数 $A$，使 $\Delta w=f(z_0+\Delta z)-f(z_0)=A\,\Delta z+o(|\Delta z|)$，则称 $f$ 在 $z_0$ **可微**，$A\,\Delta z$ 称为微分，记作 $\mathrm dw=A\,\mathrm dz$（取 $f(z)=z$ 可知 $\mathrm dz=\Delta z$）。

!!! theorem "定理（可导、可微、连续的关系）"

    $f$ 在 $z_0$ 可导 $\iff$ $f$ 在 $z_0$ 可微，且 $A=f'(z_0)$，$\mathrm dw=f'(z_0)\,\mathrm dz$；可导 $\Rightarrow$ 连续，反之不然。

    这与一元实函数一致。对比：二元实函数“偏导数存在”既推不出可微也推不出连续，而复变函数“可导”包含了“各个方向的变化率都一致”这一强得多的信息。

**求导法则**与实函数相同：

$$(f\pm g)'=f'\pm g',\qquad (fg)'=f'g+fg',\qquad \left(\frac fg\right)'=\frac{f'g-fg'}{g^2}\ (g\ne0),$$

$$[f(g(z))]'=f'(g(z))\,g'(z),\qquad \varphi'(w)=\frac{1}{f'(z)}\ \ (z=\varphi(w)\ \text{为}\ w=f(z)\ \text{的单值反函数},\ f'(z)\neq0).$$

### 2.1.2 解析函数

!!! definition "定义 2.2（解析）"

    (1) 若 $f(z)$ 在 $z_0$ **及 $z_0$ 的某个邻域内**处处可导，则称 $f$ 在 $z_0$ **解析**；

    (2) 若 $f$ 在区域 $D$ 内每一点解析，则称 $f$ 在 $D$ 内解析，或称 $f$ 是 $D$ 内的**解析函数**（全纯函数）；

    (3) 若 $f$ 在 $z_0$ 不解析，则称 $z_0$ 为 $f$ 的**奇点**。

!!! warning "易错点：“点可导”与“点解析”不同"

    - 在一点：**解析 ⇒ 可导，可导 ⇏ 解析**（解析要求一整个邻域都可导）；
    - 在区域上：**可导 ⇔ 解析**（区域是开集，每点都有邻域在区域内）。

    典型反例见下面例题 2.3：$|z|^2$ 只在 $z=0$ 一点可导，因而处处不解析。

**性质。** 在 $D$ 内解析的函数之和、差、积、商（分母不为零处）仍解析；解析函数的复合仍解析。因此**多项式在全平面解析，有理函数 $P(z)/Q(z)$ 在 $Q(z)\ne0$ 处解析**，其奇点就是分母的零点。

!!! example "例题 2.2"

    讨论 $f(z)=\dfrac{z+3}{4z^2+1}$ 的解析性并求导数。

    **解。** 分子分母都是多项式，处处解析；分母 $4z^2+1=0\iff z=\pm\frac{i}{2}$。所以 $f$ 在 $\mathbb C\setminus\{\pm\frac i2\}$ 内解析，$\pm\frac i2$ 是奇点，且

    $$f'(z)=\frac{(4z^2+1)-8z(z+3)}{(4z^2+1)^2}=\frac{-4z^2-24z+1}{(4z^2+1)^2}.$$

!!! example "例题 2.3"

    讨论 $f(z)=|z|^2$ 的可导性与解析性。

    **解。** $|z|^2=z\bar z$，

    $$\frac{\Delta f}{\Delta z}=\frac{(z+\Delta z)(\bar z+\overline{\Delta z})-z\bar z}{\Delta z}=\bar z+\overline{\Delta z}+z\,\frac{\overline{\Delta z}}{\Delta z}.$$

    $z=0$ 时它等于 $\overline{\Delta z}\to0$，所以 $f'(0)=0$；$z\ne0$ 时最后一项 $z\,\frac{\overline{\Delta z}}{\Delta z}$ 的极限不存在。故 $|z|^2$ **仅在 $z=0$ 可导，处处不解析**。

### 2.1.3 柯西-黎曼方程

逐个用定义判断太麻烦。下面的定理把“复可导”翻译成关于 $u,v$ 的条件。

!!! theorem "定理 2.1（一点可导的充要条件）"

    设 $f(z)=u(x,y)+iv(x,y)$ 在 $z=x+iy$ 的邻域内有定义。则 $f$ 在 $z$ 可导的**充要条件**是：

    1. $u(x,y),\ v(x,y)$ 在点 $(x,y)$ **可微**；
    2. 在该点满足**柯西-黎曼方程（C-R 方程）**
    $$\boxed{\ \frac{\partial u}{\partial x}=\frac{\partial v}{\partial y},\qquad \frac{\partial u}{\partial y}=-\frac{\partial v}{\partial x}\ }$$

    此时
    $$f'(z)=\frac{\partial u}{\partial x}+i\frac{\partial v}{\partial x}=\frac{\partial v}{\partial y}-i\frac{\partial u}{\partial y}=\frac{\partial u}{\partial x}-i\frac{\partial u}{\partial y}=\frac{\partial v}{\partial y}+i\frac{\partial v}{\partial x}.$$

??? proof "证明"

    **必要性。** 设 $f'(z)=a+ib$。可导即可微：$\Delta w=(a+ib)(\Delta x+i\Delta y)+o(|\Delta z|)$。分开实部虚部：

    $$\Delta u=a\Delta x-b\Delta y+o(|\Delta z|),\qquad \Delta v=b\Delta x+a\Delta y+o(|\Delta z|).$$

    故 $u,v$ 在 $(x,y)$ 可微，且 $u_x=a=v_y,\ u_y=-b=-v_x$。

    **充分性。** 由 $u,v$ 可微，$\Delta u=u_x\Delta x+u_y\Delta y+o(|\Delta z|)$，$\Delta v=v_x\Delta x+v_y\Delta y+o(|\Delta z|)$。代入 C-R 方程（$u_y=-v_x,\ v_y=u_x$）：

    $$\Delta w=\Delta u+i\Delta v=(u_x+iv_x)\Delta x+(-v_x+iu_x)\Delta y+o(|\Delta z|)=(u_x+iv_x)(\Delta x+i\Delta y)+o(|\Delta z|),$$

    于是 $\dfrac{\Delta w}{\Delta z}\to u_x+iv_x$。$\blacksquare$

**几何意义。** 在 $z_0$ 附近把 $f$ 线性化，$(\Delta x,\Delta y)\mapsto(\Delta u,\Delta v)$ 由雅可比矩阵 $\begin{pmatrix}u_x&u_y\\v_x&v_y\end{pmatrix}$ 给出。复可导要求这个线性映射就是“乘以复数 $f'(z_0)=a+ib$”，即矩阵形如 $\begin{pmatrix}a&-b\\b&a\end{pmatrix}$ = **旋转 + 伸缩**。C-R 方程正是“雅可比矩阵具有这种形状”。

<figure class="video" markdown>
<video controls preload="metadata" width="100%" poster="assets/posters/C2_CauchyRiemannGeometry.jpg">
<source src="assets/videos/C2_CauchyRiemannGeometry.mp4" type="video/mp4">
</video>
<figcaption>动画 2-1　C-R 方程的几何意义：$z^2$ 把小正方形映成（近似）小正方形；$\bar z$ 把它翻转，$x+2iy$ 把它压成长方形——后两者不满足 C-R，不可导。</figcaption>
</figure>

!!! theorem "定理 2.2（区域内解析的充要条件）"

    $f=u+iv$ 在区域 $D$ 内解析 $\iff$ $u,v$ 在 $D$ 内可微，且满足 C-R 方程。

    **推论（实用判别法）。** 若 $u_x,u_y,v_x,v_y$ 在 $D$ 内**存在且连续**，并满足 C-R 方程，则 $f$ 在 $D$ 内解析（偏导数连续 ⇒ 可微）。

!!! warning "易错点：C-R 方程只是必要条件的一半"

    只验证 C-R 方程**不够**，还需要可微。反例：$f(z)=\sqrt{|xy|}$ 在原点有 $u_x=u_y=v_x=v_y=0$，满足 C-R，但沿 $y=x$ 趋于 0 时 $\dfrac{\Delta f}{\Delta z}=\dfrac{|x|}{(1+i)x}$ 随 $x$ 的符号变化，不可导。做题时一般用推论：先说明偏导数连续，再验证 C-R。

!!! success "判别解析性的标准流程"

    1. 写出 $u(x,y),\ v(x,y)$；
    2. 求四个偏导数，说明它们连续（初等函数在定义域内自动连续）；
    3. 解 C-R 方程组，得到满足方程的点集 $E$；
    4. **结论**：$f$ 在 $E$ 上可导；若 $E$ 包含某个区域 $D$，则 $f$ 在 $D$ 内解析；若 $E$ 只是孤立点或曲线（不含任何开集），则**处处不解析**；
    5. 需要时用 $f'=u_x+iv_x$ 求导数。

!!! example "例题 2.4"

    讨论下列函数的可导性与解析性：(1) $f(z)=x^2+iy^2$；(2) $f(z)=e^x(\cos y+i\sin y)$。

    **解。** (1) $u=x^2,\ v=y^2$，$u_x=2x,\ u_y=0,\ v_x=0,\ v_y=2y$，均连续。C-R：$2x=2y$，$0=0$。所以 $f$ **仅在直线 $y=x$ 上可导**，直线不含任何开集，故**处处不解析**；在 $y=x$ 上 $f'(z)=u_x+iv_x=2x$。

    (2) $u=e^x\cos y,\ v=e^x\sin y$：

    $$u_x=e^x\cos y=v_y,\qquad u_y=-e^x\sin y=-v_x,$$

    四个偏导数连续且处处满足 C-R，所以 $f$ 在全平面解析，且 $f'(z)=u_x+iv_x=e^x(\cos y+i\sin y)=f(z)$。这个函数就是 2.3 节要定义的 $e^z$：**$(e^z)'=e^z$**。

!!! example "例题 2.5"

    设 $f(z)=x^2+axy+by^2+i(cx^2+dxy+y^2)$，问常数 $a,b,c,d$ 取何值时 $f$ 在复平面内解析？

    **解。** $u_x=2x+ay,\ u_y=ax+2by,\ v_x=2cx+dy,\ v_y=dx+2y$。C-R 方程对所有 $x,y$ 成立：

    $$2x+ay=dx+2y\ \Rightarrow\ d=2,\ a=2;\qquad ax+2by=-(2cx+dy)\ \Rightarrow\ a=-2c,\ 2b=-d.$$

    解得 $a=2,\ b=-1,\ c=-1,\ d=2$。此时 $f=(x^2+2xy-y^2)+i(-x^2+2xy+y^2)=(1-i)z^2$。

??? info "补充：极坐标形式的 C-R 方程（选读）"

    令 $z=re^{i\theta}$，$f=u(r,\theta)+iv(r,\theta)$，则 C-R 方程等价于

    $$\frac{\partial u}{\partial r}=\frac1r\frac{\partial v}{\partial\theta},\qquad \frac{\partial v}{\partial r}=-\frac1r\frac{\partial u}{\partial\theta},\qquad f'(z)=e^{-i\theta}\left(\frac{\partial u}{\partial r}+i\frac{\partial v}{\partial r}\right).$$

    对 $\ln z=\ln r+i\theta$ 用起来最方便：$u_r=\frac1r=\frac1r v_\theta$，$v_r=0=-\frac1r u_\theta$，$f'=e^{-i\theta}\frac1r=\frac1z$。

---

## 2.2 解析函数与调和函数

### 2.2.1 调和函数

!!! definition "定义 2.3（调和函数）"

    若二元实函数 $\varphi(x,y)$ 在区域 $D$ 内有连续的二阶偏导数，且满足**拉普拉斯方程**

    $$\Delta\varphi=\frac{\partial^2\varphi}{\partial x^2}+\frac{\partial^2\varphi}{\partial y^2}=0,$$

    则称 $\varphi$ 为 $D$ 内的**调和函数**。

!!! theorem "定理 2.3"

    若 $f(z)=u+iv$ 在区域 $D$ 内解析，则 $u,v$ 都是 $D$ 内的调和函数。

??? proof "证明"

    第 3 章将证明：解析函数的导数仍解析，从而 $u,v$ 具有任意阶连续偏导数（这里先承认）。由 C-R 方程 $u_x=v_y,\ u_y=-v_x$，分别对 $x$、$y$ 求偏导：

    $$u_{xx}=v_{yx},\qquad u_{yy}=-v_{xy}.$$

    由二阶混合偏导连续，$v_{xy}=v_{yx}$，相加得 $u_{xx}+u_{yy}=0$。同理 $v_{xx}+v_{yy}=0$。$\blacksquare$

### 2.2.2 共轭调和函数

!!! definition "定义 2.4（共轭调和函数）"

    设 $u,v$ 都是 $D$ 内的调和函数，且满足 C-R 方程 $u_x=v_y,\ u_y=-v_x$，则称 $v$ 是 $u$ 的**共轭调和函数**。

!!! theorem "定理 2.4"

    $f=u+iv$ 在区域 $D$ 内解析 $\iff$ 在 $D$ 内 $v$ 是 $u$ 的共轭调和函数。

!!! warning "易错点：共轭关系不对称"

    $v$ 是 $u$ 的共轭调和函数，**不能**推出 $u$ 是 $v$ 的共轭调和函数——实际上 $-u$ 才是 $v$ 的共轭调和函数（因为 $-if=v-iu$ 解析）。例如 $z^2$：$v=2xy$ 是 $u=x^2-y^2$ 的共轭调和函数，而 $v$ 的共轭调和函数是 $-(x^2-y^2)$。

**几何意义：等值线正交。** 若 $f=u+iv$ 解析且 $f'(z)\ne0$，则两族曲线 $u(x,y)=c_1$ 与 $v(x,y)=c_2$ **处处正交**：

$$\nabla u\cdot\nabla v=u_xv_x+u_yv_y\overset{\text{C-R}}{=}u_xv_x+(-v_x)u_x=0 .$$

在二维静电场中，若 $u$ 是电位，则 $u=c_1$ 是等势线，$v=c_2$ 就是电力线；$f=u+iv$ 称为**复势**。这是电磁场课程里用复变函数解平面场问题的出发点。

<figure class="video" markdown>
<video controls preload="metadata" width="100%" poster="assets/posters/C2_HarmonicOrthogonal.jpg">
<source src="assets/videos/C2_HarmonicOrthogonal.mp4" type="video/mp4">
</video>
<figcaption>动画 2-2　$z^2$ 与 $1/z$ 的实部、虚部等值线：两族曲线处处正交（等势线 ⟂ 电力线）。</figcaption>
</figure>

### 2.2.3 构造解析函数（必考）

**问题。** 已知实部 $u$（或虚部 $v$），求解析函数 $f=u+iv$，并满足给定条件（如 $f(z_0)=w_0$）。

**依据。** (1) $u$（或 $v$）本身必须是调和函数——**第一步必须验证**；(2) $u,v$ 满足 C-R 方程。

!!! success "三种方法（以已知 $u$ 为例）"

    **① 偏积分法。** 由 $v_y=u_x$ 对 $y$ 积分：$v=\int u_x\,\mathrm dy+\varphi(x)$；再代入 $v_x=-u_y$ 定出 $\varphi'(x)$，积分得 $\varphi(x)$。

    **② 全微分（曲线积分）法。** $\mathrm dv=v_x\,\mathrm dx+v_y\,\mathrm dy=-u_y\,\mathrm dx+u_x\,\mathrm dy$，它是全微分（因为 $u$ 调和），故
    $$v(x,y)=\int_{(x_0,y_0)}^{(x,y)}-u_y\,\mathrm dx+u_x\,\mathrm dy+C,$$
    通常沿折线 $(x_0,y_0)\to(x,y_0)\to(x,y)$ 积分。

    **③ 不定积分法（最快）。** $f'(z)=u_x-iu_y$。把右边写成 $z$ 的函数：**令 $y=0,\ x=z$**（解析函数由它在实轴上的值唯一确定），得 $f'(z)=U(z)$，积分得 $f(z)=\int U(z)\,\mathrm dz+C$。已知 $v$ 时用 $f'(z)=v_y+iv_x$。

    注意常数：已知 $u$ 时 $C$ 只能是**纯虚数** $ic$；已知 $v$ 时 $C$ 是实数。

!!! example "例题 2.6（教材例 2.6 改编）"

    已知 $u(x,y)=x^3-3xy^2$，求解析函数 $f(z)=u+iv$，使 $f(i)=-i$。

    **解。** (1) **验证调和**：$u_{xx}=6x,\ u_{yy}=-6x$，$\Delta u=0$。✓

    (2) **偏积分法求 $v$**：$v_y=u_x=3x^2-3y^2\Rightarrow v=3x^2y-y^3+\varphi(x)$。再由 $v_x=-u_y$：$6xy+\varphi'(x)=6xy\Rightarrow\varphi'(x)=0,\ \varphi=c$。所以 $v=3x^2y-y^3+c$。

    (3) **定常数**：$f=(x^3-3xy^2)+i(3x^2y-y^3+c)$，代入 $x=0,y=1$：$f(i)=i(-1+c)=-i\Rightarrow c=0$。

    所以 $f(z)=(x^3-3xy^2)+i(3x^2y-y^3)=z^3$。

    **用方法③验证**：$f'(z)=u_x-iu_y=(3x^2-3y^2)+6xyi$，令 $y=0,x=z$：$f'(z)=3z^2$，$f=z^3+ic$，由 $f(i)=-i+ic=-i$ 得 $c=0$。✓

!!! example "例题 2.7"

    已知 $v(x,y)=2xy+3x$，求解析函数 $f=u+iv$，使 $f(0)=0$。

    **解。** $\Delta v=0+0=0$，调和。✓ 用不定积分法：$f'(z)=v_y+iv_x=2x+i(2y+3)$，令 $y=0,x=z$：$f'(z)=2z+3i$，故 $f(z)=z^2+3iz+C$（$C$ 为实数）。由 $f(0)=0$ 得 $C=0$：

    $$f(z)=z^2+3iz,\qquad u=x^2-y^2-3y .$$

---

## 2.3 初等函数

复初等函数是实初等函数的推广：定义方式尽量保持一致，自变量取实数时与实函数重合。讨论的重点是：**定义、单值/多值、解析区域、与实函数的区别**。核心是指数函数——其余初等函数都由它定义。

### 2.3.1 指数函数

!!! definition "定义 2.5（指数函数）"

    对 $z=x+iy$，定义
    $$e^z=\exp z=e^x(\cos y+i\sin y).$$

由例题 2.4(2) 与定义直接验证，$e^z$ 有下列性质：

!!! success "指数函数的性质"

    1. **单值**，在全平面解析，$(e^z)'=e^z$；$y=0$ 时就是实指数 $e^x$；$x=0$ 时就是欧拉公式 $e^{iy}=\cos y+i\sin y$；
    2. $|e^z|=e^x$，$\operatorname{Arg}e^z=y+2k\pi$：**实部决定模，虚部决定辐角**；
    3. $e^z\ne0$（因为 $e^x>0$）；
    4. 加法定理 $e^{z_1+z_2}=e^{z_1}e^{z_2}$；
    5. **周期性**：$e^{z+2k\pi i}=e^z$，基本周期 $2\pi i$——这是实指数没有的；
    6. $\lim\limits_{z\to\infty}e^z$ **不存在**（沿正实轴 $\to\infty$，沿负实轴 $\to0$）。

!!! warning "易错点"

    $e^z$ 只是 $\exp z$ 的记号，不能理解为“$e$ 的 $z$ 次方”再按乘幂的定义去算（那样会得到多值，见 2.3.3）。$e^{z}$ 永远单值。

**映射性质。** $w=e^z$ 把竖线 $x=c$ 映成圆周 $|w|=e^c$，把横线 $y=c$ 映成射线 $\arg w=c$；因而把**带形** $0<\operatorname{Im}z<\pi$ 映成**上半平面**，把带形 $0<\operatorname{Im}z<2\pi$ 映成**去掉正实轴的平面**（第 6 章常用）。

<figure class="video" markdown>
<video controls preload="metadata" width="100%" poster="assets/posters/C2_ExpMap.jpg">
<source src="assets/videos/C2_ExpMap.mp4" type="video/mp4">
</video>
<figcaption>动画 2-3　$w=e^z$：竖线变圆周、横线变射线；带形 $0<\operatorname{Im}z<\pi$ 变上半平面；周期 $2\pi i$ 导致上下相邻的带形重复覆盖 $w$ 平面。</figcaption>
</figure>

### 2.3.2 对数函数

!!! definition "定义 2.6（对数函数）"

    满足方程 $e^w=z\ (z\ne0)$ 的函数 $w=f(z)$ 称为**对数函数**，记作 $w=\operatorname{Ln}z$。

**计算。** 令 $z=re^{i\theta}$，$w=u+iv$，由 $e^ue^{iv}=re^{i\theta}$ 得 $e^u=r$，$v=\theta+2k\pi$，即

$$\boxed{\operatorname{Ln}z=\ln|z|+i\operatorname{Arg}z=\ln|z|+i\arg z+2k\pi i,\quad k\in\mathbb Z.}$$

- 对数函数是**无穷多值**函数，各值的虚部相差 $2\pi$ 的整数倍；
- **主值**：取主辐角得到 $\ln z=\ln|z|+i\arg z$，于是 $\operatorname{Ln}z=\ln z+2k\pi i$；当 $z=x>0$ 时 $\ln z$ 就是实对数；
- 对固定的 $k$，$\ln z+2k\pi i$ 称为 $\operatorname{Ln}z$ 的一个**分支**。

!!! example "例题 2.8"

    求 $\operatorname{Ln}(-1)$、$\operatorname{Ln}(-i)$、$\operatorname{Ln}(-1+i)$ 及其主值。

    **解。**

    - $\operatorname{Ln}(-1)=\ln1+i\pi+2k\pi i=(2k+1)\pi i$，主值 $\ln(-1)=\pi i$。（**负数在复数范围内有对数**。）
    - $\operatorname{Ln}(-i)=\ln1+i\left(-\frac\pi2\right)+2k\pi i=\left(2k-\frac12\right)\pi i$，主值 $-\frac\pi2 i$。
    - $\operatorname{Ln}(-1+i)=\ln\sqrt2+i\left(\frac{3\pi}{4}+2k\pi\right)$，主值 $\frac12\ln2+\frac{3\pi}{4}i$。

!!! success "对数函数的性质"

    1. 定义域为 $z\ne0$；
    2. **运算法则（集合意义下）**：$\operatorname{Ln}(z_1z_2)=\operatorname{Ln}z_1+\operatorname{Ln}z_2$，$\operatorname{Ln}\dfrac{z_1}{z_2}=\operatorname{Ln}z_1-\operatorname{Ln}z_2$；但 $\operatorname{Ln}z^n\ne n\operatorname{Ln}z$，$\operatorname{Ln}\sqrt[n]z=\frac1n\operatorname{Ln}z$；
    3. **解析性**：主值 $\ln z$（以及每个分支）在**除去原点和负实轴**的区域 $D_0=\mathbb C\setminus(-\infty,0]$ 内解析，且
    $$(\ln z)'=\frac1z,\qquad (\operatorname{Ln}z)'=\frac1z\ (\text{各分支});$$
    在负实轴上，$\arg z$ 从 $\pi$ 跳到 $-\pi$（第 1 章例题 1.13），$\ln z$ 不连续，更不解析。负实轴称为 $\ln z$ 的**支割线**。

??? proof "$(\ln z)'=1/z$ 的证明"

    在 $D_0$ 内 $w=\ln z$ 是 $z=e^w$ 的单值连续反函数，$\frac{\mathrm dz}{\mathrm dw}=e^w\ne0$，由反函数求导法则 $\frac{\mathrm dw}{\mathrm dz}=\frac1{e^w}=\frac1z$。（也可以用 2.1.3 选读中的极坐标 C-R 方程。）

<figure class="video" markdown>
<video controls preload="metadata" width="100%" poster="assets/posters/C2_LogHelix.jpg">
<source src="assets/videos/C2_LogHelix.mp4" type="video/mp4">
</video>
<figcaption>动画 2-4　把 $\operatorname{Arg}z$ 画成高度，得到“螺旋楼梯”：$z$ 绕原点一圈，$\operatorname{Ln}z$ 的虚部升高 $2\pi$。沿负实轴剪开、只取 $-\pi<\arg z\le\pi$ 这一层，就得到单值解析的主值 $\ln z$。</figcaption>
</figure>

!!! example "例题 2.9"

    求 (1) $f(z)=\ln(1+z)$；(2) $g(z)=\ln z^2$ 的导数，并指出解析区域。

    **解。** (1) $\ln\zeta$ 在 $\zeta\notin(-\infty,0]$ 处解析，$\zeta=1+z\in(-\infty,0]\iff z\in(-\infty,-1]$。所以 $f$ 在 $\mathbb C\setminus(-\infty,-1]$ 内解析，$f'(z)=\dfrac{1}{1+z}$。

    (2) $z^2\in(-\infty,0]\iff z$ 在虚轴上（含原点）。所以 $g$ 在**除去虚轴**的区域内解析，$g'(z)=\dfrac{2z}{z^2}=\dfrac2z$。

### 2.3.3 乘幂与幂函数

!!! definition "定义 2.7（乘幂与幂函数）"

    设 $a\ne0$，$b$ 为复数，定义**乘幂** $a^b=e^{b\operatorname{Ln}a}$。当 $a=z$ 为变量、$\alpha$ 为复常数时，$w=z^\alpha=e^{\alpha\operatorname{Ln}z}\ (z\ne0)$ 称为**幂函数**；$\alpha>0$ 时补充定义 $0^\alpha=0$。

由于 $\operatorname{Ln}z$ 多值，$z^\alpha=e^{\alpha\ln z}\cdot e^{2k\pi i\alpha}$，取值个数取决于 $\alpha$：

| $\alpha$ | $z^\alpha$ 的取值 | 解析性 |
|---|---|---|
| 正整数 $n$ | $z^n$，**单值** | 全平面解析，$(z^n)'=nz^{n-1}$ |
| 负整数 $-n$ | $1/z^n$，单值 | $z\ne0$ 处解析 |
| $0$ | $z^0=1$ | — |
| 有理数 $\frac mn$（既约，$n>1$） | **$n$ 个值**，即 $\sqrt[n]{z^m}$ | 各分支在 $\mathbb C\setminus(-\infty,0]$ 内解析 |
| 无理数或虚数 | **无穷多值** | 各分支在 $\mathbb C\setminus(-\infty,0]$ 内解析 |

各分支的导数 $(z^\alpha)'=\alpha z^{\alpha-1}$（取同一分支）。

!!! example "例题 2.10"

    求 $i^i$ 与 $1^{\sqrt2}$ 的值。

    **解。** $i^i=e^{i\operatorname{Ln}i}=e^{i\cdot i\left(\frac\pi2+2k\pi\right)}=e^{-\left(\frac\pi2+2k\pi\right)}$，$k\in\mathbb Z$——**全是正实数**，主值 $e^{-\pi/2}\approx0.2079$。

    $1^{\sqrt2}=e^{\sqrt2\operatorname{Ln}1}=e^{2\sqrt2k\pi i}=\cos(2\sqrt2k\pi)+i\sin(2\sqrt2k\pi)$，无穷多个值。**不要想当然地认为 $1^{\sqrt2}=1$**（$1$ 只是 $k=0$ 时的主值）。

### 2.3.4 三角函数与双曲函数

由欧拉公式，对实数 $y$ 有 $\cos y=\frac{e^{iy}+e^{-iy}}{2},\ \sin y=\frac{e^{iy}-e^{-iy}}{2i}$。把 $y$ 换成复数 $z$ 作为定义：

!!! definition "定义 2.8（三角函数）与定义 2.10（双曲函数）"

    $$\cos z=\frac{e^{iz}+e^{-iz}}{2},\qquad \sin z=\frac{e^{iz}-e^{-iz}}{2i},\qquad \tan z=\frac{\sin z}{\cos z},\ \ \cot z=\frac{\cos z}{\sin z};$$

    $$\operatorname{ch}z=\frac{e^z+e^{-z}}{2},\qquad \operatorname{sh}z=\frac{e^z-e^{-z}}{2},\qquad \operatorname{th}z=\frac{\operatorname{sh}z}{\operatorname{ch}z}.$$

!!! success "三角函数的性质"

    1. $\sin z,\cos z$ 单值、全平面解析，$(\sin z)'=\cos z$，$(\cos z)'=-\sin z$；
    2. 周期 $2\pi$；$\sin z$ 奇、$\cos z$ 偶；**欧拉公式对复数成立**：$e^{iz}=\cos z+i\sin z$；
    3. 实函数的三角恒等式全部照搬：$\sin^2z+\cos^2z=1$、和角公式等；
    4. **零点**与实函数相同：$\sin z=0\iff z=k\pi$，$\cos z=0\iff z=k\pi+\frac\pi2$；
    5. 与双曲函数的关系：$\cos(iz)=\operatorname{ch}z$，$\sin(iz)=i\operatorname{sh}z$；于是
    $$\sin(x+iy)=\sin x\operatorname{ch}y+i\cos x\operatorname{sh}y,\qquad \cos(x+iy)=\cos x\operatorname{ch}y-i\sin x\operatorname{sh}y .$$

!!! warning "易错点：$|\sin z|\le1$ 不再成立"

    $\sin(iy)=i\operatorname{sh}y$，$|\sin(iy)|=|\operatorname{sh}y|\to\infty\ (y\to\infty)$。复正弦、余弦在复平面上**无界**。（第 3 章的刘维尔定理会告诉你：在全平面解析的非常数函数**必然**无界。）

!!! example "例题 2.11"

    求 $\cos i$ 与 $\sin(1+2i)$。

    **解。** $\cos i=\dfrac{e^{i\cdot i}+e^{-i\cdot i}}{2}=\dfrac{e^{-1}+e}{2}=\operatorname{ch}1\approx1.543>1$。

    $\sin(1+2i)=\sin1\operatorname{ch}2+i\cos1\operatorname{sh}2=\dfrac{e^2+e^{-2}}{2}\sin1+i\,\dfrac{e^2-e^{-2}}{2}\cos1$。

??? proof "零点只有实数 $k\pi$ 的证明"

    $\sin z=0\iff e^{iz}=e^{-iz}\iff e^{2iz}=1\iff 2iz=2k\pi i\iff z=k\pi$。（用到 $e^w=1\iff w=2k\pi i$。）

### 2.3.5 反三角函数与反双曲函数

!!! definition "定义 2.9（反三角函数）"

    若 $\cos w=z$，则称 $w$ 为 $z$ 的**反余弦函数**，记作 $w=\operatorname{Arccos}z$。类似定义 $\operatorname{Arcsin}z,\ \operatorname{Arctan}z$。

**计算**（反三角函数都是对数函数的复合）：由 $z=\frac{e^{iw}+e^{-iw}}{2}$ 得 $(e^{iw})^2-2ze^{iw}+1=0$，$e^{iw}=z+\sqrt{z^2-1}$（根式取双值），所以

$$\operatorname{Arccos}z=-i\operatorname{Ln}\left(z+\sqrt{z^2-1}\right),\qquad \operatorname{Arcsin}z=-i\operatorname{Ln}\left(iz+\sqrt{1-z^2}\right),\qquad \operatorname{Arctan}z=-\frac i2\operatorname{Ln}\frac{1+iz}{1-iz}.$$

反双曲函数：$\operatorname{Arsh}z=\operatorname{Ln}\left(z+\sqrt{z^2+1}\right)$，$\operatorname{Arch}z=\operatorname{Ln}\left(z+\sqrt{z^2-1}\right)$，$\operatorname{Arth}z=\frac12\operatorname{Ln}\frac{1+z}{1-z}$。

!!! example "例题 2.12"

    解方程 $\sin z=2$。

    **解。** $z=\operatorname{Arcsin}2=-i\operatorname{Ln}\left(2i+\sqrt{1-4}\right)=-i\operatorname{Ln}\left(2i\pm\sqrt3 i\right)=-i\operatorname{Ln}\left[(2\pm\sqrt3)i\right]$。

    $(2\pm\sqrt3)i$ 的模为 $2\pm\sqrt3$、辐角为 $\frac\pi2+2k\pi$，故

    $$z=-i\left[\ln(2\pm\sqrt3)+i\left(\tfrac\pi2+2k\pi\right)\right]=\frac\pi2+2k\pi-i\ln(2\pm\sqrt3),\quad k\in\mathbb Z .$$

    （注意 $\ln(2-\sqrt3)=-\ln(2+\sqrt3)$，所以也可写成 $\frac\pi2+2k\pi\pm i\ln(2+\sqrt3)$。）实数范围无解的方程，在复数范围有无穷多个解。

### 2.3.6 实初等函数 vs 复初等函数

| | 实函数 | 复函数 |
|---|---|---|
| $e^x$ / $e^z$ | 单值，恒正，无周期 | 单值，$\ne0$ 但可取负值、虚值，周期 $2\pi i$ |
| $\ln x$ / $\operatorname{Ln}z$ | 只对 $x>0$ 有定义，单值 | 对 $z\ne0$ 有定义，**无穷多值**；主值在 $\mathbb C\setminus(-\infty,0]$ 解析 |
| $x^a$ / $z^\alpha$ | 通常单值 | 一般多值，单值性取决于 $\alpha$ |
| $\sin x,\cos x$ / $\sin z,\cos z$ | 有界 $\le1$ | **无界** |
| $\arcsin x$ / $\operatorname{Arcsin}z$ | 只对 $\lvert x\rvert\le1$ 有定义 | 对任意 $z$ 有定义，多值 |

---

## 知识点归纳总结

!!! abstract "知识点归纳总结"

    ```mermaid
    flowchart TD
        L["第1章：极限与方向无关"] --> D["复导数 f'(z)<br/>可导 ⇔ 可微 ⇒ 连续"]
        D --> A["解析 = 邻域内处处可导<br/>奇点 = 不解析的点"]
        D --> CR["C-R 方程<br/>u_x = v_y, u_y = −v_x<br/>（+ 可微）"]
        CR --> J["几何：局部 = 旋转 + 伸缩"]
        CR --> H["u, v 调和<br/>等值线正交"]
        H --> K["已知 u 求 f<br/>偏积分 / 全微分 / 不定积分"]
        A --> E["exp(z)：单值，周期 2πi"]
        E --> LN["Ln z = ln|z| + i Arg z<br/>无穷多值，支割线"]
        LN --> P["z 的 α 次 = exp(α Ln z)"]
        E --> TR["sin z, cos z：无界"]
        LN --> AT["Arcsin, Arccos, Arctan"]
        J -.->|"第6章"| CF["共形映射"]
        A -.->|"第3章"| CI["柯西积分定理"]
    ```

    **核心结论**

    | 内容 | 结论 |
    |---|---|
    | 一点可导 | $u,v$ 可微 + C-R；$f'=u_x+iv_x=v_y-iu_y$ |
    | 区域解析 | 偏导数连续 + C-R（实用判别法） |
    | 调和 | 解析函数的 $u,v$ 满足 $\Delta u=\Delta v=0$ |
    | 构造 $f$ | $f'(z)=u_x-iu_y$ 或 $v_y+iv_x$，令 $y=0,x=z$ |
    | 指数 | $e^z=e^x(\cos y+i\sin y)$，$(e^z)'=e^z$，周期 $2\pi i$ |
    | 对数 | $\operatorname{Ln}z=\ln\lvert z\rvert+i\arg z+2k\pi i$，$(\ln z)'=\frac1z$ |
    | 乘幂 | $a^b=e^{b\operatorname{Ln}a}$ |
    | 三角 | $\sin z=\frac{e^{iz}-e^{-iz}}{2i}$，$\cos z=\frac{e^{iz}+e^{-iz}}{2}$ |
    | 反三角 | $\operatorname{Arcsin}z=-i\operatorname{Ln}(iz+\sqrt{1-z^2})$ 等 |

    **方法清单**

    1. 判断可导 / 解析：写 $u,v$ → 偏导连续 → 解 C-R → 点集是否含开集；
    2. 证明不可导：沿两个方向（实轴、虚轴）求差商极限，或直接说明 C-R 不成立；
    3. 构造解析函数：**先验证调和** → 不定积分法最快 → 用条件定常数（注意常数是实数还是纯虚数）；
    4. 计算 $\operatorname{Ln},\ a^b,\ \operatorname{Arcsin}$ 等：一律化为 $\operatorname{Ln}$，求模与辐角，写全 $2k\pi i$；
    5. 求复合函数的解析区域：找出使内层函数落在支割线上的点，挖掉。

    **易错点汇总**

    - 点可导 ≠ 点解析；
    - 只验证 C-R 不够，需要可微（偏导连续）；
    - $v$ 是 $u$ 的共轭调和函数 ⇏ $u$ 是 $v$ 的共轭调和函数；
    - 构造解析函数忘记验证调和、常数类型写错；
    - $\operatorname{Ln}z$、$a^b$ 漏写 $2k\pi i$；$\operatorname{Ln}z^n\ne n\operatorname{Ln}z$；
    - 认为 $|\sin z|\le1$、$1^\alpha=1$、负数没有对数——在复数范围都不对。

---

## 练习题

!!! example "练习 1"

    讨论 $f(z)=xy^2+ix^2y$ 的可导性与解析性，在可导点求导数。

??? success "解答"

    $u=xy^2,\ v=x^2y$：$u_x=y^2,\ u_y=2xy,\ v_x=2xy,\ v_y=x^2$，均连续。C-R：$y^2=x^2$ 且 $2xy=-2xy$，即 $xy=0$ 且 $|x|=|y|$，只有 $x=y=0$。

    所以 $f$ 仅在 $z=0$ 可导，$f'(0)=u_x+iv_x=0$；处处不解析。

!!! example "练习 2"

    设 $f(z)=my^3+nx^2y+i(x^3+lxy^2)$ 在复平面内解析，求 $l,m,n$，并把 $f$ 写成 $z$ 的函数。

??? success "解答"

    $u_x=2nxy,\ v_y=2lxy\Rightarrow n=l$；$u_y=3my^2+nx^2,\ -v_x=-3x^2-ly^2\Rightarrow n=-3,\ 3m=-l$。

    所以 $l=n=-3,\ m=1$，$f=(y^3-3x^2y)+i(x^3-3xy^2)=iz^3$。

    （验证：$iz^3=i[(x^3-3xy^2)+i(3x^2y-y^3)]=(y^3-3x^2y)+i(x^3-3xy^2)$ ✓。）

!!! example "练习 3"

    已知 $u=x^2-y^2+xy$，求解析函数 $f(z)=u+iv$，使 $f(i)=-1+i$。

??? success "解答"

    $\Delta u=2-2=0$，调和。$f'(z)=u_x-iu_y=(2x+y)-i(-2y+x)$，令 $y=0,x=z$：$f'(z)=(2-i)z$，故 $f(z)=\frac{2-i}{2}z^2+ic$。

    $f(i)=-\frac{2-i}{2}+ic=-1+\frac i2+ic=-1+i\Rightarrow c=\frac12$。所以 $f(z)=\left(1-\frac i2\right)z^2+\frac i2$。

!!! example "练习 4"

    求 (1) $\operatorname{Ln}(-3+4i)$；(2) $(1+i)^i$；(3) 方程 $e^z=1+\sqrt3 i$ 的全部解。

??? success "解答"

    (1) $|-3+4i|=5$，第二象限，$\arg=\pi-\arctan\frac43$：$\operatorname{Ln}(-3+4i)=\ln5+i\left(\pi-\arctan\frac43+2k\pi\right)$。

    (2) $(1+i)^i=e^{i\operatorname{Ln}(1+i)}=e^{i\left[\ln\sqrt2+i\left(\frac\pi4+2k\pi\right)\right]}=e^{-\left(\frac\pi4+2k\pi\right)}\left[\cos\left(\tfrac12\ln2\right)+i\sin\left(\tfrac12\ln2\right)\right]$。

    (3) $z=\operatorname{Ln}(1+\sqrt3 i)=\ln2+i\left(\frac\pi3+2k\pi\right)$，$k\in\mathbb Z$。

!!! example "练习 5"

    证明 $|\sin z|^2=\sin^2x+\operatorname{sh}^2y$，并由此说明 $\sin z$ 在复平面上无界。

??? success "解答"

    $\sin z=\sin x\operatorname{ch}y+i\cos x\operatorname{sh}y$，

    $$|\sin z|^2=\sin^2x\operatorname{ch}^2y+\cos^2x\operatorname{sh}^2y=\sin^2x(1+\operatorname{sh}^2y)+\cos^2x\operatorname{sh}^2y=\sin^2x+\operatorname{sh}^2y .$$

    $|\sin z|\ge|\operatorname{sh}y|\to\infty\ (y\to\infty)$，故无界。

!!! example "练习 6"

    $\operatorname{Ln}z^2=2\operatorname{Ln}z$ 成立吗？说明理由。

??? success "解答"

    不成立。$\operatorname{Ln}z^2=2\ln|z|+i(2\arg z+2k\pi)$，$2\operatorname{Ln}z=2\ln|z|+i(2\arg z+4k\pi)$。作为集合，后者只是前者的一半（缺少 $k$ 为奇数对应的那些值）。例如 $z=i$：$\operatorname{Ln}(-1)$ 包含 $\pi i$，而 $2\operatorname{Ln}i=2\cdot\left(\frac\pi2+2k\pi\right)i=(\pi+4k\pi)i$ 不包含 $3\pi i$。

!!! example "练习 7"

    指出 $f(z)=\ln(z^2+1)$（主值）的解析区域并求导数。

??? success "解答"

    需 $z^2+1\notin(-\infty,0]$。$z^2+1=t\le0\iff z^2=t-1\le-1\iff z=\pm i s,\ s\ge1$。所以 $f$ 在除去虚轴上两条射线 $\{iy:y\ge1\}$ 与 $\{iy:y\le-1\}$ 的区域内解析，$f'(z)=\dfrac{2z}{z^2+1}$。

!!! example "练习 8"

    证明：若 $f$ 在区域 $D$ 内解析且 $f'(z)\equiv0$，则 $f$ 在 $D$ 内为常数；若 $|f|$ 在 $D$ 内为常数，则 $f$ 也为常数。

??? success "解答"

    (1) $f'=u_x+iv_x=0\Rightarrow u_x=v_x=0$，由 C-R 又有 $u_y=-v_x=0,\ v_y=u_x=0$。$u,v$ 的梯度在连通区域内恒为零，故 $u,v$ 为常数。

    (2) 设 $u^2+v^2=c$。若 $c=0$ 则 $f\equiv0$。若 $c\ne0$，对 $x,y$ 求偏导：$uu_x+vv_x=0,\ uu_y+vv_y=0$。用 C-R 把 $u_y,v_y$ 换成 $-v_x,u_x$：$uu_x+vv_x=0,\ -uv_x+vu_x=0$。这是关于 $(u_x,v_x)$ 的齐次线性方程组，系数行列式 $u^2+v^2=c\ne0$，故 $u_x=v_x=0$，即 $f'\equiv0$，由 (1) $f$ 为常数。

---

**下一章预告。** 第 3 章把积分搬到复平面：沿曲线积分 $\int_C f(z)\,\mathrm dz$。你会发现，解析函数的闭路积分**只取决于回路里面有没有奇点**——这就是柯西积分定理。→ [第 3 章 复变函数的积分](ch3.md)
