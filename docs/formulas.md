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
