# Minimum Snap QP Solver

## 轨迹描述
由于每一段轨迹都是一个多项式，因此可以将轨迹表示为，
$$
f(t)=\left\{\begin{array}{cc}
f_1(t) \doteq \sum_{i=0}^N p_{1, i} t^i & T_0 \leq t \leq T_1 \\
f_2(t) \doteq \sum_{i=0}^N p_{2, i} t^i & T_1 \leq t \leq T_2 \\
\vdots & \vdots \\
f_M(t) \doteq \sum_{i=0}^N p_{M, i} t^i & T_{M-1} \leq t \leq T_M
\end{array}\right.
$$

其中由于是最小化snap，也就至少需要7次多项式，因此取$N=7$表示每段多项式包括常数项共有8个参数；$p_i$表示该段多项式中的第i个参数，其中$p_0$是常数项，$p_1$是一次项依此类推；$M$表示总共有M段轨迹，并且要求序列$\mathbf{T} = \{T_0,T_1,\cdots,T_M\}$必须都是已知的。


## 目标函数

$$
\begin{aligned}
& \min \quad J = \mathbf{P}^{\mathrm{T}} \mathbf{Q} \mathbf{P} \\
& \text{where} \quad \mathbf{P} = \begin{bmatrix}
\mathbf{p}_1 \\
\vdots \\
\mathbf{p}_{\mathrm{M}}
\end{bmatrix}, \quad \mathbf{Q} = \begin{bmatrix}
\mathbf{Q}_1 & 0 & \cdots & 0 \\
0 & \mathbf{Q}_2 & \cdots & 0 \\
\vdots & \vdots & \ddots & \vdots \\
0 & 0 & \cdots & \mathbf{Q}_{\mathrm{M}}
\end{bmatrix} \\
& \text{s.t.} \quad \mathbf{A}_{\mathrm{eq}}\mathbf{P} = \mathbf{d}_{\mathrm{eq}},\quad \mathbf{p}_i \in \mathbb{R}^8,\quad \mathbf{Q}_i \in \mathbb{R}^{8\times8} ,\quad i \in\{1,2,\cdots,M\}
\end{aligned}
$$

这是一个用来凸优化求解的模型，其中$\mathbf{p}_i$是第$i$段轨迹的参数，是8维的向量。$\mathbf{Q}_i$后续会提到。



## 代价函数

对于每一段轨迹$f(t)$都有，其中$N$为最高次数7，
$$
\begin{aligned}
& f(t)=\sum_{i=0}^{N} p_i t^i \\
& \Rightarrow f^{(4)}(t)=\sum_{i \geq 4} i(i-1)(i-2)(i-3) t^{i-4} p_i \\
& \Rightarrow\left(f^{(4)}(t)\right)^2=\sum_{i \geq 4, l \geq 4} i(i-1)(i-2)(i-3) l(l-1)(l-2)(l-3) t^{i+l-8} p_i p_l \\
& \Rightarrow J(T)=\int_{T_{j-1}}^{T_j}\left(f^4(t)\right)^2 d t=\sum_{i \geq 4, l \geq 4} \frac{i(i-1)(i-2)(i-3) j(l-1)(l-2)(l-3)}{i+l-7}\left(T_j^{i+l-7}-T_{j-1}^{i+l-7}\right) p_ip_l
\end{aligned}
$$
因此第$j$段轨迹的代价函数$J(T)$可以使用矩阵的形式表示为，
$$
J(T)=\int_{T_{j-1}}^{T_j}\left(f^4(t)\right)^2 d t = 
\left[\begin{array}{c}
\vdots \\
p_i \\
\vdots
\end{array}\right]^T

\left[\begin{array}{c}
&& \vdots &\\
&\cdots &\frac{i(i-1)(i-2)(i-3) l(l-1)(l-2)(l-3)}{i+l-7} \Delta T^{i+l-7}&\cdots\\
&&\vdots&\\
\end{array}\right]
\left[\begin{array}{c}
\vdots \\
p_l \\
\vdots
\end{array}\right] = \mathbf{p}_j^TQ_j\mathbf{p}_j
$$

其中$j$表示第$j$段轨迹，总共有$M$段；$\mathbf{p}_j\in\mathbb{R}^8$，而且对于矩阵$Q_j\in \mathbb{R}^{8\times8}$，由$J(T)$原来的表达式中可以看到，当且仅当，$i\ge 4 $且$l\ge4$的时候才有值，所以Q矩阵只有`Q(4:7,4:7)`为非零，并且在非零的部分按照$i$和$l$的值代入公式即可。

## 连续性约束

对于每一段轨迹，
$$
\mathrm{f(t)}=\mathrm{p}_7 \mathrm{t}^7+\mathrm{p}_6 \mathrm{t}^6+\mathrm{p}_5 \mathrm{t}^5+\mathrm{p}_4 \mathrm{t}^4+\mathrm{p}_3 \mathrm{t}^3+\mathrm{p}_2 \mathrm{t}^2+\mathrm{p}_1 \mathrm{t}^1+\mathrm{p}_0
$$


那么先定义一个$G(t)$，
$$
G(t) = \left[\begin{array}{cccccccc}
1 & \mathrm{t} & \mathrm{t}^2 & \mathrm{t}^3 & \mathrm{t}^4 & \mathrm{t}^5 & \mathrm{t}^6 & \mathrm{t}^7 \\
0 & 1 & 2 \mathrm{t} & 3 \mathrm{t}^2 & 4 \mathrm{t}^3 & 5 \mathrm{t}^4 & 6 \mathrm{t}^5 & 7 \mathrm{t}^6 \\
0 & 0 & 2 & 6 \mathrm{t} & 12 \mathrm{t}^2 & 20 \mathrm{t}^3 & 30 \mathrm{t}^4 & 42 \mathrm{t}^5 \\
0 & 0 & 0 & 6 & 24 \mathrm{t} & 60 \mathrm{t}^2 & 120 \mathrm{t}^3 & 210 \mathrm{t}^4
\end{array}\right]
$$
且，容易得，
$$
G{(t)}\cdot \mathbf{p} = \left[\begin{array}{cccccccc}
1 & \mathrm{t} & \mathrm{t}^2 & \mathrm{t}^3 & \mathrm{t}^4 & \mathrm{t}^5 & \mathrm{t}^6 & \mathrm{t}^7 \\
0 & 1 & 2 \mathrm{t} & 3 \mathrm{t}^2 & 4 \mathrm{t}^3 & 5 \mathrm{t}^4 & 6 \mathrm{t}^5 & 7 \mathrm{t}^6 \\
0 & 0 & 2 & 6 \mathrm{t} & 12 \mathrm{t}^2 & 20 \mathrm{t}^3 & 30 \mathrm{t}^4 & 42 \mathrm{t}^5 \\
0 & 0 & 0 & 6 & 24 \mathrm{t} & 60 \mathrm{t}^2 & 120 \mathrm{t}^3 & 210 \mathrm{t}^4
\end{array}\right] 
\cdot
\left[\begin{array}{l}
\mathrm{p}_0 \\
\mathrm{p}_1 \\
\mathrm{p}_2 \\
\mathrm{p}_3 \\
\mathrm{p}_4 \\
\mathrm{p}_5 \\
\mathrm{p}_6 \\
\mathrm{p}_7
\end{array}\right] =\left[\begin{array}{c}
f(t) \\
f^{(1)}(t) \\
f^{(2)}(t) \\
f^{(3)}(t)
\end{array}\right]=\left[\begin{array}{l}
\mathrm{x} \\
\mathrm{v} \\
\mathrm{a} \\
\mathrm{j}
\end{array}\right]
$$


导数约束要求满足边界条件以及经过中间的waypoints。那么先满足边界条件，
$$
\mathrm{G}(0) \bullet \mathrm{p}_{\mathbf{1}}=\left[\begin{array}{l}
\mathrm{x}_{\text {start }} \\
\mathrm{v}_{\text {start }} \\
\mathrm{a}_{\text {start }} \\
\mathrm{j}_{\text {start }}
\end{array}\right],\quad 
\mathbf{G}\left(\mathbf{T}_{\mathbf{M}}\right) \bullet \mathbf{p}_{\mathbf{M}}=\left[\begin{array}{l}
\mathrm{x}_{\text {stop }} \\
\mathrm{v}_{\text {stop }} \\
\mathrm{a}_{\text {stop }} \\
\mathrm{j}_{\text {stop }}
\end{array}\right]
$$


然后中间条件是， 
$$
\begin{aligned}
& \mathrm{G}\left(\mathrm{T}_1\right) \cdot \mathrm{p}_1=\mathrm{x}_1 \\
& \mathrm{G}\left(\mathrm{T}_2\right) \cdot \mathrm{p}_2=\mathrm{x}_2	\\
& \quad\quad\quad\quad \vdots\\
& \mathrm{G}\left(\mathrm{T}_{M-1}\right) \cdot \mathrm{p}_{M-1}=\mathrm{x}_{M-1}	\\
\end{aligned}
$$
综上所述可以构建为$\mathbf{A}_d\mathbf{p}_d = \mathbf{d}_d$的形式，其中$G_0(t)$表示$G(t)$的第一行，
$$
\mathbf{A}_{d}\mathbf{p}_{d} = \left[\begin{array}{ccc}
\mathrm{G}(0)_{4 \times 8} & \mathbf{0}_{\mathbf{4} \times \mathbf{8}} & \cdots &\mathbf{0}_{\mathbf{4} \times \mathbf{8}} \\

\mathbf{0}_{\mathbf{4} \times \mathbf{8}} & \mathbf{0}_{\mathbf{4} \times \mathbf{8}}& \cdots  & \mathrm{G}\left(\mathrm{T}_M\right)_{4 \times 8} \\

\mathrm{G}_0\left(\mathrm{T}_1\right)_{1 \times 8} & \mathbf{0}_{\mathbf{1} \times \mathbf{8}}&  \cdots  & \mathbf{0}_{\mathbf{1} \times \mathbf{8}} \\

\mathbf{0}_{\mathbf{1 \times 8}} & \mathrm{G}_0\left(\mathrm{T}_2\right)_{1 \times 8}& \cdots & \mathbf{0}_{\mathbf{1} \times \mathbf{8}}\\

\vdots & \vdots & \vdots & \vdots 
\end{array}\right]

\left[\begin{array}{c}
\mathrm{p}_{1,0} \\
\vdots \\
\mathrm{p}_{1,7} \\
\mathrm{p}_{2,0} \\
\vdots \\
\mathrm{p}_{2,7}\\
\vdots\\
\vdots\\
\mathrm{p}_{M,0} \\
\vdots \\
\mathrm{p}_{M,7} \\
\end{array}\right]=\left[\begin{array}{c}
\mathrm{x}_{\text {start }} \\
\mathrm{v}_{\text {start }} \\
\mathrm{a}_{\text {start }} \\
\mathrm{j}_{\text {start }} \\
\mathrm{x}_{\text {stop }} \\
\mathrm{v}_{\text {stop }} \\
\mathrm{a}_{\text {stop }} \\
\mathrm{j}_{\text {stop }} \\
\mathrm{x}_1 \\
\mathrm{x}_2 \\
\vdots\\
\mathrm{x}_{M-1}\\
\end{array}\right]
$$

其中，$\mathbf{p} \in \mathcal{R}^{8M\times 1}$，边缘条件的也就上面中的前两行矩阵大小是$8\times 8M$，然后中间点有M-1个，所以下面有M-1行，因此$A_d \in \mathcal{R}^{(M+7)\times8M}$。


## 导数约束

对于第$j$个点来说，其上一段轨迹的终点和下一段轨迹的起点的状态应该是要相同的，也就是说，要满足，
$$
\left[\begin{array}{l}
\mathrm{x}_j\left(\mathrm{~T}_j\right) \\
\mathrm{v}_j\left(\mathrm{~T}_j\right) \\
\mathrm{a}_j\left(\mathrm{~T}_j\right) \\
\mathrm{j}_j\left(\mathrm{~T}_j\right)
\end{array}\right]=\left[\begin{array}{l}
\mathrm{x}_{j+1}(0) \\
\mathrm{v}_{j+1}(0) \\
\mathrm{a}_{j+1}(0) \\
\mathrm{j}_{j+1}(0)
\end{array}\right]
$$


再引用之前定义过的$G(t)$函数，可以简化为，
$$
\left[\begin{array}{ll}
\mathrm{G}_{j}\left(\mathrm{T}_j\right) &  - \mathrm{G}_{j+1}(0)
\end{array}\right]\left[\begin{array}{l}
\mathbf{p}_j \\
\mathbf{p}_{j+1}
\end{array}\right]=\mathbf{0}
$$


因此可以将所有轨迹的连续性约束补充合并为$\mathbf{A}_c\mathbf{p}_c = \mathbf{d_c}$的形式，
$$
\mathbf{A}_{c}\mathbf{p}_{c} =\left[\begin{array}{ccc}
\mathrm{G}_1\left(\mathrm{T}_1\right)_{4 \times 8} & -\mathrm{G}_2(0)_{4 \times 8} & \mathbf{0}_{4 \times 8} &\cdots & \mathbf{0}_{4 \times 8}\\

\mathbf{0}_{4 \times 8} & \mathrm{G}_2\left(\mathrm{T}_2\right)_{4 \times 8} & -\mathrm{G}_3(0)_{4 \times 8}& \cdots & \mathbf{0}_{4 \times 8}\\

 \mathbf{0}_{4 \times 8} & \mathbf{0}_{4 \times 8}& \mathrm{G}_3\left(\mathrm{T}_3\right)_{4 \times 8} &\cdots & \mathbf{0}_{4 \times 8} \\
 
 \vdots & \vdots & \vdots & \ddots & \vdots\\
 
  \mathbf{0}_{4 \times 8}  & \mathbf{0}_{4 \times 8}  & \cdots &\mathrm{G}_{M-1}\left(\mathrm{T}_{M-1}\right)_{4 \times 8} &-\mathrm{G}_M\left(0\right)_{4 \times 8}
\end{array}\right]\left[\begin{array}{c}
\mathrm{p}_{1,0} \\
\vdots \\
\mathrm{p}_{1,7} \\
\mathrm{p}_{2,0} \\
\vdots \\
\mathrm{p}_{2,7} \\ 
\vdots\\
\vdots\\
\mathrm{p}_{M,0} \\
\vdots \\
\mathrm{p}_{M,7}
\end{array}\right]=[\mathbf{0}]_{4M\times 1}
$$


## 合并所有约束

其中导数约束总共有M+1个点，那么状态矩阵$\mathbf{A}_d \in \mathbb{R}^{4(M+7) \times 8M}$，而对于连续性约束，总共有M段轨迹，那么状态矩阵$\mathbf{A}_c \in \mathbb{R}^{4M \times 8M}$，因此可以纵向合并状态矩阵为$\mathbf{A}$，和$\mathbf{p}$，因此可以得到总的约束为，
$$
\mathbf{A}\mathbf{p } = \mathbf{d}
$$


## 模型求解

根据得到的最小化的代价函数$J(t)$以及约束条件$\mathbf{A p}=\mathbf{d}$，就可以放到`qp`求解器中直接算到最优的状态向量$\mathbf{p}$，表示的是每段轨迹的多项式参数，然后根据参数画出曲线即可得到minimum snap的轨迹。