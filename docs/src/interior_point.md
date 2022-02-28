# Algorithm

### Overview
To simulate the system forward in time, we need to solve a Nonlinear Complementarity Problem (NCP) at each time step. To efficiently and reliably satisfy the NCP, we developed a custom primal-dual interior-point solver for NCPs with cone constraints and quaternions. The algorithm is largely based upon Mehrotra's predictor-coorector algorithm, while borrowing practical numerical features from CVXOPT to handle cones and non-Euclidean optimization to handle quaternions. We also introduce heuristics that further improve reliability and overall performance of the solver for our simulation-step NCPs.

The primary advantages of this algorithm are the correction to the classic Newton step, which can greatly reduce the iterations required by the solver (often halving the total number of iterations), and feedback on the problem's central-path parameter that helps avoid premature ill-conditioning and adaptively drives the complementarity violation to zero in order to reliably simulate hard contact.

### Problem formulation
The solver aims to satisfy instantiations of the following problem:

$$\text{find}\quad  x, y, z \\
\text{subject to}\quad  c(x, y, z; \theta) = 0, \\
y^{(i)} \circ z^{(i)} = \kappa \mathbf{e}, \quad i = 1,\dots,n, \\
y^{(i)}, z^{(i)} \in \mathcal{K}, \quad i = 1,\dots, n,$$

with decision variables $x \in \mathbf{R}^k$ and $y, z \in \mathbf{R}^m$, equality-constraint set $c : \mathbf{R}^k \times \mathbf{R}^m \times \mathbf{R}^m \times \mathbf{R}^l \rightarrow \mathbf{R}^h$, problem data $\theta \in \mathbf{R}^l$; and where $\mathcal{K}$ is the Cartesian product of $n$ total positive-orthant and second-order cones. The variables are partitioned: $x = (x^{(1)}, \dots, x^{(p)})$, where $i = 1$ are Euclidean variables and $i = 2, \dots, p$ are each quaternion variables; and $y = (y^{(1)}, \dots, y^{(n)})$, $z = (z^{(1)}, \dots, z^{(n)})$, where $j = 1$ is the positive-orthant and the remaining $j = 2, \dots, n$ are second-order cones. For convenience, we denote $w = (x, y, z)$.

The algorithm aims to satisfy a sequence of relaxed problems with $\kappa > 0$ and $\kappa \rightarrow 0$ in order to reliably converge to a solution of the original problem (i.e., $\kappa = 0$). This continuation approach helps avoid premature ill-conditioning and is the basis for numerous convex and non-convex general-purpose interior-point solvers.

### Violation metrics:
Two metrics are used to measure progress:
The constraint violation,

$$r_{\text{vio}} = \| c(w; \theta) \|_{\infty},$$

and complementarity violation,

$$b_{\text{vio}} = {\text{max}}_i \{\| y^{(i)} \circ z^{(i)} \|_{\infty}\}.$$

The NCP is considered solved when $r_{\text{vio}} < r_{\text{tol}}$ and $b_{\text{vio}} < b_{\text{tol}}$.
!!! info "solver options"
    Both `r_tol` and `b_tol` are options that can easily be accessed and modified via [`SolverOptions`](@ref).

### Newton Steps
The main loop of the solver performs Newton's method on the equality-constraint set $c$ and the bilinear constraints. The solver typically converges in about 10 iterations.
!!! info "solver options"
    The maximal number of Newton's iterations `max_iter` can be set via [`SolverOptions`](@ref).

### Line Search
Newton's method provides a search direction, then we perform a line search along this search direction to detemine the step length $\alpha$. We use a backtracking line search that accepts the step whenever it decreases $c_{\text{vio}}$ or $b_{\text{vio}}$.
The line search starts with a step $\alpha=1$, if the step acceptance conditions are not met the step is decreased geometrically:
$$\alpha \leftarrow \alpha \times s $$. The line search takes at most `max_ls` backtracking steps.

!!! info "solver options"
    The scaling parameter $s$ is called `ls_scale` and the maximum number of line search iteration `max_ls` can be set via [`SolverOptions`](@ref).
