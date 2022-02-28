# Nonlinear friction

### Mathematical Model

In contrast to the LCP approach, we utilize the optimality conditions of in a form amenable to a primal-dual interior-point solver. The associated cone program is,

$$\text{minimize}  \quad [0 \: v^T] \beta\\
\text{subject to} \quad \beta_{(1)} = \mu \gamma, \\
\beta \in \mathcal{Q}^3,$$

where subscripts indicate vector indices and the $n$-dimensional second-order cone $\mathcal{Q}^n$ is defined by:

$$\mathcal{Q}^n = \{(a_{(1)}, a_{(2:n)}) \in \mathbf{R} \times \mathbf{R}^{n-1}\, | \, \|a_{(2:n)}\|_2 \leq a_{(1)} \}$$

The relaxed optimality conditions for the above problem in interior-point form are:

$$v - \eta_{(2:3)} = 0, \\
\beta_{(1)} - \mu \gamma = 0, \\
\beta \circ \eta =  \kappa \mathbf{e}, \\
\beta, \eta \in \mathcal{Q}^3,$$

with dual variable $\eta \in \mathcal{Q}^3$ associated with the second-order-cone constraints, and central-path parameter, $\kappa \in \mathbf{R}_{+}$.
The second-order-cone product is:

$$\beta \circ \eta = (\beta^T \eta, \beta_{(1)} \eta_{(2:n)} + \eta_{(1)} \beta_{(2:n)}),$$

and,

$$\mathbf{e} = (1, 0, \dots, 0),$$

is its corresponding identity element. Friction is recovered from the solution: $b = \beta^*_{(2:3)}$. The benefits of this model are increased physical fidelity and fewer optimization variables, without substantial increase in computational cost.

Nonlinear complementarity problem: To simulate a system represented in maximal coordinates that experiences contact, a solver aims to satisfy the following relaxed feasibility problem:

$$\text{find} \quad z_{+}, w, \gamma, \beta^{(1:P)}, \eta^{(1:P)}, s\\
\text{s.t.} \quad f(z_{-}, z, z_{+}, w) + B(z) u + C(z)^T \lambda = 0 \\
s - \phi(z_{+}) = 0, \\
\gamma \circ s = \kappa \textbf{1},\\
\beta^{(i)} \circ \eta^{(i)} = \kappa \mathbf{e}, \quad i = 1, \dots, P,\\
v^{(i)}(z, z_{+}) - \eta_{(2:3)}^{(i)} = 0, \quad i = 1, \dots, P, \\
\beta^{(i)}_{(1)} - \mu^{(i)} \gamma^{(i)} = 0, \quad i = 1, \dots, P,\\
\gamma, s \geq 0,\\
\beta^{(i)}, \eta^{(i)} \in \mathcal{Q}^3, \quad i = 1, \dots, P,$$

where $u \in \mathbf{R}^m$ is the control input at the current time step, $\lambda = (\beta^{(1)}_{(2:3)}, \gamma^{(1)}, \dots, \beta^{(P)}_{(2:3)}, \gamma^{(P)}) \in \mathbf{\Lambda}$ is the concatenation of impact and friction impulses, $B : \mathbf{Z} \rightarrow \mathbf{R}^{6N \times m}$ is the input Jacobian mapping control inputs into maximal coordinates, $C : \mathbf{Z} \rightarrow \mathbf{R}^{\text{dim}(\mathbf{\Lambda}) \times 6N}$ is a contact Jacobian mapping between maximal coordinates and contact surfaces, $s \in \mathbf{R}^P$ is a slack variable introduced for convenience, and $v^{(i)} : \mathbf{Z} \times \mathbf{Z} \rightarrow \mathbf{R}^2$ is the tangential velocity at contact point $i$. Joint limits and internal friction are readily incorporated into this problem formulation.

### Constructor
```@docs
NonlinearContact
```
