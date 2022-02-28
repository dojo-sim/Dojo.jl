# Linearized friction

### Mathematical Model
Coulomb friction instantaneously maximizes the dissipation of kinetic energy between two objects in contact. For a single contact point, this physical phenomenon can be modeled by the following optimization problem,

$$\text{minimize}_{b} \quad v^T b \\
\text{subject to} \quad \|b\|_2v\leq \mu \gamma,$$

where $v \in \mathbf{R}^{2}$ is the tangential velocity at the contact point, $b \in \mathbf{R}^2$ is the friction force, and $\mu \in \mathbf{R}_{+}$ is the coefficient of friction between the two objects.



This problem is naturally a convex second-order cone program, and can be efficiently and reliably solved. However, classically, an approximate version of the above problem:

$$\text{minimize}_{\beta} \quad [v^T  -v^T] \beta, \\
\text{subject to} \quad \beta^T \mathbf{1} \leq \mu \gamma, \\
\beta \geq 0,$$

which satisfies the LCP formulation, is instead solved. Here, the friction cone is linearized and the friction vector, $\beta \in \mathbf{R}^{4}$, is correspondingly overparameterized and subject to additional non-negative constraints.

The optimality conditions of the above problem and constraints used in the LCP are:

$$[v^T  -v^T]^T + \psi \mathbf{1} - \eta = 0, \\
\mu \gamma -\beta^T \textbf{1} \geq 0,\\
\psi \cdot (\mu \gamma - \beta^T \textbf{1}) = 0, \\
\beta \circ \eta = 0, \\
\beta, \psi, \eta \geq 0,$$

where $\psi \in \mathbf{R}$ and $\eta \in \mathbf{R}^{4}$ are the dual variables associated with the friction cone and positivity constraints, respectively, and $\textbf{1}$ is a vector of ones.

### Constructor
[`LinearContact`](@ref)
