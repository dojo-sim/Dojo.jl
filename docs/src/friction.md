# Friction

## Mathematical Model
Coulomb friction instantaneously maximizes the dissipation of kinetic energy between two objects in contact. For a single contact point, this physical phenomenon can be modeled by the following optimization problem,

$$\text{minimize}_{b} v^T b \\
\text{subject to} \|b\|_2v\leq \mu \gamma,$$

where $v \in \mathbf{R}^{2}$ is the tangential velocity at the contact point, $b \in \mathbf{R}^2$ is the friction force, and $\mu \in \mathbf{R}_{+}$ is the coefficient of friction between the two objects.



This problem is naturally a convex second-order cone program, and can be efficiently and reliably solved. However, classically, an approximate version of the above problem:

$$\text{minimize}_{\beta} [v^T  -v^T] \beta, \\
\text{subject to} \beta^T 1 \leq \mu \gamma, \\
\beta \geq 0,$$

which satisfies the LCP formulation, is instead solved. Here, the friction cone is linearized (Fig. \ref{friction_cones}) and the friction vector, $\beta \in \mathbf{R}^{4}$, is correspondingly overparameterized and subject to additional non-negative constraints \cite{stewart1996implicit}.

The optimality conditions of \eqref{mdp_linear} and constraints used in the LCP are:
\begin{align}
& \begin{bmatrix} v^T & -v^T\end{bmatrix}^T + \psi \mathbf{1} - \eta = 0, \label{mdp_eq}\\
& \mu \gamma -\beta^T \textbf{1} \geq 0, \label{mdp_cone}\\
& \psi \cdot (\mu \gamma - \beta^T \textbf{1}) = 0, \label{mdp_cone_comp}\\
& \beta \circ \eta = 0, \label{mdp_friction_comp}\\
& \beta, \psi, \eta \geq 0, \label{mdp_ineq}
\end{align}
where $\psi \in \mathbf{R}$ and $\eta \in \mathbf{R}^{4}$ are the dual variables associated with the friction cone and positivity constraints, respectively, and $\textbf{1}$ is a vector of ones.


## Linearized friction cone


## Non-linear friction cone
