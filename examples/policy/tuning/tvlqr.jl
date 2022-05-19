using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
Pkg.instantiate()

# ## visualizer
vis = Visualizer()
open(vis)

# ## setup
using Dojo
using IterativeLQR
using LinearAlgebra
using FiniteDiff
using DojoEnvironments
using JLD2

# ## scripts
include(joinpath(module_dir(), "examples/policy/methods/continuation.jl"))
include(joinpath(module_dir(), "examples/policy/methods/tvlqr.jl"))
include(joinpath(module_dir(), "DojoEnvironments/src",
    "quadruped/methods/template.jl"))

################################################################################
# Load trajectory
################################################################################
file = JLD2.jldopen(joinpath(@__DIR__, "../data/trotting_forward.jl"))
file["x"]
file["u"]
JLD2.close(file)





maximum(abs.(vcat(u_sol...)))

dynamics_jacobian_state(dx, env, x, u, w)
dynamics_jacobian_input(du, env, x, u, w)
A = [zeros(n,n) for i = 1:T-1]
B = [zeros(n,m) for i = 1:T-1]
for i = 1:T-1
    dynamics_jacobian_state(A[i], env, x_sol[i], u_sol[i], zeros(0))
    dynamics_jacobian_input(B[i], env, x_sol[i], u_sol[i], zeros(0))
end
qt = [0.3; 0.05; 0.05;
    5e-1 * ones(3);
    1e-6 * ones(3);
    1e-6 * ones(3);
    fill([2, 1e-6], 12)...]

rt = timestep * 100.0 * ones(m)

Q = [Diagonal(qt) for i = 1:T]
R = [Diagonal(rt) for i = 1:T-1]

K_tv, P_tv = tvlqr(A,B,Q,R)
