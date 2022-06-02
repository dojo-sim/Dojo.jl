# using Pkg
# Pkg.develop(path=joinpath(@__DIR__, "../../DojoEnvironments"))
# Pkg.develop(path=joinpath(@__DIR__, "../.."))
# Pkg.activate(joinpath(@__DIR__, ".."))
# Pkg.instantiate()

# ## visualizer
vis = Visualizer()
open(vis)

# ## setup
using Dojo
using IterativeLQR
using LinearAlgebra
using FiniteDiff
using DojoEnvironments
using Plots
using StaticArrays


################################################################################
# Continuation
################################################################################
function reset!(env::Environment; rtol=1e-4, btol=1e-3, undercut=2.0)
    env.opts_step.rtol = rtol
    env.opts_step.btol = btol
    env.opts_step.undercut = 5.0
    env.opts_grad.rtol = rtol
    env.opts_grad.btol = btol
    env.opts_grad.undercut = undercut
    return nothing
end

function continuation_callback!(solver::Solver, env::Environment; ρ=1.5)
    # contact smoothness continuation
    env.opts_step.rtol = max(1e-6, env.opts_step.rtol/ρ)
    env.opts_step.btol = max(1e-4, env.opts_step.btol/ρ)
    env.opts_grad.rtol = max(1e-6, env.opts_grad.rtol/ρ)
    env.opts_grad.btol = max(1e-4, env.opts_grad.btol/ρ)

    # visualize current policy
    ū = solver.problem.actions
    x̄ = IterativeLQR.rollout(model, x1, ū)
    DojoEnvironments.visualize(env, x̄)

    println("r_tol $(scn(env.opts_grad.rtol))  " *
        "κ_tol $(scn(env.opts_grad.btol))")
    return nothing
end

################################################################################
# ## system
################################################################################
gravity = -9.81
timestep = 0.01
friction_coefficient = 0.4
damper = 1.0*1.5
spring = 0.0*5.0
env = DojoEnvironments.get_environment(:quadruped,
    representation=:minimal,
    timestep=timestep,
    contact_body=false,
    gravity=gravity,
    friction_coefficient=friction_coefficient,
    # damper=damper,
    # spring=spring,
    infeasible_control=true,
    vis=vis)

# ## dimensions
n = env.num_states
m = env.num_inputs
nu_infeasible = 6

# ## template
include(joinpath(@__DIR__, "../../DojoEnvironments/src",
    "quadruped/methods/template.jl"))


################################################################################
# ## simulation test
################################################################################
mech = get_mechanism(:quadruped,
    contact_body=false,
    timestep=timestep,
    gravity=gravity,
    friction_coefficient=friction_coefficient,
    damper=damper,
    limits=false,
    spring=spring)

initialize!(mech, :quadruped, body_position=[0,0,0.0])
u_hover = 0.0*[0.02;0;0.6; 0;0;0; zeros(12)]
function ctrl!(m, k; u=u_hover)
    nu = input_dimension(m)
    set_input!(m, SVector{nu}(u))
end

storage = simulate!(mech, 0.27, ctrl!,
    record=true,
    verbose=true,
    opts=SolverOptions(rtol=1e-5, btol=1e-4, undercut=5.0, verbose=false),
    )
Dojo.visualize(mech, storage, vis=env.vis)

# z_rest = get_maximal_state(mech)


################################################################################
# ## reference trajectory
################################################################################
N = 1
initialize!(env.mechanism, :quadruped)
xref = quadruped_trajectory(env.mechanism,
    r=0.00,
    z=0.25;
    Δx=-0.01,
    Δfront=0.05,
    width_scale=0.0,
    height_scale=1.0,
    N=1,
    Ncycles=N)
zref = [minimal_to_maximal(env.mechanism, x) for x in xref]
DojoEnvironments.visualize(env, xref)

# ## horizon
T = length(zref)

################################################################################
# ## ILQR problem
################################################################################
# ## model
dyn = IterativeLQR.Dynamics(
    (y, x, u, w) -> dynamics(y, env, x, u, w),
    (dx, x, u, w) -> dynamics_jacobian_state(dx, env, x, u, w),
    (du, x, u, w) -> dynamics_jacobian_input(du, env, x, u, w),
    n, n, m)

model = [dyn for t = 1:T-1]

# ## rollout
x1 = deepcopy(xref[1])
ū = [u_hover for t = 1:T-1]

x̄ = IterativeLQR.rollout(model, x1, ū)
DojoEnvironments.visualize(env, x̄)

# ## objective
############################################################################
qt = [0.3; 0.05; 0.05;
    5e-2 * ones(3);
    1e-0 * ones(3);
    1e-0 * ones(3);
    fill([2, 1e-0], 12)...]
ots = [(x, u, w) -> transpose(x - xref[t]) * Diagonal(timestep * qt) * (x - xref[t]) +
# transpose(u - u_hover) * Diagonal(timestep * 0.01 * ones(m)) * (u - u_hover) for t = 1:T-1]
    transpose(u) * Diagonal(timestep * 1.0 * ones(m)) * u for t = 1:T-1]
oT = (x, u, w) -> transpose(x - xref[end]) * Diagonal(timestep * qt) * (x - xref[end])


cts = [IterativeLQR.Cost(ot, n, m) for ot in ots]
cT = IterativeLQR.Cost(oT, n, 0)
obj = [cts..., cT]


# ## constraints
############################################################################
ul = -1.0 * 1e-3*ones(nu_infeasible)
uu = +1.0 * 1e-3*ones(nu_infeasible)

function contt(x, u, w)
    [
        1e-1 * (ul - u[1:nu_infeasible]);
        1e-1 * (u[1:nu_infeasible] - uu);
    ]
end

function goal(x, u, w)
    Δ = 1e-0 * (x - xref[end])
    return Δ
end

con_policyt = IterativeLQR.Constraint(contt, n, m, indices_inequality=collect(1:2nu_infeasible))
con_policyT = IterativeLQR.Constraint(goal, n, 0)

cons = [[con_policyt for t = 1:T-1]..., con_policyT]


# ## solver
options = Options(line_search=:armijo,
        max_iterations=50,
        max_dual_updates=12,
        min_step_size=1e-5,
        objective_tolerance=1e-3,
        lagrangian_gradient_tolerance=1e-3,
        constraint_tolerance=1e-4,
        initial_constraint_penalty=1e-1,
        scaling_penalty=10.0,
        max_penalty=1e4,
        verbose=true)

s = IterativeLQR.Solver(model, obj, cons, options=options)

IterativeLQR.initialize_controls!(s, ū)
IterativeLQR.initialize_states!(s, x̄)


# ## solve
local_callback!(solver::IterativeLQR.Solver) = continuation_callback!(solver, env)
reset!(env, rtol=1e-3, btol=1e-3)
@time IterativeLQR.constrained_ilqr_solve!(s, augmented_lagrangian_callback! = local_callback!)

# ## solution
x_sol, u_sol = IterativeLQR.get_trajectory(s)

# ## visualize
x_view = [[x_sol[1] for t = 1:15]..., x_sol..., [x_sol[end] for t = 1:15]...]
DojoEnvironments.visualize(env, x_view)

################################################################################
# TVLQR
################################################################################
N = 50
x_tv = [fill(x_sol[1:end-1], N)...; [x_sol[end]]]
u_tv = [fill(u_sol, N)...;]

env.opts_grad = SolverOptions(rtol=1e-4, btol=1e-3)
K_tv, P_tv = tvlqr(x_tv, u_tv, env;
        q_tracking=[10;10;10;
            1e0 * ones(3);
            1e0 * ones(3);
            1e0 * ones(3);
            fill([1e0, 1e0], 12)...],
        r_tracking=env.mechanism.timestep * 50 * ones(length(u_tv[1])))
env.opts_grad = SolverOptions(rtol=1e-5, btol=1e-5)


# K_tv, P_tv = tvlqr(x_tv, u_tv, env;
#         q_tracking=[0.3; 0.05; 0.05;
#             5e-1 * ones(3);
#             1e-0 * ones(3);
#             1e-0 * ones(3);
#             fill([2, 1e-0], 12)...],
#         r_tracking=env.mechanism.timestep * 1e1 * ones(length(u_tv[1])))

nu = input_dimension(mech)
nx = minimal_dimension(mech)
plot(hcat([reshape(K, nu*nx) for K in K_tv]...)', legend=false)

eigvals(P_tv[1])
cond(P_tv[1])
K_tv[1][18,35]
K_tv[1][18,36]

K_tv[1][17,33]
K_tv[1][17,34]

plot(hcat(x_sol...)')
plot(hcat(u_sol...)')

plot(hcat(x_tv...)')
plot(hcat(u_tv...)')

################################################################################
# Test Policy
################################################################################

# initialize!(mech, :quadruped, body_position=[0,0,0.])
# set_maximal_state!(mech, z_rest)
set_maximal_state!(mech, deepcopy(zref[1]))
function ctrl!(mechanism, k)
    nu = input_dimension(mechanism)
    x = get_minimal_state(mechanism)
    u = u_sol[1]/1.0 + 1.0*K_tv[1] * (x_sol[1] - x)
    u = [zeros(6); u[7:18]]
    set_input!(mechanism, SVector{nu}(u))
end

storage = simulate!(mech, 5.0, ctrl!,
    record=true,
    verbose=true,
    opts=SolverOptions(rtol=1e-5, btol=1e-4, undercut=5.0, verbose=false),
    )
Dojo.visualize(mech, storage, vis=env.vis, build=false)

################################################################################
# CIMPC compat
################################################################################

mutable struct TVLQRPolicy114{T}
    K::Vector{Matrix{T}}
    x::Vector{Vector{T}}
    u::Vector{Vector{T}}
    timestep::T
    H::Int
end

policy = TVLQRPolicy114(K_tv[1:T-1], x_sol[1:T-1], u_sol[1:T-1], timestep, T-1)

function exec_policy(p::TVLQRPolicy114{T}, x::Vector{T}, t::T) where {T}
    timestep = p.timestep
    i = floor(t/timestep) % H + 1
    u = p.u[i] + p.K[i] * (p.x[i] - x)
    return u / timestep # force
end

JLD2.jldsave(joinpath(@__DIR__, "tvlqr_policies", "standing_tvlqr_policy.jld2"),
    policy=policy,
    K=K_tv[1:T-1],
    x=x_sol[1:T-1],
    u=u_sol[1:T-1],
    timestep=timestep,
    H=T-1)

file = JLD2.jldopen(joinpath(@__DIR__, "tvlqr_policies", "standing_tvlqr_policy.jld2"))
policy = file["policy"]
K = file["K"]
x = file["x"]
u = file["u"]
timestep = file["timestep"]
H = file["H"]
JLD2.close(file)

# H = 4
# for t in 0:0.0035:0.3
#     @show floor(t/timestep) % H + 1
# end
