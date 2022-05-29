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
timestep = 0.02
friction_coefficient = 0.8
damper = 0.5
spring = 5.0
env = get_environment(:quadruped,
    representation=:minimal,
    timestep=timestep,
    contact_body=false,
    gravity=gravity,
    friction_coefficient=friction_coefficient,
    damper=damper,
    spring=spring,
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
    spring=spring)

initialize!(mech, :quadruped, body_position=[0,0,0.0])
u_hover = 0.0*[0.02;0;0.6; 0;0;0; zeros(12)]
function ctrl!(m, k; u=u_hover)
    nu = input_dimension(m)
    set_input!(m, SVector{nu}(u))
end

storage = simulate!(mech, 1.0, ctrl!,
    record=true,
    verbose=true,
    opts=SolverOptions(rtol=1e-5, btol=1e-4, undercut=5.0, verbose=false),
    )
Dojo.visualize(mech, storage, vis=env.vis)

################################################################################
# ## reference trajectory
################################################################################
N = 1
initialize!(env.mechanism, :quadruped)
xref = quadruped_trajectory(env.mechanism,
    r=0.00,
    z=0.29;
    Δx=-0.04,
    Δfront=0.10,
    width_scale=0.0,
    height_scale=1.0,
    N=2,
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
    1e-3 * ones(3);
    1e-3 * ones(3);
    fill([2, 1e-3], 12)...]
ots = [(x, u, w) -> transpose(x - xref[t]) * Diagonal(timestep * qt) * (x - xref[t]) +
# transpose(u - u_hover) * Diagonal(timestep * 0.01 * ones(m)) * (u - u_hover) for t = 1:T-1]
    transpose(u) * Diagonal(timestep * 0.01 * ones(m)) * u for t = 1:T-1]
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
        1e-3 * (ul - u[1:nu_infeasible]);
        1e-3 * (u[1:nu_infeasible] - uu);
    ]
end

function goal(x, u, w)
    Δ = 1e-2 * (x - xref[end])[[1:6;13:2:36]]
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
