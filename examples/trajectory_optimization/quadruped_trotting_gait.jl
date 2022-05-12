using Pkg
# Pkg.develop(path=joinpath(@__DIR__, "../../DojoEnvironments"))
Pkg.activate(joinpath(@__DIR__, ".."))
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

# ## system
gravity = -9.81
timestep = 0.01
friction_coefficient = 0.8
damper = 5.0
spring = 0.0
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

# ## template
include(joinpath(@__DIR__, "../../DojoEnvironments/src",
    "quadruped/methods/template.jl"))


# ## dimensions
n = env.num_states
m = env.num_inputs
nu_infeasible = 6

# ## reference trajectory
N = 1
initialize!(env.mechanism, :quadruped)
xref = quadruped_trajectory(env.mechanism,
    r=0.08,
    z=0.29;
    Δx=-0.04,
    Δfront=0.10,
    width_scale=0.0,
    height_scale=1.0,
    N=32,
    Ncycles=N)
zref = [minimal_to_maximal(env.mechanism, x) for x in xref]
DojoEnvironments.visualize(env, xref)


# ## gravity compensation
## TODO: solve optimization problem instead
mech = get_mechanism(:quadruped,
    timestep=timestep,
    gravity=gravity,
    friction_coefficient=friction_coefficient,
    damper=damper,
    spring=spring)

initialize!(mech, :quadruped)
storage = simulate!(mech, 1.0,
    record=true,
    verbose=false)

Dojo.visualize(mech, storage, vis=env.vis)
u_control = [0;0;1; 0;0;0; zeros(12)]


# ## horizon
T = length(zref) + 1

# ## model
dyn = IterativeLQR.Dynamics(
    (y, x, u, w) -> dynamics(y, env, x, u, w),
    (dx, x, u, w) -> dynamics_jacobian_state(dx, env, x, u, w),
    (du, x, u, w) -> dynamics_jacobian_input(du, env, x, u, w),
    n, n, m)

model = [dyn for t = 1:T-1]

# ## rollout
x1 = xref[1]
ū = [u_control for t = 1:T-1]
x̄ = IterativeLQR.rollout(model, x1, ū)
DojoEnvironments.visualize(env, x̄)

# ## objective
############################################################################
qt = [0.3; 0.05; 0.05; 0.01 * ones(3); 0.01 * ones(3); 0.01 * ones(3); fill([0.2, 0.001], 12)...]
ots = [(x, u, w) -> transpose(x - xref[t]) * Diagonal(timestep * qt) * (x - xref[t]) +
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
        1e-0 * (ul - u[1:nu_infeasible]);
        1e-0 * (u[1:nu_infeasible] - uu);
    ]
end

function goal(x, u, w)
    Δ = x - xref[end]
    return Δ
end

con_policyt = IterativeLQR.Constraint(contt, n, m, indices_inequality=collect(1:2nu_infeasible))
con_policyT = IterativeLQR.Constraint(goal, n, 0)

cons = [[con_policyt for t = 1:T-1]..., con_policyT]


# ## solver
options = Options(line_search=:armijo,
        max_iterations=100,
        max_dual_updates=30,
        objective_tolerance=1e-3,
        lagrangian_gradient_tolerance=1e-3,
        constraint_tolerance=1e-3,
        initial_constraint_penalty=1e-1,
        scaling_penalty=10.0,
        max_penalty=1e4,
        verbose=true)
s = IterativeLQR.Solver(model, obj, cons, options=options)

IterativeLQR.initialize_controls!(s, ū)
IterativeLQR.initialize_states!(s, x̄)

# ## solve
@time IterativeLQR.solve!(s)

# ## solution
x_sol, u_sol = IterativeLQR.get_trajectory(s)

# ## visualize
open(env.vis)
x_view = [[x_sol[1] for t = 1:15]..., x_sol..., [x_sol[end] for t = 1:15]...]
DojoEnvironments.visualize(env, x_view)


mech = get_mechanism(:quadruped,
    timestep=timestep,
    contact_body=false,
    gravity=gravity,
    friction_coefficient=friction_coefficient,
    damper=damper,
    spring=spring)

x_sol
u_sol
for i = 1:T-1
    @show i
    set_minimal_state!(mech, x_sol[i])
    set_input!(mech, u_sol[i])
    z = minimal_to_maximal(mech, x_sol[i])
    u = u_sol[i]
    step!(mech, z, u, opts=SolverOptions(rtol=1e-6, btol=2e-4))
end

full_vector(mech.system)

minimal_to_maximal_jacobian(mech, rand(1))

mech.bodies[1]


z = minimal_to_maximal(mech, x_sol[1])
u = u_sol[1]
step!(mech, z, u, opts=SolverOptions(rtol=1e-6, btol=2e-4))

get_next_state(mech)
function action(mechanism::Mechanism)
    a = 0.0

    return a
end

a = action(mech)
