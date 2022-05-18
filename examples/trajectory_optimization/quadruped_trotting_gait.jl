using Pkg
Pkg.develop(path=joinpath(@__DIR__, "../../DojoEnvironments"))
Pkg.develop(path=joinpath(@__DIR__, "../.."))
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


################################################################################
# Continuation
################################################################################
function reset!(env::Environment; rtol=1e-4, btol=1e-3, undercut=5.0)
    env.opts_step.rtol = rtol
    env.opts_step.btol = btol
    env.opts_step.undercut = undercut
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
# u_hover = [0.05;0;1.13; 0;-0.01;0; zeros(12)]
u_hover = [0.02;0;0.6; 0;0;0; zeros(12)]
function ctrl!(m, k; u=u_hover)
    nu = input_dimension(m)
    set_input!(m, SVector{nu}(u))
end

# Main.@elapsed storage = simulate!(mech, 0.6, ctrl!,
Main.@profiler storage = simulate!(mech, 0.6, ctrl!,
    record=true,
    verbose=true,
    opts=SolverOptions(rtol=1e-4, btol=1e-2, undercut=5.0, verbose=false),
    )
Dojo.visualize(mech, storage, vis=env.vis)


# xtest = deepcopy(xref[1])
# xtest[14] += 1
# ztest = minimal_to_maximal(mech, xtest)
# set_robot(vis, mech, ztest)

################################################################################
# ## reference trajectory
################################################################################
N = 1
initialize!(env.mechanism, :quadruped)
xref = quadruped_trajectory(env.mechanism,
    r=0.08,
    z=0.29;
    Δx=-0.04,
    Δfront=0.10,
    width_scale=0.0,
    height_scale=1.0,
    N=30,
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
    0.001 * ones(3);
    0.001 * ones(3);
    0.001 * ones(3);
    fill([0.2, 0.0001], 12)...]
ots = [(x, u, w) -> transpose(x - xref[t]) * Diagonal(timestep * qt) * (x - xref[t]) +
    transpose(u - u_hover) * Diagonal(timestep * 0.01 * ones(m)) * (u - u_hover) for t = 1:T-1]
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
        min_step_size=1e-4,
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
reset!(env)
@time IterativeLQR.constrained_ilqr_solve!(s, augmented_lagrangian_callback! = local_callback!)

# ## solution
x_sol, u_sol = IterativeLQR.get_trajectory(s)

# ## visualize
x_view = [[x_sol[1] for t = 1:15]..., x_sol..., [x_sol[end] for t = 1:15]...]
DojoEnvironments.visualize(env, x_view)




x1 = deepcopy(xref[1])
ū = vcat(fill(u_sol, 5)...)
x̄ = IterativeLQR.rollout(fill(dyn, length(ū)), x1, ū)
DojoEnvironments.visualize(env, x̄)










mech = get_mechanism(:quadruped,
    timestep=timestep,
    contact_body=false,
    gravity=gravity,
    friction_coefficient=friction_coefficient,
    damper=damper,
    spring=spring)


z_sim = get_maximal_state(storage)
x_sim = [maximal_to_minimal(mech, z) for z in z_sim]
H = length(z_sim)
u_sim = [zeros(m) for i=1:H-1]

for i = 1:H-1
    @show i
    set_minimal_state!(mech, x_sim[i])
    set_input!(mech, u_sim[i])
    z = minimal_to_maximal(mech, x_sim[i])
    u = u_sim[i]
    J = minimal_to_maximal_jacobian(mech, x_sim[i])
    # step!(mech, z, u, opts=SolverOptions(rtol=1e-6, btol=2e-4))
end

full_vector(mech.system)

z = z_sim[1]
u = u_sim[1]
step!(mech, z, u, opts=SolverOptions(rtol=1e-8, btol=1e-6))
z20 = get_maximal_state(mech)
z30 = get_next_state(mech)
x20 = maximal_to_minimal(mech, z20)
x30 = maximal_to_minimal(mech, z30)

function get_initial_configurations(mechanism::Mechanism, z::Vector{T}) where T
    configuration_indices = [1,2,3,4,5,6,13,15,17,19,21,23,25,27,29,31,33,35]
    set_maximal_state!(mechanism, z)

    # current configuration
    x = maximal_to_minimal(mechanism, z)
    q2 = x[configuration_indices]

    # previous configuration
    for body in mechanism.bodies
        x, q = previous_configuration(body.state)
        set_maximal_configurations!(body, x=x, q=q)
    end
    z = get_maximal_state(mechanism)
    x = maximal_to_minimal(mechanism, z)
    q1 = x[configuration_indices]
    return q1, q2
end


configuration_indices = [1,2,3,4,5,6,13,15,17,19,21,23,25,27,29,31,33,35]
q10, q20 = get_initial_configurations(mech, z20)
q30 = x30[configuration_indices]

u0 = zeros(m)
w0 = zeros(0)
h0 = 0.0
μ0 = 0.0

θ0 = [q10, q20, u0, w0, μ0, h0]
norm(full_vector(mech.system))



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




function constraint_jacobian_configuration(relative::Symbol, joint::Joint{T,Nλ,Nb},
        xa::AbstractVector, qa::Quaternion,
        xb::AbstractVector, qb::Quaternion,
        η) where {T,Nλ,Nb}

    ∇comp = szeros(T,Nb,7)
    ∇mincoord = minimal_coordinates_jacobian_configuration(relative, joint, xa, qa, xb, qb, attjac=false)
    ∇unlim = joint_constraint_jacobian_configuration(relative, joint, xa, qa, xb, qb, η)

    # return [∇comp; ∇mincoord; -∇mincoord; ∇unlim]
end

relative = :parent
joint = mech.joints[4].rotational
xa = srand(3)
xb = srand(3)
qa = Quaternion(normalize(srand(4))...)
qb = Quaternion(normalize(srand(4))...)
η = 0.0

constraint_jacobian_configuration(relative, joint, xa, qa, xb, qb, η)
@benchmark $constraint_jacobian_configuration($relative, $joint, $xa, $qa, $xb, $qb, $η)




function impulse_map(relative::Symbol, joint::Joint{T,Nλ,0},
    xa::AbstractVector, qa::Quaternion,
    xb::AbstractVector, qb::Quaternion,
    η) where {T,Nλ}
    J = constraint_jacobian_configuration(relative, joint, xa, qa, xb, qb, η)
    G = [[Diagonal(sones(T,3)) szeros(T,3,3)]; [szeros(4,3) LVᵀmat(relative == :parent ? qa : qb)]]
    return Diagonal([sones(3); 0.5 * sones(3)]) * transpose(J * G)
end

relative = :parent
joint = mech.joints[4].translational
xa = srand(3)
xb = srand(3)
qa = Quaternion(normalize(srand(4))...)
qb = Quaternion(normalize(srand(4))...)
η = 0.0


impulse_map(relative, joint, xa, qa, xb, qb, η)
@benchmark $impulse_map($relative, $joint, $xa, $qa, $xb, $qb, $η)






function joint_constraint_jacobian_configuration(relative::Symbol, joint::Joint{T},
    xa::AbstractVector, qa::Quaternion,
    xb::AbstractVector, qb::Quaternion, η) where {T}
    X, Q = displacement_jacobian_configuration(relative, joint, xa, qa, xb, qb, attjac=false)
    return constraint_mask(joint) * [X Q]
end


relative = :parent
joint = mech.joints[4].translational
xa = srand(3)
xb = srand(3)
qa = Quaternion(normalize(srand(4))...)
qb = Quaternion(normalize(srand(4))...)
η = 0.0
joint_constraint_jacobian_configuration(relative, joint, xa, qa, xb, qb, η)
@benchmark $joint_constraint_jacobian_configuration($relative, $joint, $xa, $qa, $xb, $qb, $η)




function displacement_jacobian_configuration(relative::Symbol, joint::Rotational{T},
        xa::AbstractVector, qa::Quaternion,
        xb::AbstractVector, qb::Quaternion;
        attjac::Bool=true, vmat=true) where T
    X = szeros(T, 3, 3)
    if relative == :parent
		Q = Lᵀmat(joint.axis_offset) * Rmat(qb) * Tmat()
		attjac && (Q *= LVᵀmat(qa))
    elseif relative == :child
		Q = Lᵀmat(joint.axis_offset) * Lᵀmat(qa)
		attjac && (Q *= LVᵀmat(qb))
	end
	# vmat && (Q = Vmat() * Q)
	return X, Q
end


relative = :parent
joint = mech.joints[4].rotational
xa = srand(3)
xb = srand(3)
qa = Quaternion(normalize(srand(4))...)
qb = Quaternion(normalize(srand(4))...)
η = 0.0
displacement_jacobian_configuration(relative, joint, xa, qa, xb, qb)
Main.@code_warntype displacement_jacobian_configuration(relative, joint, xa, qa, xb, qb)
@benchmark $displacement_jacobian_configuration($relative, $joint, $xa, $qa, $xb, $qb)




function displacement_jacobian_configuration(relative::Symbol, joint::Translational{T},
    xa::AbstractVector, qa::Quaternion,
    xb::AbstractVector, qb::Quaternion;
    attjac=true) where T

    vertices = joint.vertices

    if relative == :parent
        d = xb + vector_rotate(vertices[2], qb) - (xa + vector_rotate(vertices[1], qa)) # in the world frame
        X = -rotation_matrix(inv(qa))
        Q = -rotation_matrix(inv(qa)) * ∂rotation_matrix∂q(qa, vertices[1])
        Q += ∂rotation_matrix_inv∂q(qa, d)
        attjac && (Q *= LVᵀmat(qa))
    elseif relative == :child
        X = rotation_matrix(inv(qa))
        Q = rotation_matrix(inv(qa)) * ∂rotation_matrix∂q(qb, vertices[2])
        attjac && (Q *= LVᵀmat(qb))
    end

    return X, Q
end




displacement_jacobian_configuration

joint_constraint_jacobian_configuration
minimal_coordinates_jacobian_configuration
minimal_coordinates_jacobian_configuration
minimal_velocities_jacobian_configuration
minimal_velocities_jacobian_velocity
