################################################################################
# Development
################################################################################
using Dojo
using Test

# Open visualizer
vis = Visualizer()
open(vis)

# pendulum
mechanism = get_mechanism(:npendulum, timestep = 0.01, gravity = -9.81 * 0.0, Nb=5)#, spring = 100.0, damper = 5.0)
Random.seed!(100)
ϕ1 = 0.0#0.3π
initialize!(mechanism, :npendulum, ϕ1 = ϕ1)

# # sphere
# mechanism = get_mechanism(:sphere, timestep = 0.01, gravity = -9.81)
# initialize!(mechanism, :sphere)

# half cheetah
mechanism = get_mechanism(:halfcheetah, timestep=0.01, gravity=-9.81)
initialize!(mechanism, :halfcheetah)

# # atlas
# mechanism = get_mechanism(:atlas, timestep=0.01, gravity=-9.81, cf=0.5, damper=100.0, spring=1.0, contact=true)
# initialize_atlasstance!(mechanism, tran=[0,0,0.5], rot=[0.0,0.0,0.0])

storage = simulate!(mechanism, 1.0, record = true, verbose = false)
visualize(mechanism, storage, vis = vis)

## Maximal gradients 
maximal_dimension(mechanism)
minimal_dimension(mechanism)
z = get_maximal_state(mechanism)
x = get_minimal_state(mechanism)
u = zeros(control_dimension(mechanism))
maximal_to_minimal(mechanism, z) - x
minimal_to_maximal(mechanism, x) - z

M_fd = maximal_to_minimal_jacobian(mechanism, z)
M_a = maximal_to_minimal_jacobian_analytical(mechanism, z)
@test size(M_fd) == size(M_a)
norm(M_fd - M_a, Inf)

ibody = 1
N_fd = minimal_to_maximal_jacobian(mechanism, x)#[(ibody-1) * 13 .+ (1:13), :]
N_a = minimal_to_maximal_jacobian_analytical(mechanism, x)#[(ibody-1) * 13 .+ (1:13), :]
sum(N_fd)
sum(N_a)
@test size(N_fd) == size(N_a)
@test norm(N_fd - N_a, Inf) < 1.0e-6

M_fd * N_fd
M_a * N_a
@test abs(sum(diag(M_fd * N_fd)) - minimal_dimension(mechanism)) < 1.0e-6
@test abs(sum(diag(M_a * N_a)) - minimal_dimension(mechanism)) < 1.0e-6

N_a * inv(N_a' * N_a)# - M_a

id = reverse(mechanism.system.dfs_list)[1]
joint = get_node(mechanism, id)
off = 0#(mechanism.joints[1])
Ne = length(mechanism.joints)
nu = control_dimension(joint)
idx = collect(off .+ (1:(2nu)))

currentvals = minimal_coordinates(mechanism)
currentvels = minimal_velocities(mechanism)

function joint_position_velocity(mechanism, joint, θ) 
    nu = control_dimension(joint)
    x, q = set_joint_position!(mechanism, joint, θ[1:nu]) 
    v, ω = set_velocity!(mechanism, joint, θ[nu .+ (1:nu)])
    return [x; v; vector(q); ω]
end

function joint_positions_velocities(mechanism, joint, θ) 
    # get minimal 
    currentvals = minimal_coordinates(mechanism)
    currentvels = minimal_velocities(mechanism)

    # set minimal 
    zp = joint_position_velocity(mechanism, joint, θ) 

    # child joints 
    child_joints = get_child_joints(mechanism, joint)

    for node in child_joints 
        set_joint_position!(mechanism, node, currentvals[node.id]) 
        set_velocity!(mechanism, node, currentvels[node.id])
    end 
end

function joint_position_velocity_jacobian(mechanism, joint, θ)
    FiniteDiff.finite_difference_jacobian(y -> joint_position_velocity(mechanism, joint, y), θ) 
end


function joint_position_velocity(mechanism, joint, z, θ) 
    nu = control_dimension(joint)
    body_parent = get_body(mechanism, joint.parent_id)

    xp = z[1:3] 
    vp = z[4:6]
    qp = UnitQuaternion(z[7:10]..., false)
    ϕp = z[11:13]

    if body_parent.name != :origin
        set_position!(body_parent, x=xp, q=qp)
        set_velocity!(body_parent, v=vp, ω=ϕp)
    end

    x, q = set_joint_position!(mechanism, joint, θ[1:nu]) 
    v, ω = set_velocity!(mechanism, joint, θ[nu .+ (1:nu)])

    return [x; v; vector(q); ω]
end

G = zeros(maximal_dimension(mechanism), minimal_dimension(mechanism))

# root 
zp = joint_position_velocity(mechanism, joint, x[idx])
∂z∂θ = joint_position_velocity_jacobian(mechanism, joint, x[idx])
G[(joint.child_id - Ne - 1) * 13 .+ (1:13), idx] = ∂z∂θ
norm(N_fd[(joint.child_id - Ne - 1) * 13 .+ (1:13), idx] - ∂z∂θ)

child_joints = get_child_joints(mechanism, joint)
∂z∂z = Dict()
da = copy(∂z∂θ)
for node in child_joints 
    d = FiniteDiff.finite_difference_jacobian(b -> joint_position_velocity(mechanism, node, b, [currentvals[node.id]; currentvels[node.id]]), zp)
    da = d * da
    push!(∂z∂z, "$(node.child_id)_$(node.parent_id)" => d)
    zp = joint_position_velocity(mechanism, node, zp, [currentvals[node.id]; currentvels[node.id]])

    G[(node.child_id - Ne - 1) * 13 .+ (1:13), idx] = da
end
N_fd[:, idx]
da
norm(G[:, idx] - N_fd[:, idx], Inf)








