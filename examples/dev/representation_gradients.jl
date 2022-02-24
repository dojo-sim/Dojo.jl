################################################################################
# Development
################################################################################
using Dojo
using Test

# Open visualizer
vis = Visualizer()
open(vis)

# pendulum
mech = get_mechanism(:npendulum, timestep = 0.01, gravity = -9.81 * 0.0, Nb=5)#, spring = 100.0, damper = 5.0)
Random.seed!(100)
ϕ1 = 0.0#0.3π
initialize!(mech, :npendulum, ϕ1 = ϕ1)

[id for id in reverse(mech.system.dfs_list) if get_node(mech, id) isa JointConstraint]

mech.joints[1].parent_id
mech.joints[1].child_id
mech.joints[5].parent_id
mech.joints[5].child_id
mech.joints[6].parent_id
mech.joints[6].child_id
mech.joints[7].parent_id
mech.joints[7].child_id

mech.joints[2].parent_id
mech.joints[2].child_id
mech.joints[3].parent_id
mech.joints[3].child_id
mech.joints[4].parent_id
mech.joints[4].child_id		

# sphere
mech = get_mechanism(:sphere, timestep = 0.01, gravity = -9.81)
initialize!(mech, :sphere)

# half cheetah
mech = get_mechanism(:halfcheetah, timestep=0.01, gravity=-9.81)
initialize!(mech, :halfcheetah)

# atlas
mech = get_mechanism(:atlas, timestep=0.01, gravity=-9.81, friction_coefficient=0.5, damper=100.0, spring=1.0, contact=true)
initialize_atlasstance!(mech, tran=[0,0,0.5], rot=[0.0,0.0,0.0])

joint = mech.joints[1]
joint.name
n = control_dimension(joint)
off = 0
idx = collect(off .+ (1:(2n)))
child_ids = [id for id in recursivedirectchildren!(mech.system, joint.id) if get_node(mech, id) isa JointConstraint]

child_joints = unique([get_node(mech, id) for id in recursivedirectchildren!(mech.system, joint.id) if get_node(mech, id) isa JointConstraint])
child_joints = unique([id for id in recursivedirectchildren!(mech.system, joint.id) if get_node(mech, id) isa JointConstraint])
reverse(mech.system.dfs_list)
# mech = get_ant(;timestep=0.01, gravity=0.0)
# initialize_ant!(mech)

storage = simulate!(mech, 1.0, record = true, verbose = false)
visualize(mech, storage, vis = vis)

child_ids = [id for id in recursivedirectchildren!(mech.system, joint.id) if get_node(mech, id) isa JointConstraint]
child_joints = [get_node(mech, id) for id in recursivedirectchildren!(mech.system, joint.id) if get_node(mech, id) isa JointConstraint]
child_joints = get_child_joints(mech, joint)

joint.id
children(mech.system, 1)
mech.joints[1].parent_id
mech.joints[1].child_id
mech.joints[2].parent_id
mech.joints[2].child_id
mech.joints[3].parent_id
mech.joints[3].child_id
mech.joints[4].parent_id
mech.joints[4].child_id
mech.joints[5].parent_id

joint.child_id
for j in child_joints 
    @show j.id
    # @show j.parent_id 
    # @show j.child_id 
end

## Maximal gradients 
maximal_dimension(mech)
minimal_dimension(mech)
z = get_maximal_state(mech)
x = get_minimal_state(mech)
u = zeros(control_dimension(mech))
maximal_to_minimal(mech, z) - x
minimal_to_maximal(mech, x) - z

M_fd = maximal_to_minimal_jacobian(mech, z)
M_a = maximal_to_minimal_jacobian_analytical(mech, z)
@test size(M_fd) == size(M_a)
norm(M_fd - M_a, Inf)

ibody = 1
N_fd = minimal_to_maximal_jacobian(mech, x)#[(ibody-1) * 13 .+ (1:13), :]
N_a = minimal_to_maximal_jacobian_analytical(mech, x)#[(ibody-1) * 13 .+ (1:13), :]
sum(N_fd)
sum(N_a)
@test size(N_fd) == size(N_a)
@test norm(N_fd - N_a, Inf) < 1.0e-6

M_fd * N_fd
M_a * N_a
@test abs(sum(diag(M_fd * N_fd)) - minimal_dimension(mech)) < 1.0e-6
@test abs(sum(diag(M_a * N_a)) - minimal_dimension(mech)) < 1.0e-6

N_a * inv(N_a' * N_a)# - M_a

# # bodies
# joint = mech.joints[1]
# pbody = get_body(mech, joint.parent_id)
# cbody = get_body(mech, joint.child_id)

# # translational delta in pbody frame
# xθ = [1.0, 1.0]
# Δx = get_position_delta(joint.translational, pbody, cbody, xθ[SUnitRange(joint.minimal_index[1][1], joint.minimal_index[1][2])]) 

# # rotational delta in cbody frame
# Δq = get_position_delta(joint.rotational, pbody, cbody, xθ[SUnitRange(joint.minimal_index[2][1], joint.minimal_index[2][2])])


#######
off = 0#(mech.joints[1])
mechanism = mech
length(mech.joints)
id = reverse(mechanism.system.dfs_list)[1]

joint = mechanism.joints[id]
n = control_dimension(joint)
idx = collect(off .+ (1:(2n)))

child_joints = get_child_joints(mech, joint)

for j in child_joints 
    @show j.parent_id 
    @show j.child_id
end

child_ids = [id for id in recursivedirectchildren!(mech.system, joint.id) if get_node(mech, id) isa JointConstraint]
child_joints = [get_node(mech, id) for id in child_ids]
# joint.child_id
# child_joints[1].parent_id
# child_joints[1].child_id

# for j in child_joints 
#     @show j.parent_id 
#     @show j.child_id 
# end

function position_velocity(y) 
    mech = mechanism

    currentvals = minimal_coordinates(mech)
    currentvels = minimal_velocities(mech)

    set_joint_position!(mech, joint, y[1:n]) 
    set_minimal_velocities!(mech, joint, y[n .+ (1:n)])

    for node in child_joints
        set_joint_position!(mech, node, currentvals[node.id])
        set_minimal_velocities!(mech, node, currentvels[node.id])
    end

    return get_maximal_state(mech)
end 

ichild = joint.child_id - length(mech.joints)
norm(minimal_to_maximal_jacobian(mech, x)[(ichild - 1) * 13 .+ (1:13), idx] - FiniteDiff.finite_difference_jacobian(position_velocity, x[idx])[(ichild - 1) * 13 .+ (1:13), idx])

currentvals = minimal_coordinates(mech)
currentvels = minimal_velocities(mech)

function joint_position_velocity(mech, joint, θ) 
    n = control_dimension(joint)
    x, q = set_joint_position!(mech, joint, θ[1:n]) 
    v, ω = set_minimal_velocities!(mech, joint, θ[n .+ (1:n)])
    return [x; v; vector(q); ω]
end

D1 = FiniteDiff.finite_difference_jacobian(a -> joint_position_velocity(mech, joint, a), x[idx])
norm(D1 - minimal_to_maximal_jacobian(mech, x)[(ichild - 1) * 13 .+ (1:13), 1:2n])

xa, qa = set_joint_position!(mech, joint, x[idx][1:n]) 
va, ωa = set_minimal_velocities!(mech, joint, x[idx][n .+ (1:n)])
zp = [xa; va; vector(qa); ωa]

for node in child_joints
    set_joint_position!(mech, node, currentvals[node.id])
    set_minimal_velocities!(mech, node, currentvels[node.id])
end

# root 
∂z∂θ = joint_position_velocity_jacobian(mech, joint, x[idx])
Dz = []
Da = []
y = ∂z∂θ
for node in child_joints	
    ∂z∂z = FiniteDiff.finite_difference_jacobian(b -> joint_position_velocity(mech, node, b, [currentvals[node.id]; currentvels[node.id]]), zp)
    y = ∂z∂z * y
    push!(Da, y)
    push!(Dz, ∂z∂z)
   
    zp = joint_position_velocity(mech, node, zp, [currentvals[node.id]; currentvels[node.id]])
end

ichild + length(mech.joints)
child_joints[1].parent_id

norm(minimal_to_maximal_jacobian(mech, x)[(ichild - 1) * 13 .+ (1:13), 1:2n] - ∂z∂θ)



minimal_to_maximal_jacobian(mech, x)[(child_joints[1].child_id - length(mech.joints) - 1) * 13 .+ (1:13), 1:2n]
Da[1]
norm(minimal_to_maximal_jacobian(mech, x)[(child_joints[1].child_id - length(mech.joints) - 1) * 13 .+ (1:13), 1:2n] - Dz[1] * ∂z∂θ)
norm(minimal_to_maximal_jacobian(mech, x)[(child_joints[1].child_id - length(mech.joints) - 1) * 13 .+ (1:13), 1:2n] - Da[1])

norm(minimal_to_maximal_jacobian(mech, x)[(child_joints[2].child_id - length(mech.joints) - 1) * 13 .+ (1:13), 1:2n] - Dz[2] * Dz[1] * ∂z∂θ)
norm(minimal_to_maximal_jacobian(mech, x)[(child_joints[2].child_id - length(mech.joints) - 1) * 13 .+ (1:13), 1:2n] - Da[2])

norm(minimal_to_maximal_jacobian(mech, x)[(child_joints[3].child_id - length(mech.joints) - 1) * 13 .+ (1:13), 1:2n] - Dz[3] * Dz[2] * Dz[1] * ∂z∂θ)
norm(minimal_to_maximal_jacobian(mech, x)[(child_joints[3].child_id - length(mech.joints) - 1) * 13 .+ (1:13), 1:2n] - Da[3])

function joint_position_velocity(mech, joint, z, θ) 
    n = control_dimension(joint)

    body_parent = get_body(mech, joint.parent_id)
    xp = z[1:3] 
    vp = z[4:6]
    qp = UnitQuaternion(z[7:10]..., false)
    ϕp = z[11:13]

    set_maximal_configuration!(body_parent, x=xp, q=qp)
    set_maximal_velocity!(body_parent, v=vp, ω=ϕp)

    x, q = set_joint_position!(mech, joint, θ[1:n]) 
    v, ω = set_minimal_velocities!(mech, joint, θ[n .+ (1:n)])

    return [x; v; vector(q); ω]
end

function joint_position_velocity_jacobian(mech, joint, θ)
    FiniteDiff.finite_difference_jacobian(y -> joint_position_velocity(mech, joint, y), θ) 
end

function position_velocity_jacobian(θ) 
    G = zeros(maximal_dimension(mechanism), length(θ)) 

    mech = mechanism

    currentvals = minimal_coordinates(mech)
    currentvels = minimal_velocities(mech)

    x, q = set_joint_position!(mech, joint, θ[1:n]) 
    v, ω = set_minimal_velocities!(mech, joint, θ[n .+ (1:n)])
    zp = [x; v; vector(q); ω]

    for node in child_joints
        set_joint_position!(mech, node, currentvals[node.id])
        set_minimal_velocities!(mech, node, currentvels[node.id])
    end

    # root 
    ∂z∂θ = joint_position_velocity_jacobian(mech, joint, θ)
    G[(joint.child_id - 1 - Ne) * 13 .+ (1:13), :] = ∂z∂θ

    for node in child_joints		
        # ∂z∂θi = joint_position_velocity_jacobian(mech, node, [currentvals[node.id]; currentvels[node.id]])
        # ∂z∂z = FiniteDiff.finite_difference_jacobian(b -> joint_position_velocity(mech, node, b, [currentvals[node.id]; currentvels[node.id]]), zp)
        # # @show ∂z∂z
        # ∂z∂θ = ∂z∂z * (∂z∂θi)
        # G[(node.child_id - 1 - Ne) * 13 .+ (1:13), :] = ∂z∂θ

        # zp = joint_position_velocity(mech, node, zp, θ)
    end

    return G
end

J[:, idx] = position_velocity_jacobian(a[idx])