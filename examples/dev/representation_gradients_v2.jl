################################################################################
# Development
################################################################################
using Dojo
using Test

# Open visualizer
vis=visualizer()
open(vis)

# pendulum
mechanism = get_mechanism(:npendulum, timestep=0.01, gravity=-9.81 * 0.0, Nb=10)#, spring = 100.0, damper = 5.0)
Random.seed!(100)
ϕ1 = 0.3π
initialize!(mechanism, :npendulum, ϕ1 = ϕ1)

# sphere
mechanism = get_mechanism(:sphere, timestep=0.01, gravity=-9.81)
initialize!(mechanism, :sphere)

# half cheetah
mechanism = get_mechanism(:halfcheetah, timestep=0.01, gravity=-9.81)
initialize!(mechanism, :halfcheetah)

# atlas
mechanism = get_mechanism(:atlas, timestep=0.01, gravity=-9.81, friction_coefficient=0.5, damper=100.0, spring=1.0, contact=true)
initialize_atlasstance!(mechanism, tran=[0,0,0.5], rot=[0.0,0.0,0.0])

function ctrl!(mechanism, k) 
    set_input!(mechanism, 0.1 * srand(input_dimension(mechanism)))
end


storage = simulate!(mechanism, 1.0, ctrl!, record=true, verbose=false)
# visualize(mechanism, storage, vis=vis)

## Maximal gradients 
maximal_dimension(mechanism)
minimal_dimension(mechanism)
z = get_maximal_state(mechanism)
x = get_minimal_state(mechanism)
u = zeros(input_dimension(mechanism))
maximal_to_minimal(mechanism, z) - x
minimal_to_maximal(mechanism, x) - z
Nb = length(mechanism.bodies)
G = attitude_jacobian(z, Nb)[1:13Nb,1:12Nb]
Array(G)
M_fd = maximal_to_minimal_jacobian(mechanism, z) * G
M_a = maximal_to_minimal_jacobian_analytical(mechanism, z)
@test size(M_fd) == size(M_a)
norm(M_fd - M_a, Inf)

ibody = 1
N_fd = minimal_to_maximal_jacobian(mechanism, x)#[(ibody-1) * 13 .+ (1:13), :]
N_a = minimal_to_maximal_jacobian_analytical(mechanism, x)#[(ibody-1) * 13 .+ (1:13), :]
sum(N_fd)
sum(N_a)
@test size(N_fd) == size(N_a)
@test norm(N_fd - N_a, Inf) < 1.0e-5

M_fd * G' * N_fd
M_a * G' * N_a
@test abs(sum(diag(M_fd * G' * N_fd)) - minimal_dimension(mechanism)) < 1.0e-6
@test abs(sum(diag(M_a * G' * N_a)) - minimal_dimension(mechanism)) < 1.0e-6

N_a * inv(N_a' * N_a)# - M_a

id = reverse(mechanism.system.dfs_list)[1]
get_node(mechanism, id) isa JointConstraint
joint = get_node(mechanism, id)
off = 0#2 * sum(input_dimension.(mechanism.joints[1:1]))

Ne = length(mechanism.joints)
nu = input_dimension(joint)
idx = collect(off .+ (1:(2nu)))

cv = minimal_coordinates_velocities(mechanism)

function joint_position_velocity(mechanism, joint, θ) 
    mechanism = deepcopy(mechanism)
    set_minimal_coordinates_velocities!(mechanism, joint, xmin=θ)
    x, v, q, ω = initial_configuration_velocity(get_body(mechanism, joint.child_id).state)
    [x; v; vector(q); ω]
end

function joint_positions_velocities(mechanism, joint, θ) 
    mechanism = deepcopy(mechanism)
    # get minimal
    cv = minimal_coordinates_velocities(mechanism)

    # set minimal 
    zp = joint_position_velocity(mechanism, joint, θ) 

    # child joints 
    child_joints = [get_node(mechanism, id) for id in recursivedirectchildren!(mechanism.system, joint.id) if get_node(mechanism, id) isa JointConstraint]

    for node in child_joints 
        zp = joint_position_velocity(mechanism, node, cv[node.id]) 
    end 
end

function joint_position_velocity_jacobian(mechanism, joint, θ)
    FiniteDiff.finite_difference_jacobian(y -> joint_position_velocity(mechanism, joint, y), θ) 
end

function joint_position_velocity(mechanism, joint, z, θ) 
    mechanism = deepcopy(mechanism)
    body_parent = get_body(mechanism, joint.parent_id)

    xp = z[1:3] 
    vp = z[4:6]
    qp = UnitQuaternion(z[7:10]..., false)
    ϕp = z[11:13]

    if body_parent.name != :origin
        set_maximal_coordinates!(body_parent, x=xp, q=qp)
        set_maximal_velocities!(body_parent, v=vp, ω=ϕp)
    end

    joint_position_velocity(mechanism, joint, θ)
end

G = zeros(maximal_dimension(mechanism), minimal_dimension(mechanism))

# root 
∂z∂θ = joint_position_velocity_jacobian(mechanism, joint, x[idx])
G[(joint.child_id - Ne - 1) * 13 .+ (1:13), idx] = ∂z∂θ
norm(N_fd[(joint.child_id - Ne - 1) * 13 .+ (1:13), idx] - ∂z∂θ)

child_ids = [id for id in recursivedirectchildren!(mechanism.system, joint.id) if get_node(mechanism, id) isa JointConstraint]
child_joints = [get_node(mechanism, id) for id in child_ids]
# child_joints = get_child_joints(mechanism, joint) 
# sort([j.id for j in child_joints])


# get_child_joints(mechanism, joint)
∂z∂z = Dict()
∂a∂z = Dict() 
push!(∂a∂z, "$(joint.child_id)" => ∂z∂θ)
# da = copy(∂z∂θ)

for node in child_joints
    haskey(∂a∂z, "$(node.child_id)") && continue
    zp = z[(node.parent_id - Ne - 1) * 13 .+ (1:13)]

    d = FiniteDiff.finite_difference_jacobian(b -> joint_position_velocity(mechanism, node, b, cv[node.id]), zp)
    push!(∂z∂z, "$(node.child_id)_$(node.parent_id)" => d)
    push!(∂a∂z, "$(node.child_id)" => d * ∂a∂z["$(node.parent_id)"])

    G[(node.child_id - Ne - 1) * 13 .+ (1:13), idx] = ∂a∂z["$(node.child_id)"]
end
# ∂z∂z["10_11"]
# ∂z∂z["11_8"]

# ∂a∂z
# ∂z∂z["10_11"]
# ∂z∂z["11_8"]

# N_fd[:, idx]
# G[:, idx]
# norm(G[:, idx] - N_fd[:, idx], Inf)

# [j.id for j in get_child_joints(mechanism, joint)] 
# [j.child_id for j in get_child_joints(mechanism, joint)]

i = joint.parent_id
i = joint.child_id
i = child_joints[1].parent_id
i = child_joints[1].child_id

i = child_joints[2].parent_id
i = child_joints[2].child_id

i = child_joints[3].parent_id
i = child_joints[3].child_id

i = child_joints[4].parent_id
i = child_joints[4].child_id

i = child_joints[5].parent_id
i = child_joints[5].child_id



i = child_joints[end].parent_id
i = child_joints[end].child_id

N_fd[(i - Ne - 1) * 13 .+ (1:13), idx]
G[(i - Ne - 1) * 13 .+ (1:13), idx]
norm(N_fd[(i - Ne - 1) * 13 .+ (1:13), idx] - G[(i - Ne - 1) * 13 .+ (1:13), idx], Inf)

norm(N_fd[:, idx] - G[:, idx], Inf)