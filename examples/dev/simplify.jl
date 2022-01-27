using Dojo
using MeshCat

vis = Visualizer()
open(vis)


mech = getmechanism(:humanoid, contact=true, timestep=0.05, g=-9.81, spring=30.0, damper=5.0)
initialize!(mech, :humanoid, rot=[0.1,0,0], tran=[0,0,1.5])

function ctrl!(mechanism, k)
    nu = control_dimension(mechanism)
    u = szeros(nu)
    set_control!(mechanism, u)
    return
end

storage = simulate!(mech, 2.3, ctrl!, record=true, verbose=false)
visualize(mech, storage, vis=vis)


l = 1.0
m = 1.0
joint_axis = [1.0; 0; 0]
width, depth = 0.1, 0.1
p2 = [0; 0; l/2] # joint connection point

# Links
origin = Origin{Float64}()
body1 = Box(width, depth, l, m)

# Constraints
joint_between_origin_and_body1 = JointConstraint(Revolute(origin, body1,
    joint_axis; p2=p2, spring = 0, damper = 0,
    rot_joint_limits = [SVector{1}([0.25 * π]), SVector{1}([π])]
    ))
bodies = [body1]
eqcs = [joint_between_origin_and_body1]

length(joint_between_origin_and_body1)
length(joint_between_origin_and_body1.constraints[1])
length(joint_between_origin_and_body1.constraints[2])















mech.origin.id
getfield.(mech.joints.values, :id)
getfield.(mech.bodies.values, :id)
getfield.(mech.contacts.values, :id)



full_vector(system) = vcat(getfield.(system.vector_entries,:value)...)
mech.system
eqcs = mech.joints.values
bodies = mech.bodies.values
ineqcs = mech.contacts.values
A, dims = adjacency_matrix(eqcs, bodies, ineqcs)


full_vector(system) = vcat(getfield.(system.vector_entries,:value)...)
mech.system.matrix_entries.rowval

full_matrix(mech.system)

















ineqc0
bound0 = ineqc0.constraints[1]
λ = nothing
s = SVector{4}([1.0; rand(3)/10])
γ = SVector{4}([1.0; rand(3)/10])
g(bound0, s, γ, x3, q3, v25, ϕ25)
jacv0 = constraint_jacobian_velocity(bound0, x3, q3, x2, v25, q2, ϕ25, λ, timestep)
jacv1 = FiniteDiff.finite_difference_jacobian(
    v -> g(bound0, s, γ, next_position(x2, SVector{3}(v[1:3]), timestep),
        next_orientation(q2, SVector{3}(v[4:6]), timestep), SVector{3}(v[1:3]), SVector{3}(v[4:6])), [v25; ϕ25])

jacz0 = constraint_jacobian_configuration(bound0, x3, q3, x2, v25, q2, ϕ25, λ, timestep)
jacz1 = FiniteDiff.finite_difference_jacobian(
    z -> g(bound0, s, γ, SVector{3}(z[1:3]),
        UnitQuaternion(z[4:7],false), SVector{3}(v[1:3]), SVector{3}(v[4:6])), [x3; vector(q3)])

norm(jacv0 - jacv1)
jacv0 - jacv1
norm(jacz0 - jacz1)
jacz0 - jacz1



# using Symbolics
# @variables j[1:6]
# @variables p[1:3]
# J = SMatrix{3,3,Num,9}([j[1] j[2] j[3];
#      j[2] j[4] j[5];
#      j[3] j[5] j[6]])
# Symbolics.jacobian(J * p, j)
# lift_inertia(flatten_inertia(J)) - J
# flatten_inertia(J)







@variables xa[1:3], qa[1:4], xb[1:3], qb[1:4], p[1:3], p2[1:3]
qa_ = UnitQuaternion(qa..., false)
qb_ = UnitQuaternion(qb..., false)
∂qrotation_matrix(qa_, p) * LVᵀmat(qa_)

Xt = rotation_matrix(qa_)
cbpb_w = rotation_matrix(qb_) * p2
Qt = rotation_matrix(inv(qb_)) * skew(cbpb_w) * rotation_matrix(qa_)
Symbolics.jacobian([Xt; Qt]*p, [xb; qb])




xa = rand(3)
qa = UnitQuaternion(rand(4)...)
xb = rand(3)
qb = UnitQuaternion(rand(4)...)
tra0 = eqc0.constraints[1]
tra0.vertices = (srand(3), srand(3))
Gb(tra0, xa, qa, xb, qb, 0)




using Symbolics
@variables q[1:4], p[1:3], p4[1:4]
q_ = UnitQuaternion(q..., false)
Symbolics.jacobian(VRᵀmat(q_)*p4, q)

@show Symbolics.jacobian(rotation_matrix(inv(q_))*p, q) - ∂qrotation_matrix(inv(q_), p) * Tmat()




VLmat(qa_) * RᵀVᵀmat(qa_)
inv(qb_)
Tmat() * qb

LVᵀmat(inv(qb_))
LVᵀmat(qb_)
Tmat()'

@variables p4[1:4]
Symbolics.jacobian(LVᵀmat(inv(qb_)) * p, qb) - ∂qLVᵀmat(p) * Tmat()
Symbolics.jacobian(VRᵀmat(inv(qb_)) * p4, qb)
- ∂qVRᵀmat(p4)





function ∂Ga(joint::Translational{T,Nλ,0}, λ::AbstractVector) where {T,Nλ} # ∂(G'*λ)/∂z with z = [x3, q3]
    p = constraintmat(joint)' * λ

    return
end


constraintmat(eqc0.constraints[1])
constraintmat(eqc0.constraints[2])


################################################################################
# Symbolics
################################################################################
using Symbolics

@variables xa[1:3], qa[1:4], xb[1:3], qb[1:4], p[1:3], p2[1:3]
qa_ = UnitQuaternion(qa..., false)
qb_ = UnitQuaternion(qb..., false)

Xt = rotation_matrix(qa_)
Xt - VLmat(qa_) * RᵀVᵀmat(qa_)

j = Symbolics.jacobian(Xt * p, qa)
@show j

X = transpose(rotation_matrix(qa_))
cbpb_w = rotation_matrix(qb_) * p2 # body b kinematics point
Q = transpose(rotation_matrix(inv(qb_)) * skew(cbpb_w) * rotation_matrix(qa_))
rotation_matrix(inv(qb_))

inv(qb_)

Symbolics.jacobian([X'; Q'] * p, qa)
Symbolics.jacobian([X'; Q'] * p, xb)
j = Symbolics.jacobian([X'; Q'] * p, qb)




code = build_function(j, [xa; qa; xb; qb; p; p2])[2]
fj = eval(code)
xa = rand(3)
qa = vector(UnitQuaternion(rand(4)...))
xb = rand(3)
qb = vector(UnitQuaternion(rand(4)...))
p = rand(3)
p2 = rand(3)
j = zeros(6,4)
using BenchmarkTools
@benchmark fj(j, [xa; qa; xb; qb; p; p2])

a = 10
a = 10
a = 10
a = 10
a = 10
a = 10
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################



mech = getpendulum()
mech = getdzhanibekov();
body0 = mech.bodies[1]
joint0 = mech.joints[1]

tra0 = joint0.constraints[1]
rot0 = joint0.constraints[2]

xa = rand(3)
qa = UnitQuaternion(rand(4)...)
xb = rand(3)
qb = UnitQuaternion(rand(4)...)
λ2 = rand(2)
λ3 = rand(3)

impulse_map_child(tra0, xa, qa, xb, qb, λ2)
impulse_map_parent(tra0, xa, qa, xb, qb, λ2)

impulse_map_child(rot0, xa, qa, xb, qb, λ2)
impulse_map_parent(rot0, xa, qa, xb, qb, λ2)

impulse_map_parent(mech, joint0, body0)
impulse_projector(joint0.constraints[1])
impulse_projector(joint0.constraints[2])
transpose([nullspace_mask(tra0); constraint_mask(tra0)])

impulse_map_parent(mech, joint0, body0) * szeros(0)
szeros(Float64, 6, 0) * szeros(0)

vcat(ones(6,3), zeros(6,2))
vcat(ones(6,3)', zeros(6,2)')
hcat(ones(6,3), zeros(6,2))'

@generated function impulse_map_parent(mechanism, joint::JointConstraint{T,N,Nc}, body::Body) where {T,N,Nc}
    vec = [:(impulse_map_parent(joint.constraints[$i], body, get_body(mechanism, joint.childids[$i]), joint.childids[$i], joint.λsol[2][λindex(joint,$i)], mechanism.timestep)) for i = 1:Nc]
    return :(hcat($(vec...)))
end
