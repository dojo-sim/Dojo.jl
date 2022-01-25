using Dojo
using MeshCat

vis = Visualizer()
open(vis)


mech = getmechanism(:humanoid, contact=true, Δt=0.05, g=-9.81, spring=30.0, damper=5.0)
initialize!(mech, :humanoid, rot=[0.1,0,0], tran=[0,0,1.5])

function ctrl!(mechanism, k)
    nu = controldim(mechanism)
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
getfield.(mech.eqconstraints.values, :id)
getfield.(mech.bodies.values, :id)
getfield.(mech.ineqconstraints.values, :id)



full_vector(system) = vcat(getfield.(system.vector_entries,:value)...)
mech.system
eqcs = mech.eqconstraints.values
bodies = mech.bodies.values
ineqcs = mech.ineqconstraints.values
A, dims = adjacencyMatrix(eqcs, bodies, ineqcs)


full_vector(system) = vcat(getfield.(system.vector_entries,:value)...)
mech.system.matrix_entries.rowval

full_matrix(mech.system)

















ineqc0
bound0 = ineqc0.constraints[1]
λ = nothing
s = SVector{4}([1.0; rand(3)/10])
γ = SVector{4}([1.0; rand(3)/10])
g(bound0, s, γ, x3, q3, v25, ϕ25)
jacv0 = ∂g∂v(bound0, x3, q3, x2, v25, q2, ϕ25, λ, Δt)
jacv1 = FiniteDiff.finite_difference_jacobian(
    v -> g(bound0, s, γ, getx3(x2, SVector{3}(v[1:3]), Δt),
        getq3(q2, SVector{3}(v[4:6]), Δt), SVector{3}(v[1:3]), SVector{3}(v[4:6])), [v25; ϕ25])

jacz0 = ∂g∂z(bound0, x3, q3, x2, v25, q2, ϕ25, λ, Δt)
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
