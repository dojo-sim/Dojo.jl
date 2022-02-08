using Dojo
using MeshCat

vis = Visualizer()
open(vis)

################################################################################
# Analytical Jacobian
################################################################################

function test_data_system(model::Symbol; ϵ::T=1.0e-6, tsim::T=0.1, ctrl::Any=(m,k)->nothing,
        timestep::T=0.01, gravity=[0.0; 0.0; -9.81], verbose::Bool=false, kwargs...) where T

    # mechanism
    mechanism = get_mechanism(model, timestep=timestep, gravity=gravity; kwargs...)
    initialize!(mechanism, model)

    # simulate
    storage = simulate!(mechanism, tsim, ctrl,
        record=true, verbose=false, opts=SolverOptions(rtol=ϵ, btol=ϵ))

    # Set data
    Nb = data_dim(mechanism)
    data = getdata(mechanism)
    setdata!(mechanism, data)
    sol = getsolution(mechanism)

    return nothing
end

# test_data_system(:snake, Nb=3)
using Test

include("utils.jl")
include("data.jl")
include("data_gradients.jl")
include("finite_difference.jl")

mech = get_pendulum(timestep=0.05, damper=3.0, spring=1.0);
joint0 = mech.joints[1]
body0 = mech.bodies[1]
initialize!(mech, :pendulum, ϕ1=0.2, ω1=-0.3)
simulate!(mech, 0.30, verbose=false)

# Finite Difference
Nd = data_dim(mech, attjac=false)
data0 = get_data0(mech)# + 0.05*rand(Nd)
sol0 = get_solution0(mech)
datajac0 = finitediff_data_jacobian(mech, data0, sol0)
attjac0 = data_attitude_jacobian(mech)
datajac0 *= attjac0
plot(Gray.(1e10*abs.(datajac0)))
plot(Gray.(1e0*abs.(datajac0)))


# Analytical
D = create_data_matrix(mech.joints, mech.bodies, mech.contacts)
jacobian_data!(D, mech)

nodes = [mech.joints; mech.bodies; mech.contacts]
dimrow = length.(nodes)
dimcol = data_dim.(nodes)
datajac1 = full_matrix(D, dimrow, dimcol)
plot(Gray.(1e10 .* abs.(datajac1)))
plot(Gray.(1e0 .* abs.(datajac1)))

plot(Gray.(1e10 .* abs.(datajac0)))
plot(Gray.(1e10 .* abs.(datajac1)))
plot(Gray.(1e6 .* abs.(datajac0 - datajac1)))
plot(Gray.(1e0 .* abs.(datajac0 - datajac1)))

norm((datajac0 - datajac1)[1:5,1:3])
norm((datajac0 - datajac1)[1:5,4:4])
norm((datajac0 - datajac1)[1:5,5:10])
norm((datajac0 - datajac1)[1:5,11:16])
norm((datajac0 - datajac1)[1:5,17:19])
norm((datajac0 - datajac1)[1:5,20:22])
norm((datajac0 - datajac1)[6:11,1:1])
norm((datajac0 - datajac1)[6:11,2:2])
norm((datajac0 - datajac1)[6:11,3:3])
norm((datajac0 - datajac1)[6:11,4:4])
norm((datajac0 - datajac1)[6:11,5:10])
norm((datajac0 - datajac1)[6:11,11:13])
norm((datajac0 - datajac1)[6:11,14:16])
norm((datajac0 - datajac1)[6:11,17:19], Inf)
norm((datajac0 - datajac1)[6:11,20:22], Inf)

datajac0[6:11,17:19]
datajac1[6:11,17:19]

datajac0[6:11,20:22]
datajac1[6:11,20:22]



λ10 = srand(length(joint0.constraints[1]))
λ20 = srand(length(joint0.constraints[2]))
impulse_map_child_jacobian_child(joint0.constraints[1], mech.origin, body0, λ10)
impulse_map_child_jacobian_child(joint0.constraints[2], mech.origin, body0, λ20)



datajac0[6:11,17:19]

datajac1[6:11,17:19]




datajac0[6:11,20:22]

datajac1[6:11,20:22]



data_system.matrix_entries[joint0.id, joint0.id].value
data_system.matrix_entries[joint0.id, body0.id].value
data_system.matrix_entries[body0.id, joint0.id].value
data_system.matrix_entries[body0.id, body0.id].value






data_system = create_data_system(mech.eqconstraints.values,
    mech.bodies.values, mech.ineqconstraints.values);
plot(Gray.(1e10 .* abs.(full_matrix(data_system))))
∂eqc_data!(data_system, mech)
plot(Gray.(1e10 .* abs.(full_matrix(data_system)[:,1:4])))
plot(Gray.(1e10 .* abs.(datajac0[:,1:4])))



datajac1 = full_matrix(data_system)
plot(Gray.(1e10 .* abs.(datajac1)))



#
# ∂eqc_data!(data_system, mech)
# plot(Gray.(1e10 .* abs.(full_matrix(data_system))))
#
# ∂body_data!(data_system, mech)
# plot(Gray.(1e10 .* abs.(full_matrix(data_system))))
#
# ∂ineqc_data!(data_system, mech)
# plot(Gray.(1e10 .* abs.(full_matrix(data_system))))
# plot(log.(10, abs.(sum(full_matrix(data_system), dims=1)[1,:])))




initialize!(mech, :snake, x=[0,0,1.0])
storage = simulate!(mech, 1.35, ctrl!, record=true, verbose=false)
visualize(mech, storage, vis=vis)
contact0 = mech.contacts.values[1]
joint0 = mech.joints.values[2]
body0 = mech.bodies.values[1]

∇0 = ∂joint∂joint_data(mech, joint0)
∇0 = ∂joint∂body_data(mech, joint0, body0)
∇0 = ∂contact∂body_data(mech, contact0, body0)
∇0 = ∂contact∂contact_data(mech, contact0, body0)
∇0 = ∂body∂joint_data(mech, joint0, body0)
∇0 = ∂body∂body_data(mech, body0)
∇0 = ∂body∂contact_data(mech, contact0, body0)

data_system = create_data_system(mech.joints.values,
    mech.bodies.values, mech.contacts.values);

∂contact_data!(data_system, mech)
plot(Gray.(1e10 .* abs.(full_matrix(data_system))))

∂body_data!(data_system, mech)
plot(Gray.(1e10 .* abs.(full_matrix(data_system))))

∂joint_data!(data_system, mech)
plot(Gray.(1e10 .* abs.(full_matrix(data_system))))
plot(log.(10, abs.(sum(full_matrix(data_system), dims=1)[1,:])))




full_matrix(data_system)


mech = get_pendulum()
mech = get_halfcheetah()
A = data_adjacency_matrix(mech.joints, mech.bodies, mech.contacts)
sum(A)
plot(Gray.(A))
typeof(mech.joints[1]) <: JointConstraint
typeof(mech.joints[1]) <: JointConstraint

A = data_adjacency_matrix(mech.joints, mech.bodies, mech.contacts)
D = create_data_matrix(mech.joints, mech.bodies, mech.contacts)



# mech = get_snake(jointtype=:PlanarAxis, spring=2.0, damper=4.0)
mech = get_snake(jointtype=:Prismatic, spring=0.0, damper=4.0, gravity=0.0)
initialize!(mech, :snake, v=[0,0,0.], ω=[0,0,0.])
function ctrl!(mech, k)
    nu = control_dimension(mech)
    set_control!(mech, [szeros(6); 0.01sones(nu-6)]*mech.timestep)
end
storage = simulate!(mech, 5.3, ctrl!, record=true, verbose=false)
visualize(mech, storage, vis=vis)


joint1 = mech.joints[1]
joint2 = mech.joints[2]
body1 = mech.bodies[1]
body2 = mech.bodies[2]
x0, v0, q0, ω0 = current_configuration_velocity(mech.origin.state)
x1, v1, q1, ω1 = current_configuration_velocity(body1.state)
x2, v2, q2, ω2 = current_configuration_velocity(body2.state)



Fτ1 = SVector{3}(-1,-2,3.0)
Fτ2 = SVector{3}(1,2,3.0)

position_error(joint1.constraints[1], x0, q0, x1, q1)
position_error(joint2.constraints[1], x1, q1, x2, q2)
spring_parent(joint1.constraints[1], x0, q0, x1, q1, unitary=true)
spring_parent(joint2.constraints[1], x1, q1, x2, q2)

spring_child(joint1.constraints[1], x0, q0, x1, q1, unitary=true)
spring_child(joint2.constraints[1], x1, q1, x2, q2)

minimal_velocities(joint1.constraints[1], x0, q0, v0, ω0, x1, q1, v1, ω1)
minimal_velocities(joint2.constraints[1], x0, q0, v0, ω0, x1, q1, v1, ω1)

damper_force(joint1.constraints[1], x0, q0, v0, ω0, x1, q1, v1, ω1)
damper_parent(joint1.constraints[1], x0, q0, v0, ω0, x1, q1, v1, ω1, unitary=true)
damper_parent(joint2.constraints[1], x1, q1, v1, ω1, x2, q2, v2, ω2)

damper_child(joint1.constraints[1], x0, q0, v0, ω0, x1, q1, v1, ω1, unitary=true)
damper_child(joint2.constraints[1], x1, q1, v1, ω1, x2, q2, v2, ω2)





apply_input(joint1.constraints[1], Fτ1, x0, q0, x1, q1)
apply_input(joint2.constraints[1], Fτ2, x1, q1, x2, q2)

input_jacobian_control_parent(joint1.constraints[1], mech.origin.state, body2.state)
input_jacobian_control_parent(joint2.constraints[1], body1.state, body2.state)

input_jacobian_control_child(joint1.constraints[1], mech.origin.state, body2.state)
input_jacobian_control_child(joint2.constraints[1], body1.state, body2.state)

input_jacobian_configuration_parent(joint1.constraints[1], mech.origin.state, body2.state)
input_jacobian_configuration_parent(joint2.constraints[1], body1.state, body2.state)

input_jacobian_configuration_child(joint1.constraints[1], mech.origin.state, body2.state)
input_jacobian_configuration_child(joint2.constraints[1], body1.state, body2.state)


apply_input(joint1.constraints[2], Fτ1, x0, q0, x1, q1)
apply_input(joint2.constraints[2], Fτ2, x1, q1, x2, q2)
