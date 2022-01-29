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

mech = getpendulum(timestep=0.05, damper=0.0, spring=0.0);
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
data_system = create_data_system(mech.joints, mech.bodies, mech.contacts);
jacobian_data!(data_system, mech)
datajac1 = full_matrix(data_system)
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

datajac0[1:5,17:19]
datajac1[1:5,17:19]

datajac0[1:5,20:22]
datajac1[1:5,20:22]


datajac0[6:11,17:19]
datajac1[6:11,17:19]


datajac0[6:11,20:22]

datajac1[6:11,20:22]


joint0.id
body0.id


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
