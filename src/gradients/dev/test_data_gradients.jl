using Dojo
using MeshCat

vis = Visualizer()
open(vis)

################################################################################
# Analytical Jacobian
################################################################################

function test_data_system(model::Symbol; ϵ::T=1.0e-6, tsim::T=0.1, ctrl::Any=(m,k)->nothing,
        Δt::T=0.01, g::T=-9.81, verbose::Bool=false, kwargs...) where T

    # mechanism
    mechanism = getmechanism(model, Δt=Δt, g=g; kwargs...)
    initialize!(mechanism, model)

    # simulate
    storage = simulate!(mechanism, tsim, ctrl,
        record=true, verbose=false, opts=InteriorPointOptions(rtol=ϵ, btol=ϵ))

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

mech = getpendulum(Δt=0.05, damper=0.0, spring=0.0);
eqc0 = mech.eqconstraints.values[1]
body0 = mech.bodies.values[1]
initialize!(mech, :pendulum, ϕ1=0.2, ω1=-0.3)
simulate!(mech, 0.30, verbose=false)

# Finite Difference
Nd = data_dim(mech, attjac=false)
data0 = get_data(mech)# + 0.05*rand(Nd)
sol0 = get_solution(mech)
datajac0 = finitediff_data_jacobian(mech, data0, sol0)
attjac0 = data_attitude_jacobian(mech)
datajac0 *= attjac0
plot(Gray.(1e10*abs.(datajac0)))
plot(Gray.(1e0*abs.(datajac0)))

# Analytical
data_system = create_data_system(mech.eqconstraints.values,
    mech.bodies.values, mech.ineqconstraints.values);
∂data!(data_system, mech)
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

datajac0[6:11,17:19]
datajac1[6:11,17:19]


datajac0[6:11,20:22]

datajac1[6:11,20:22]


eqc0.id
body0.id


data_system.matrix_entries[eqc0.id, eqc0.id].value
data_system.matrix_entries[eqc0.id, body0.id].value
data_system.matrix_entries[body0.id, eqc0.id].value
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
ineqc0 = mech.ineqconstraints.values[1]
eqc0 = mech.eqconstraints.values[2]
body0 = mech.bodies.values[1]

∇0 = ∂eqc∂eqc_data(mech, eqc0)
∇0 = ∂eqc∂body_data(mech, eqc0, body0)
∇0 = ∂ineqc∂body_data(mech, ineqc0, body0)
∇0 = ∂ineqc∂ineqc_data(mech, ineqc0, body0)
∇0 = ∂body∂eqc_data(mech, eqc0, body0)
∇0 = ∂body∂body_data(mech, body0)
∇0 = ∂body∂ineqc_data(mech, ineqc0, body0)

data_system = create_data_system(mech.eqconstraints.values,
    mech.bodies.values, mech.ineqconstraints.values);

∂ineqc_data!(data_system, mech)
plot(Gray.(1e10 .* abs.(full_matrix(data_system))))

∂body_data!(data_system, mech)
plot(Gray.(1e10 .* abs.(full_matrix(data_system))))

∂eqc_data!(data_system, mech)
plot(Gray.(1e10 .* abs.(full_matrix(data_system))))
plot(log.(10, abs.(sum(full_matrix(data_system), dims=1)[1,:])))


full_matrix(data_system)
