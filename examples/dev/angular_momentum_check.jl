# Utils
function module_dir()
    return joinpath(@__DIR__, "..", "..")
end

# Activate package
using Pkg
Pkg.activate(module_dir())

# Load packages
using Plots
using Random
using MeshCat

# Open visualizer
vis = Visualizer()
open(vis)

# Include new files
include(joinpath(module_dir(), "examples", "loader.jl"))

################################################################################
# snake
################################################################################
include("conservation_test.jl")
timestep_ = 0.01
Nb_ = 2

function controller!(mechanism, k)
    for (i,joint) in enumerate(mechanism.joints)
        nu = input_dimension(joint)
        if 5 >= nu >= 1
            if k ∈ (1:1)
                u = 0.0e-0 * timestep_ * [1.0; zeros(nu-1)]
            else
                u = 0.0 * [1.0; zeros(nu-1)]
            end
            set_input!(joint, SA[u...])
        end
    end
    return
end

Random.seed!(100)
ω_ = 1.0 *[1.0; 0.0; 0.0]#rand(3)
v_ = 0.0*[1.0; 0.0; 0.0]#rand(3)
Δv_ = 0.0*[0.0; 0.0; 0.0]#rand(3)
Δω_ = 1.0*[0.0; 1.0; 0.0]#rand(3)
ϕ1_ = 0.0
jointtype = :Spherical
mech = getmechanism(:snake, timestep = timestep_, g = 0.0, contact = false, spring = 0.0, damper = 0.0, Nb = Nb_, jointtype = jointtype)
initialize!(mech, :snake, ϕ1 = ϕ1_, v=v_, ω=ω_, Δv = Δv_, Δω = Δω_)
# mech = getmechanism(:npendulum, timestep = 0.01, g = 0.0 * -9.81, Nb = Nb_)
# initialize!(mech, :npendulum, ϕ1 = 0.0 * π, v=v_, ω=ω_, Δv = Δv_, Δω = Δω_)

storage = simulate!(mech, 1.0, record = true, solver = :mehrotra!, verbose = false)
m0 = momentum(mech)
e0 = mechanical_energy(mech)
mech = getmechanism(:snake, timestep = timestep_, g = 0.00, contact = false, spring = 0.0, damper = 0.0, Nb = Nb_, jointtype = jointtype)
initialize!(mech, :snake, ϕ1 = ϕ1_, v=v_, ω=ω_, Δv = Δv_, Δω = Δω_)
# mech = getmechanism(:npendulum, timestep = 0.01, g = 0.0 * -9.81, Nb = Nb_)
# initialize!(mech, :npendulum, ϕ1 = 0.0 * π, v=v_, ω=ω_, Δv = Δv_, Δω = Δω_)

storage = simulate!(mech, 5.00, record = true, solver = :mehrotra!, verbose = false)
m1 = momentum(mech)
e1 = mechanical_energy(mech)

abs(e1 - e0)
norm((m1 - m0)[1:3], Inf)
norm((m1 - m0)[4:6], Inf)
norm(m1[1:3])
norm(m0[1:3])
norm(m1[4:6])
norm(m0[4:6])

eqc = mech.joints[2]
f1 = (zerodimstaticadjoint(∂g∂ʳpos(mech, eqc, mech.bodies[3])) * eqc.λsol[2])[1:3]
f2 = (zerodimstaticadjoint(∂g∂ʳpos(mech, eqc, mech.bodies[4])) * eqc.λsol[2])[1:3]
norm(f1 + f2)
norm(f1) - norm(f2)
t1 = (zerodimstaticadjoint(∂g∂ʳpos(mech, eqc, mech.bodies[3])) * eqc.λsol[2])[4:6]
t2 = (zerodimstaticadjoint(∂g∂ʳpos(mech, eqc, mech.bodies[4])) * eqc.λsol[2])[4:6]
norm(t1 + t2)
norm(t1) - norm(t2)

zerodimstaticadjoint(∂g∂ʳpos(mech, eqc, mech.bodies[3])) * eqc.λsol[2]


∂g∂ʳpos(mech, eqc, mech.bodies[3])
∂g∂ʳpos(mech, eqc, mech.bodies[4])
mech.bodies[3].state.q2
mech.bodies[4].state.q2

visualize(mech, storage, vis = vis)

plot(hcat(Vector.(storage.x[1])...)')
plot!(hcat(Vector.(storage.x[2])...)')

plot(hcat([[q.w, q.x, q.y, q.z] for q in storage.q[1]]...)', width=2.0, color=:black)
plot!(hcat([[q.w, q.x, q.y, q.z] for q in storage.q[2]]...)', width=1.0, color=:red)

plot(hcat(Vector.(storage.ω[1])...)', width=2.0, color=:black, label="")
plot!(hcat(Vector.(storage.ω[2])...)', width=1.0, color=:red, label="")


# test solmat
data = get_data(mech)
set_data!(mech, data)
sol = get_solution(mech)
Nb = length(collect(mech.bodies))
attjac = attitude_jacobian(data, Nb)

# IFT
set_entries!(mech)
datamat = full_data_matrix(deepcopy(mech))
solmat = full_matrix(mech.system)
sensi = - (solmat \ datamat)

# finite diff
fd_datamat = finitediff_data_matrix(deepcopy(mech), data, sol, δ = 1e-5) * attjac
@test norm(fd_datamat + datamat, Inf) < 1e-8

fd_solmat = finitediff_sol_matrix(mech, data, sol, δ = 1e-5)
@test norm(fd_solmat + solmat, Inf) < 1e-8

###
################################################################################
# snake initial velocity
################################################################################
include("conservation_test.jl")
n = 1
timestep_ = 0.01
timestep_ /= n
Nb_ = 2

Random.seed!(100)
ω_ = 1.0 * -1.0*2.0*[0,0,1.0] #* 100 * timestep_
v_ = 0.0*rand(3)
Δv_ = 0.0*rand(3)
Δω_ = 1.0 * 1.0*4.0*[0,0,1.0] #* 100 * timestep_

function controller!(mechanism, k)
    for (i,joint) in enumerate(mechanism.joints)
        nu = input_dimension(joint)
        if nu <= 5
            if k ∈ (10:10 + 100n)
                u = 1.0 * 3e-2 * timestep_ * [1.0, 0.0, 0.0] #[0.0; 1.0; zeros(nu-2)]
            else
                u = zeros(nu)
            end
            set_input!(joint, SA[u...])
        end
    end
    return
end
mech = getmechanism(:snake, timestep = timestep_, g = 0.0, spring = 0.0, damper = 0.05, contact = false, Nb = Nb_, jointtype = :Spherical)
initialize!(mech, :snake, ω = ω_, v = v_, Δv = Δv_, Δω = Δω_)
storage = simulate!(mech, 200.0, controller!, record = true, solver = :mehrotra!, verbose = false)
m0 = momentum(mech)


visualize(mech, storage, vis = vis)
