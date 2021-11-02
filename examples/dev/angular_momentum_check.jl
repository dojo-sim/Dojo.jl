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
include(joinpath(module_dir(), "examples", "dev", "loader.jl"))



################################################################################
# snake
################################################################################
include("conservation_test.jl")
Δt_ = 0.01
Nlink_ = 2

function controller!(mechanism, k)
    for (i,joint) in enumerate(mechanism.eqconstraints)
        nu = getcontroldim(joint)
        if 5 >= nu >= 1
            if k ∈ (1:1)
                u = 0.0e-0 * Δt_ * [1.0; zeros(nu-1)]
            else
                u = 0.0 * [1.0; zeros(nu-1)]
            end
            setForce!(mechanism, joint, SA[u...])
        end
    end
    return
end

Random.seed!(100)
ω_ = 1.0*[1.0; 0.0; 0.0]#rand(3)
v_ = 0.0*[1.0; 1.0; 1.0]#rand(3)
Δv_ = 0.0*[0.0; 0.0; 0.0]#rand(3)
Δω_ = 1.0*[1.0; 1.0; 1.0]#rand(3)
ϕ1_ = 0.0
# jointtype = :Fixed
# mech = getmechanism(:snake, Δt = Δt_, g = 0.00, contact = false, spring = 0.0, damper = 0.0, Nlink = Nlink_, jointtype = jointtype)
# initialize!(mech, :snake, ϕ1 = ϕ1_, v=v_, ω=ω_, Δv = Δv_, Δω = Δω_)
mech = getmechanism(:npendulum, Δt = 0.01, g = 0.0 * -9.81, Nlink = Nlink_)
initialize!(mech, :npendulum, ϕ1 = 0.0 * π, v=v_, ω=ω_, Δv = Δv_, Δω = Δω_)

storage = simulate!(mech, 1.00, controller!, record = true, solver = :mehrotra!, verbose = false)
m0 = momentum(mech)
e0 = mechanicalEnergy(mech)

# mech = getmechanism(:snake, Δt = Δt_, g = 0.00, contact = false, spring = 0.0, damper = 0.0, Nlink = Nlink_, jointtype = jointtype)
# initialize!(mech, :snake, ϕ1 = ϕ1_, v=v_, ω=ω_, Δv = Δv_, Δω = Δω_)
mech = getmechanism(:npendulum, Δt = 0.01, g = 0.0 * -9.81, Nlink = Nlink_)
initialize!(mech, :npendulum, ϕ1 = 0.0 * π, v=v_, ω=ω_, Δv = Δv_, Δω = Δω_)

storage = simulate!(mech, 10.00, controller!, record = true, solver = :mehrotra!, verbose = false)
m1 = momentum(mech)
e1 = mechanicalEnergy(mech)

norm((m1 - m0)[1:3], Inf)
norm((m1 - m0)[4:6], Inf)
norm(m1[1:3])
norm(m0[1:3])
norm(m1[4:6])
norm(m0[4:6])

m1[4:6]
m0[4:6]
visualize(mech, storage, vis = vis)

plot(hcat(Vector.(storage.x[1])...)')
plot!(hcat(Vector.(storage.x[2])...)')

plot(hcat([[q.w, q.x, q.y, q.z] for q in storage.q[1]]...)', width=2.0, color=:black)
plot!(hcat([[q.w, q.x, q.y, q.z] for q in storage.q[2]]...)', width=1.0, color=:red)

plot(hcat(Vector.(storage.ω[1])...)', width=2.0, color=:black, label="")
plot!(hcat(Vector.(storage.ω[2])...)', width=1.0, color=:red, label="")
