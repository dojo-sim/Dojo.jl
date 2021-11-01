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
            if k ∈ (1:40)
                u = 3e-0 * Δt_ * [1.0; zeros(nu-1)]
            else
                u = 0.0 * [1.0; zeros(nu-1)]
            end
            setForce!(mechanism, joint, SA[u...])
        end
    end
    return
end

Random.seed!(300)
ω_ = 1.0*rand(3)
v_ = 1.0*rand(3)
Δv_ = 1.0*rand(3)
Δω_ = 1.0*rand(3)
ϕ1_ = pi/2
jointtype = :Prismatic
mech = getmechanism(:snake, Δt = Δt_, g = 0.00, contact = false, spring = 0.0, damper = 0.3, Nlink = Nlink_, jointtype = jointtype)
initialize!(mech, :snake, ϕ1 = ϕ1_, Δv = Δv_, Δω = Δω_)
storage = simulate!(mech, 1.00, controller!, record = true, solver = :mehrotra!, verbose = false)
m0 = momentum(mech)

mech = getmechanism(:snake, Δt = Δt_, g = 0.00, contact = false, spring = 0.0, damper = 0.3, Nlink = Nlink_, jointtype = jointtype)
initialize!(mech, :snake, ϕ1 = ϕ1_, Δv = Δv_, Δω = Δω_)
storage = simulate!(mech, 2.00, controller!, record = true, solver = :mehrotra!, verbose = false)
m1 = momentum(mech)
norm((m1 - m0)[4:6], Inf)
(m1 - m0)[5]
m1 - m0
visualize(mech, storage, vis = vis)







################################################################################
# humanoid
################################################################################

include("conservation_test.jl")
Random.seed!(100)
Δt_ = 0.01
function controller!(mechanism, k)
    for (i,joint) in enumerate(mechanism.eqconstraints)
        if getcontroldim(joint) == 1
            if k ∈ (1:100)
                u = 1e0 * Δt_
            else
                u = 0.0
            end
            setForce!(mechanism, joint, SA[u])
        end
    end
    return
end

mech = getmechanism(:humanoid, Δt = Δt_, g = 0.00, contact = false, spring = 0.0, damper = 1.0)
initialize!(mech, :humanoid)
storage = simulate!(mech, 2.00, controller!, record = true, solver = :mehrotra!, verbose = false)
m0 = momentum(mech)

mech = getmechanism(:humanoid, Δt = Δt_, g = 0.00, contact = false, spring = 0.0, damper = 1.0)
initialize!(mech, :humanoid)
storage = simulate!(mech, 60.00, controller!, record = true, solver = :mehrotra!, verbose = false)
m1 = momentum(mech)
norm((m1 - m0)[4:6], Inf)
(m1 - m0)[5]
m1 - m0
visualize(mech, storage, vis = vis)
