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
Δt_ = 0.02
Nlink_ = 2

Random.seed!(50)
ω_ = 0.0*rand(3)
v_ = 0.0*rand(3)
Δv_ = 0.0*rand(3)
Δω_ = 11.0*rand(3)
ϕ1_ = pi/2

function controller!(mechanism, k)
    for (i,joint) in enumerate(mechanism.eqconstraints)
        if getcontroldim(joint) == 1
            if k ∈ (1:50)
                u = 1e0 * Δt_
            else
                u = 0.0
            end
            setForce!(mechanism, joint, SA[u])
        end
    end
    return
end

function get_momentum(h)
    mech = getmechanism(:snake, Δt = Δt_, g = 0.00, contact = false, spring = 0.0,
        damper = 0.0, Nlink = Nlink_, jointtype = jointtype)
    initialize!(mech, :snake, ϕ1 = ϕ1_, v = v_, ω = ω_, Δv = Δv_, Δω = Δω_)
    storage = simulate!(mech, h, controller!, record = true, solver = :mehrotra!, verbose = false)
    m0 = momentum(mech)
    return norm(m0[4:6])
end

plot([get_momentum(1 + 0.05i) for i = 1:50])

jointtype = :Revolute
mech = getmechanism(:snake, Δt = Δt_, g = 0.00, contact = false, spring = 0.0,
    damper = 0.0, Nlink = Nlink_, jointtype = jointtype)
initialize!(mech, :snake, ϕ1 = ϕ1_, v = v_, ω = ω_, Δv = Δv_, Δω = Δω_)
storage = simulate!(mech, 1.00, controller!, record = true, solver = :mehrotra!, verbose = false)
m0 = momentum(mech)

mech = getmechanism(:snake, Δt = Δt_, g = 0.00, contact = false, spring = 0.0,
    damper = 0.0, Nlink = Nlink_, jointtype = jointtype)
initialize!(mech, :snake, ϕ1 = ϕ1_, v = v_, ω = ω_, Δv = Δv_, Δω = Δω_)
storage = simulate!(mech, 2.00, controller!, record = true, solver = :mehrotra!, verbose = false)
m1 = momentum(mech)
norm((m1 - m0)[4:6], Inf)
(m1 - m0)[5]
m1 - m0
norm(m0[1:3]) - norm(m1[1:3])
norm(m0[4:6]) - norm(m1[4:6])

visualize(mech, storage, vis = vis)

plot(hcat([Vector(x) for x in storage.ω[1]]...)')

ϕ = 0.05
q = UnitQuaternion(RotZ(ϕ))
v = [1, 0, 0]
vrotate(v, q)


joint1.vertices



Δt = 0.01
λ = [0, 0, 0, 1, 1.]
eqc2 = collect(mech.eqconstraints)[2]
joint1 = eqc2.constraints[1]
joint2 = eqc2.constraints[2]
body1 = collect(mech.bodies)[1]
body2 = collect(mech.bodies)[2]
Xa, Qa = ∂g∂posa(joint2, posargsk(body1.state)..., posargsk(body2.state)...)
Xb, Qb = ∂g∂posb(joint2, posargsk(body1.state)..., posargsk(body2.state)...)

_, qa = posargsk(body1.state)
_, qb = posargsk(body2.state)

norm([Xa Qa] + [Xb Qb])
norm([Xa Qa * LVᵀmat(qa)] + [Xb Qb * LVᵀmat(qb)])


λa = zerodimstaticadjoint(∂g∂ʳpos(mech, eqc2, body1)) * λ
λb = zerodimstaticadjoint(∂g∂ʳpos(mech, eqc2, body2)) * λ
λa_ = zerodimstaticadjoint(constraintmat(joint2) * [Xa Qa * LVᵀmat(qa)]) * λ[4:5]
λb_ = zerodimstaticadjoint(constraintmat(joint2) * [Xb Qb * LVᵀmat(qb)]) * λ[4:5]

norm(∂g∂ʳposa(joint2, body1, body2, body2.id, Δt) - constraintmat(joint2) * [Xa Qa * LVᵀmat(qa)])
norm(∂g∂ʳposb(joint2, body1, body2, body2.id, Δt) - constraintmat(joint2) * [Xb Qb * LVᵀmat(qb)])

vrotate(λa[4:6], qa) + vrotate(λb[4:6], qb)


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




function get_momentum(h)
    mech = getmechanism(:humanoid, Δt = Δt_, g = 0.00, contact = false, spring = 0.0, damper = 1.0)
    initialize!(mech, :humanoid)
    storage = simulate!(mech, h, controller!, record = true, solver = :mehrotra!, verbose = false)
    m0 = momentum(mech)
    return norm(m0[4:6])
end

plot([get_momentum(2 + 0.05i) for i = 1:10])
