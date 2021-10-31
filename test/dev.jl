
vis = Visualizer()
open(vis)

# mech = getmechanism(:snake, contact = false, spring = 000.0)
# initialize!(mech, :snake)
# storage = simulate!(mech, 10.1, record = true, solver = :mehrotra!, verbose = false)
# visualize(mech, storage, vis = vis)


# body1 = collect(mech.bodies)[1]
# eqc1 = collect(mech.eqconstraints)[1]
# eqc2 = collect(mech.eqconstraints)[2]
# eqc1.parentid
# eqc2.parentid
# eqc2.constraints[1].spring
include("conservation_test.jl")

potentialEnergy(mech)
kineticEnergy(mech)
mechanicalEnergy(mech)

function get_energy(tsim, Δt)
    # mech = getmechanism(:slider, Δt = Δt, g = -1.0, spring = 100.0)
    # initialize!(mech, :slider, z1 = 2.0)
    mech = getmechanism(:nslider, Δt = Δt, g = -9.81, spring = 10.0, Nlink = 5)
    initialize!(mech, :nslider, Δz = 1.0)
    telap = @elapsed simulate!(mech, tsim, record = false, solver = :mehrotra!, verbose = false)

    ET = kineticEnergy(mech)
    EV = potentialEnergy(mech)
    EM = mechanicalEnergy(mech)
    return [ET, EV, EM, tsim / telap]
end

ts = [0.11 * i for i = 1:50]
energies = [get_energy(t, 0.02) for t in ts]
plt = plot()
plot!(ts, hcat(energies...)[1,:], linewidth = 3.0, label = "kinetic")
plot!(ts, hcat(energies...)[2,:], linewidth = 3.0, label = "potential")
plot!(ts, hcat(energies...)[3,:], linewidth = 3.0, label = "mechanical")
plt = plot()
energies[1][3]
err = abs.((hcat(energies...)[3,:] .- energies[1][3]) / energies[1][3])
plot(ts, log.(10, err), linewidth = 3.0, label = "mechanical", ylims = (-8, 2))


Δts = [1e-6*exp(1/4 * log(10))^(i-1) for i = 1:20]
energies = [get_energy(Δt*897, Δt) for Δt in Δts]
plt = plot()
err = abs.((hcat(energies...)[3,:] .- energies[10][3]) / energies[10][3])
tratio = hcat(energies...)[4,:]
scatter(log.(10, tratio), log.(10, err), linewidth = 3.0, label = "mechanical", ylims = (-10, 2))

mech = getmechanism(:pendulum)
initialize!(mech, :pendulum, ϕ1 = π/2)
storage = simulate!(mech, 2.0, record = true, solver = :mehrotra!, verbose = false)
visualize(mech, storage, vis = vis)


mech = getmechanism(:nslider, Δt = 0.01, g = -10.0, spring = 100.0, Nlink = 5)
initialize!(mech, :nslider, z1 = 0.0, Δz = 0.5)
storage = simulate!(mech, 2.5, record = true, solver = :mehrotra!, verbose = false)
visualize(mech, storage, vis = vis)


body1 = collect(mech.bodies)[1]
body1.state
eqc1 = collect(mech.eqconstraints)[1]
tra1 = eqc1.constraints[1]
tra1.spring
potentialEnergy(mech, eqc1)
0.5 * 100 * z*2

gv = mech.g
m = body1.m
x, q = posargsk(body1.state)
z = x[3]
EV = potentialEnergy(mech, body1)
potentialEnergy(mech, body1) + m * gv * z


################################################################################
# atlas
################################################################################

include("conservation_test.jl")
Random.seed!(100)
Δt_ = 0.01
mech = getmechanism(:atlas, Δt = Δt_, g = 0.0, spring = 0.0, damper = 0.05, contact = false)
initialize!(mech, :atlas, tran = [0, 0, 2.])

function controller!(mechanism, k)
    for (i,joint) in enumerate(mechanism.eqconstraints)
        if getcontroldim(joint) == 1
            if k ∈ (1:50)
                u = 1e0 * (1.0 - 0.5) * Δt_
            elseif k ∈ (51:100)
                u = -1e0 * (1.0 - 0.5) * Δt_
            else
                u = 0.0
            end
            setForce!(mechanism, joint, SA[u])
        end
    end
    return
end

storage = simulate!(mech, 1.0, controller!, record = true, solver = :mehrotra!, verbose = false)
visualize(mech, storage, vis = vis)

# mech.eqconstraints[1].constraints[1].spring
# mech.eqconstraints[1].constraints[2].spring
# mech.eqconstraints[1].constraints[1].damper
# mech.eqconstraints[1].constraints[2].damper
momentum(mech)

[eqc.name for eqc in mech.eqconstraints]
[(body.name, body.m) for body in mech.bodies]

################################################################################
# snake
################################################################################
include("conservation_test.jl")
n = 1000
Δt_ = 0.1 / n
Nlink_ = 2
mech = getmechanism(:snake, Δt = Δt_, g = 0.0, spring = 0.0, damper = 0.0, contact = false, Nlink = Nlink_)
initialize!(mech, :snake, Δv = zeros(3), Δω = zeros(3))

function controller!(mechanism, k)
    for (i,joint) in enumerate(collect(mechanism.eqconstraints)[2:end])
        if  i == Nlink_-1
            if getcontroldim(joint) == 1
                # u = 5e-1 * (rand() - 0.5) * Δt_ * (k < 100)
                if k ∈ (1:n)
                    u = 1e-0 * Δt_
                elseif k ∈ (n+1:2n)
                    u = -1e-0 * Δt_
                else
                    u = 0.0
                end
                setForce!(mechanism, joint, SA[u])
            end
        end
    end
    return
end

storage = simulate!(mech, 1.0, controller!, record = true, solver = :mehrotra!, verbose = true)
visualize(mech, storage, vis = vis)

momentum(mech, collect(mech.bodies)[1])[4]
momentum(mech, collect(mech.bodies)[2])[4]
momentum(mech)[4]

plot([ω[1] for ω in storage.ω[1]])
plot([ω[1] for ω in storage.ω[2]])

collect(mech.bodies)[1].m

bodies = collect(mech.bodies)
eqcs = collect(mech.eqconstraints)
body1 = bodies[1]
eqc1 = eqcs[1]
tra1 = eqc1.constraints[1]
rot1 = eqc1.constraints[2]
eqc2 = eqcs[2]
tra2 = eqc2.constraints[1]
rot2 = eqc2.constraints[2]

rot1.qoffset
rot2.qoffset


################################################################################
# slider
################################################################################
include("conservation_test.jl")
Δt_ = 0.001
mech = getmechanism(:slider, Δt = Δt_, g = 0.0, spring = 0.0, damper = 0.0)
initialize!(mech, :slider)

function controller!(mechanism, k)
    for (i,joint) in enumerate(collect(mechanism.eqconstraints)[1:end])
        if getcontroldim(joint) == 1
            # u = 5e-1 * (rand() - 0.5) * Δt_ * (k < 100)
            if k ∈ (1:500)
                u = 5e-1 * Δt_
            elseif k ∈ (501:1000)
                u = -5e-1 * Δt_
            else
                u = 0.0
            end
            setForce!(mechanism, joint, SA[u])
        end
    end
    return
end

storage = simulate!(mech, 1.0, controller!, record = true, solver = :mehrotra!, verbose = false)
visualize(mech, storage, vis = vis)

momentum(mech)

plot([x[3] for x in storage.x[1]])
plot([ω[1] for ω in storage.ω[1]])

################################################################################
# snake
################################################################################
include("conservation_test.jl")
Δt_ = 0.01
Nlink_ = 5
mech = getmechanism(:snake, Δt = Δt_, g = 0.0, spring = 0.0, damper = 10.0, contact = false, Nlink = Nlink_)
initialize!(mech, :snake, Δv = zeros(3), Δω = zeros(3))

function controller!(mechanism, k)
    for (i,joint) in enumerate(collect(mechanism.eqconstraints)[2:end])
        if k ∈ (1:500)
            u = 5e0 * Δt_
        elseif k ∈ (501:1000)
            u = 5e0 * Δt_
        else
            u = 0.0
        end
        # setForce!(mechanism, joint, SA[0, 0, 0, u, 0, 0])
        setForce!(mechanism, joint, SA[u])
    end
    return
end

storage = simulate!(mech, 1.0, controller!, record = true, solver = :mehrotra!, verbose = false)
visualize(mech, storage, vis = vis)

momentum(mech)[4]

plot([ω[1] for ω in storage.ω[1]])


################################################################################
# npendulum
################################################################################
include("conservation_test.jl")
Δt_ = 0.01
Nlink_ = 2
mech = getmechanism(:npendulum, Δt = Δt_, g = 0.0, spring = 0.0, damper = 0.0, Nlink = Nlink_)
initialize!(mech, :npendulum, ϕ1 = 0.)#, Δv = zeros(3), Δω = zeros(3))

function controller!(mechanism, k)
    n = 50
    for (i,joint) in enumerate(collect(mechanism.eqconstraints)[2:end])
        if  i == Nlink_-1
            if getcontroldim(joint) == 1
                # u = 5e-1 * (rand() - 0.5) * Δt_ * (k < 100)
                if k ∈ (1:n)
                    u = 1e1 * Δt_
                elseif k ∈ (n+1:2n)
                    u = -1e1 * Δt_
                else
                    u = 0.0
                end
                setForce!(mechanism, joint, SA[u])
            end
        end
    end
    return
end

qq = UnitQuaternion(rand(4)...)
τ0 = rand(3)
τ1 = vrotate(vrotate(τ0, qq), inv(qq))
τ0 - τ1


# storage = simulate!(mech, 50*2Δt_, controller!, record = true, solver = :mehrotra!, verbose = true)
storage = simulate!(mech, 2.0, controller!, record = true, solver = :mehrotra!, verbose = true)
visualize(mech, storage, vis = vis)

momentum(mech, collect(mech.bodies)[1])[4]
momentum(mech, collect(mech.bodies)[2])[4]

plot([ω[1] for ω in storage.ω[1]])
plot([ω[1] for ω in storage.ω[2]])

bodies = collect(mech.bodies)
eqcs = collect(mech.eqconstraints)
body1 = bodies[1]
eqc1 = eqcs[1]
tra1 = eqc1.constraints[1]
rot1 = eqc1.constraints[2]
eqc2 = eqcs[2]
tra2 = eqc2.constraints[1]
rot2 = eqc2.constraints[2]

rot1.qoffset
rot2.qoffset




vis = Visualizer()
open(vis)
