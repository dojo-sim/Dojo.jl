
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




include("conservation_test.jl")
Δt_ = 0.01
mech = getmechanism(:atlas, Δt = Δt_, g = 0.0, spring = 000.0, damper = 400.0, contact = false)
initialize!(mech, :atlas, tran = [0, 0, 2.])

function controller!(mechanism, k)
    for (i,joint) in enumerate(collect(mechanism.eqconstraints)[2:end])
        if getcontroldim(joint) == 1
            u = 5e-1 * (rand() - 0.5) * Δt_ * (k < 100)
            @show u
            setForce!(mechanism, joint, SA[u])
        end
    end
    return
end

storage = simulate!(mech, 10.0, controller!, record = true, solver = :mehrotra!, verbose = false)
visualize(mech, storage, vis = vis)


momentum(mech)
