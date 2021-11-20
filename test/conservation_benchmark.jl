# visualizer
vis = Visualizer()
open(vis)

const Dojo = Main

# Data
ϵ0 = 1e-14
Δt0 = 1e-2
start0 = Int(floor(1/Δt0)) + 1

# Controller
function controller!(mechanism, k; U = 0.10, Δt = Δt0)
    N = Int(floor(1/Δt))
    for (i,joint) in enumerate(mechanism.eqconstraints)
        nu = controldim(joint)
        u = (nu <= 5 && k ∈ (1:N)) * U * Δt * sones(nu)
        setForce!(mechanism, joint, u)
    end
    return
end
nocontrol!(mechanism, k) = controller!(mechanism, k, U = 0.0)

################################################################################
#  ENERGY DRIFT: HUMANOID
################################################################################
# multiple bodies
# initial linear velocities
# no gravity
# no spring and damper
# no control
################################################################################

function energy_drift(vis::Visualizer; Δt = 1e-2, tsim = 1.0, g = 0.0,
        spring = 0.0, damper = 0.0, seed::Int = 100, ϵ = 1e-14, display::Bool = true)
    start = Int(floor(1/Δt)) + 1
    # Test mechanical energy conservation after one second of simulation
    Random.seed!(seed)
    mech = getmechanism(:humanoid, Δt = Δt, g = g, spring = spring, damper = damper, contact = false)
    initialize!(mech, :humanoid)
    # Initialize bodies with random 1m/s velocities
    for body in mech.bodies
        v = rand(3) .- 0.5
        v ./= norm(v)
        setVelocity!(body, v = v)
    end

    tcompute = @elapsed storage = simulate!(mech, tsim + 1.0, nocontrol!, record = true, solver = :mehrotra!, verbose = false, ϵ = ϵ)
    display && visualize(mech, downsample(storage, 1), vis = vis)

    ke = Dojo.kineticEnergy(mech, storage)[start:end]
    pe = Dojo.potentialEnergy(mech, storage)[start:end]
    me = Dojo.mechanicalEnergy(mech, storage)[start:end]

    display && plot([(i-1)*Δt for i in 1:length(ke)], ke .- ke[1])
    display && plot([(i-1)*Δt for i in 1:length(pe)], pe .- pe[1])
    display && plot([(i-1)*Δt for i in 1:length(me)], me .- me[1])

    return ke, pe, me, (tsim + 1.0)/tcompute
end


function energy_drift_experiment(vis::Visualizer; Δts = exp.(Vector(-1.4:-0.2:-4.0) .* log(10)))
    N = length(Δts)
    tratios = zeros(N)
    drifts = zeros(N)
    for (i,Δt) ∈ enumerate(Δts)
        println("$i/$N")
        ke, pe, me, tratio = energy_drift(vis, Δt = Δt, display = false)
        drift = abs(me[end] .- me[1]) / mean(me)
        tratios[i] = tratio
        drifts[i] = drift
    end
    return tratios, drifts
end

ke, pe, me, tratio = energy_drift(vis, Δt = 1e-2)
norm((me .- me[1]) ./ mean(me), Inf)
abs(me[end] .- me[1]) / mean(me)
tratio


dojo_tratios, dojo_drifts = energy_drift_experiment(vis)
mujoco_tratios = [0.5, 10, 100, 500]
mujoco_drifts = [2e-9, 1e-6, 1e-4, 2e-3]

plt = plot(xaxis = :log, yaxis = :log, yflip = true, ylims = (1e-10, 1e-0),
    title = "Energy Drift",
    xlabel = "Faster than real-time multiplier",
    ylabel = "Energy Drift",
    )
scatter!(plt, dojo_tratios, dojo_drifts, label = "Dojo")
scatter!(plt, mujoco_tratios, mujoco_drifts, label = "MuJoCo")



################################################################################
#  MOMENTUM DRIFT: HUMANOID
################################################################################
# multiple bodies
# initial linear velocities
# no gravity
# no spring and damper
# woth control
################################################################################

function momentum_drift(vis::Visualizer; Δt = 1e-2, tsim = 1.0, g = 0.0,
        spring = 0.0, damper = 0.0, seed::Int = 100, ϵ = 1e-14, display::Bool = true)
    start = Int(floor(1/Δt)) + 1
    # Test mechanical energy conservation after one second of simulation
    Random.seed!(seed)
    mech = getmechanism(:humanoid, Δt = Δt, g = g, spring = spring, damper = damper, contact = false)
    initialize!(mech, :humanoid)
    # Initialize bodies with random 1m/s velocities
    for body in mech.bodies
        v = rand(3) .- 0.5
        v ./= norm(v)
        setVelocity!(body, v = v)
    end

    tcompute = @elapsed storage = simulate!(mech, tsim + 1.0, controller!, record = true, solver = :mehrotra!, verbose = false, ϵ = ϵ)
    display && visualize(mech, downsample(storage, 1), vis = vis)

    m = Dojo.momentum(mech, storage)[1:end]
    mlin = [Vector(mi)[1:3] for mi in m]
    mang = [Vector(mi)[4:6] for mi in m]

    display && plot([(i-1)*Δt for i in 1:length(m)], hcat(mlin...)')
    display && plot([(i-1)*Δt for i in 1:length(m)], hcat(mang...)')

    return mlin, mang, (tsim + 1.0)/tcompute
end

function momentum_drift_experiment(vis::Visualizer; Δts = exp.(Vector(-1.4:-0.2:-4.0) .* log(10)))
    N = length(Δts)
    tratios = zeros(N)
    linear_drifts = zeros(N)
    angular_drifts = zeros(N)
    for (i,Δt) ∈ enumerate(Δts)
        println("$i/$N")
        mlin, mang, tratio = momentum_drift(vis, Δt = Δt, display = false)
        linear_drift = norm(mlin[end] .- mlin[1]) / mean(norm.(mlin))
        angular_drift = norm(mang[end] .- mang[1]) / mean(norm.(mang))
        tratios[i] = tratio
        linear_drifts[i] = linear_drift
        angular_drifts[i] = angular_drift
    end
    return tratios, linear_drifts, angular_drifts
end

# vis = Visualizer()
# open(vis)

mlin, mang, tratio = momentum_drift(vis, Δt = 1e-2)
norm(mlin[end] .- mlin[1]) / mean(norm.(mlin))
norm(mang[end] .- mang[1]) / mean(norm.(mang))
tratio

# dojo_tratios, dojo_linear_drifts, dojo_angular_drifts = momentum_drift_experiment(vis)
physx_tratios = [0.1, 1, 10, 60]
physx_linear_drifts = [4e-5, 2e-5, 8e-6, 2e-6]
brax_tratios = [1, 10, 100, 400]
brax_angular_drifts = [5e-4, 1e-3, 2e-3, 4e-3]

# Linear Momentum Drift
plt = plot(xaxis = :log, yaxis = :log, yflip = true, ylims = (1e-17, 1e-0),
    title = "Linear Momentum Drift",
    xlabel = "Faster than real-time multiplier",
    ylabel = "Linear Momentum Drift",
    )
scatter!(plt, dojo_tratios, dojo_linear_drifts, label = "Dojo")
scatter!(plt, physx_tratios, physx_linear_drifts, label = "PhysX")


# Angular Momentum Drift
plt = plot(xaxis = :log, yaxis = :log, yflip = true, ylims = (1e-17, 1e-0),
    title = "Angular Momentum Drift",
    xlabel = "Faster than real-time multiplier",
    ylabel = "Angular Momentum Drift",
    )
scatter!(plt, dojo_tratios, dojo_angular_drifts, label = "Dojo")
scatter!(plt, brax_tratios, brax_angular_drifts, label = "Brax")




################################################################################
#  ENERGY DRIFT: HUMANOID
################################################################################
# multiple bodies
# initial linear velocities
# no gravity
# no spring and damper
# no control
################################################################################

function energy_drift(vis::Visualizer; Δt = 1e-2, tsim = 1.0, g = 0.0,
        spring = 0.0, damper = 0.0, seed::Int = 100, ϵ = 1e-14, display::Bool = true)
    start = Int(floor(1/Δt)) + 1
    # Test mechanical energy conservation after one second of simulation
    Random.seed!(seed)
    mech = getmechanism(:humanoid, Δt = Δt, g = g, spring = spring, damper = damper, contact = false)
    initialize!(mech, :humanoid)
    # Initialize bodies with random 1m/s velocities
    for body in mech.bodies
        v = rand(3) .- 0.5
        v ./= norm(v)
        setVelocity!(body, v = v)
    end

    tcompute = @elapsed storage = simulate!(mech, tsim + 1.0, nocontrol!, record = true, solver = :mehrotra!, verbose = false, ϵ = ϵ)
    display && visualize(mech, downsample(storage, 1), vis = vis)

    ke = Dojo.kineticEnergy(mech, storage)[start:end]
    pe = Dojo.potentialEnergy(mech, storage)[start:end]
    me = Dojo.mechanicalEnergy(mech, storage)[start:end]

    display && plot([(i-1)*Δt for i in 1:length(ke)], ke .- ke[1])
    display && plot([(i-1)*Δt for i in 1:length(pe)], pe .- pe[1])
    display && plot([(i-1)*Δt for i in 1:length(me)], me .- me[1])

    return ke, pe, me, (tsim + 1.0)/tcompute
end


function energy_drift_experiment(vis::Visualizer; Δts = exp.(Vector(-1.4:-0.2:-4.0) .* log(10)))
    N = length(Δts)
    tratios = zeros(N)
    drifts = zeros(N)
    for (i,Δt) ∈ enumerate(Δts)
        println("$i/$N")
        ke, pe, me, tratio = energy_drift(vis, Δt = Δt, display = false)
        drift = abs(me[end] .- me[1]) / mean(me)
        tratios[i] = tratio
        drifts[i] = drift
    end
    return tratios, drifts
end

ke, pe, me, tratio = energy_drift(vis, Δt = 1e-2)
norm((me .- me[1]) ./ mean(me), Inf)
abs(me[end] .- me[1]) / mean(me)
tratio


dojo_tratios, dojo_drifts = energy_drift_experiment(vis)
mujoco_tratios = [0.5, 10, 100, 500]
mujoco_drifts = [2e-9, 1e-6, 1e-4, 2e-3]

plt = plot(xaxis = :log, yaxis = :log, yflip = true, ylims = (1e-10, 1e-0),
    title = "Energy Drift",
    xlabel = "Faster than real-time multiplier",
    ylabel = "Energy Drift",
    )
scatter!(plt, dojo_tratios, dojo_drifts, label = "Dojo")
scatter!(plt, mujoco_tratios, mujoco_drifts, label = "MuJoCo")

################################################################################
# LONG TERM ENERGY BEHAVIOR: DRIFT VS OSCILLATIONS
################################################################################

function long_term_experiment(vis::Visualizer; tsim = 300.0, Δts = exp.(Vector(-1.2:-0.4:-2.0) .* log(10)))
    N = length(Δts)
    mes = []
    for (i,Δt) ∈ enumerate(Δts)
        println("$i/$N")
        ke, pe, me, tratio = energy_drift(vis, tsim = tsim, Δt = Δt, display = false)
        push!(mes, me)
    end
    return Δts, mes
end

# Δts, mes = long_term_experiment(vis)
N = length(Δts)

plt = plot(ylims = (0.0, Inf),
    title = "Energy",
    layout = (1,N),
    )

for j = 1:N
    plot!(plt[1,j], [(i-1)*Δts[j] for i = 1:length(mes[j])], mes[j] ./ mean(mes[j]),
        label = "Dojo", xlabel = "Freq. " * string(Int(floor(1/Δts[j]))) * "Hz")
end
display(plt)
