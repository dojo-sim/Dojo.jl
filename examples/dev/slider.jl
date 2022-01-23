# Utils
function module_dir()
    return joinpath(@__DIR__, "..", "..")
end

# Activate package
using Pkg
Pkg.activate(module_dir())

# Load packages
using MeshCat

# Open visualizer
vis = Visualizer()
open(vis)

# Include new files
include(joinpath(module_dir(), "examples", "loader.jl"))


mech = getmechanism(:slider, Δt = 0.0005, g = -0.00, spring = 10.0, damper = 0.3);
initialize!(mech, :slider, z1 = 0.0)
function ctrl!(mechanism, k)
    u = [10]*mechanism.Δt
    set_control!(mechanism, u)
    return nothing
end

@elapsed storage = simulate!(mech, 4.00, ctrl!, record = true, verbose = false,
    opts=InteriorPointOptions(verbose=false, btol = 1e-6))
visualize(mech, storage, vis = vis)

# plt = plot()
alt = [s[3] for s in storage.x[1]]
tim = Vector(0:mech.Δt:10.00)[1:length(alt)]
plot!(plt, tim, alt, label = mech.Δt)



mech = getmechanism(:pendulum, Δt = 0.1, g = -0.00, spring = 10.0, damper = 0.3);
initialize!(mech, :pendulum, ϕ1 = 0.0)
function ctrl!(mechanism, k)
    u = [5]*mechanism.Δt
    set_control!(mechanism, u)
    return nothing
end

@elapsed storage = simulate!(mech, 4.00, ctrl!, record = true, verbose = false,
    opts=InteriorPointOptions(verbose=false, btol = 1e-6))
visualize(mech, storage, vis = vis)

# plt = plot()
alt = [s[3] for s in storage.x[1]]
tim = Vector(0:mech.Δt:10.00)[1:length(alt)]
plot!(plt, tim, alt, label = mech.Δt)
