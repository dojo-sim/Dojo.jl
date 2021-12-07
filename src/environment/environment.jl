import Base.reset

mutable struct Environement12{M}
    mechanism::M
end

function make(model::String; kwargs...)
    mechanism = getmechanism(Symbol(model); kwargs...)
    env = Environement12(mechanism)
    reset(env)
    return env
end

function reset(env::Environement12{M}; kwargs...) where {M}
    initialize!(mechanism; kwargs...)
    return nothing
end



# script
env = make("dice", Δt = 0.05, g = -9.81);
reset(env)





################################################################################
# Development
################################################################################
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

mech = getmechanism(:dice, Δt = 0.01, g = -9.81, cf = 0.2, contact = true, mode=:box, conetype = :soc)
# mech = getmechanism(:dice, Δt = 0.01, g = -9.81, cf = 0.2, contact = true, mode=:box, conetype = :linear)
# mech = getmechanism(:dice, Δt = 0.01, g = -9.81, contact = true, mode=:box, conetype = :impact)
Random.seed!(100)
ω = 0.0 * (rand(3) .- 0.5) * 1
x = [0, 0, 1.0]
v = 1.0 * [4, 1, 0.0]
initialize!(mech, :dice, x = x, v = v, ω = ω)
storage = simulate!(mech, 1.3, record = true, solver = :mehrotra!, verbose = false)
visualize(mech, storage, vis = vis)
