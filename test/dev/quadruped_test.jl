# Utils
function module_dir()
    return joinpath(@__DIR__, "..")
end

# Activate package
using Pkg
Pkg.activate(module_dir())

using MeshCat
# using IterativeLQR

# Open visualizer
vis = Visualizer()
open(vis)

# Include new files
include(joinpath(module_dir(), "examples", "loader.jl"))



mech = getmechanism(:quadruped, Δt = 0.05, g = -9.81, cf = 0.8, damper = 0.1, spring = 0*300.0)
initialize!(mech, :quadruped, tran = [0,0,0.04], rot = 2*[0.1,0.2,0.3])
@elapsed storage = simulate!(mech, 0.50, record = true, solver = :mehrotra!, verbose = true)
visualize(mech, storage, vis = vis)

ss = get_sdf(mech, storage)
plot(hcat(ss...))

# force = 1.0 / 0.05
# m = 0.20 * force
# # I * a = m
# a = m/ 0.003

# γ1 = [1,0,0,0.0]
# s1 = [1,0,0,0.0]
# γ0, s0 = initial_state_ort(γ, s; ϵ = 1e-20)
# γ0, s0 = initial_state_soc(γ, s; ϵ = 1e-20)
#
# ineqcs = mech.ineqconstraints
# ineqc1 = mech.ineqconstraints.values[1]
# ineqc2 = mech.ineqconstraints.values[2]
# ineqc3 = mech.ineqconstraints.values[3]
# ineqc4 = mech.ineqconstraints.values[4]
# resetVars!.(ineqcs, scaling = 0.05)
# ineqc1.ssol[2]
# ineqc1.γsol[2]
# initial_state!.(ineqcs.values)
# ineqc1.ssol[2]
# ineqc1.γsol[2]
