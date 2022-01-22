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



Δt_ = 0.05
mech = getmechanism(:humanoid, contact = true, Δt = Δt_, g = -9.81, spring = 500.0, damper = 50.)
initialize!(mech, :humanoid, rot = [0.1,0,0], tran = [0,0,1.5])
eqcs = collect(mech.eqconstraints)

function controller!(mechanism, k)
    for (i,eqc) in enumerate(collect(mechanism.eqconstraints)[2:end])
        pbody = getbody(mech, eqc.parentid)
        minJ = minimum(diag(pbody.J))
        for (i,joint) in enumerate(eqc.constraints)
            cbody = getbody(mech, eqc.childids[i])
            minJ = min(minJ, minimum(diag(cbody.J)))
        end
        nu = controldim(eqc)
        u = 10 * minJ * (ones(nu) .- 0.2) * Δt_
        setForce!(mechanism, eqc, SVector{nu}(u))
    end
    return
end

@elapsed storage = simulate!(mech, 2.3, controller!, record = true, solver = :mehrotra!, verbose = false)
visualize(mech, storage, vis = vis)
