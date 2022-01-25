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



mech = getmechanism(:orbital, spring = 0.0, damper = 1.0, Nb = 2)
initialize!(mech, :orbital, ϕx = 3π/4, ϕy = 3π/4)

function controller!(mechanism, k)
    for (i,eqc) in enumerate(collect(mechanism.eqconstraints)[2:end])
        # pbody = get_body(mech, eqc.parentid)
        # minJ = minimum(diag(pbody.J))
        # for (i,joint) in enumerate(eqc.constraints)
        #     cbody = get_body(mech, eqc.childids[i])
        #     minJ = min(minJ, minimum(diag(cbody.J)))
        # end
        minJ = 0.01
        nu = control_dimension(eqc)
        u = 10 * minJ * (rand(nu) .- 0.2) * Δt_
        set_input!(eqc, SVector{nu}(u))
    end
    return
end

@elapsed storage = simulate!(mech, 5, controller!, record = true, solver = :mehrotra!, verbose = false)
visualize(mech, storage, vis = vis)

collect(mech.eqconstraints)[1]
collect(mech.eqconstraints)[1].isdamper
collect(mech.eqconstraints)[2]
collect(mech.eqconstraints)[2].isdamper
collect(mech.eqconstraints)[2].constraints[1].damper
