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



Δt_ = 0.01
mech = getmechanism(:humanoid, contact = true, Δt = Δt_, g = -9.81, spring = 100.0, damper = 10.)
initialize!(mech, :humanoid, rot = [0,0,0.1], tran = [0,0,1.3])
eqcs = collect(mech.eqconstraints)

function controller!(mechanism, k)
    for (i,eqc) in enumerate(collect(mechanism.eqconstraints)[2:end])
        pbody = getbody(mech, eqc.parentid)
        minJ = minimum(diag(pbody.J))
        for (i,joint) in enumerate(eqc.constraints)
            cbody = getbody(mech, eqc.childids[i])
            minJ = min(minJ, minimum(diag(cbody.J)))
        end
        nu = getcontroldim(eqc)
        u = 10 * minJ * (ones(nu) .- 0.2) * Δt_
        setForce!(mechanism, eqc, SVector{nu}(u))
    end
    return
end

@elapsed storage = simulate!(mech, 1.5, controller!, record = true, solver = :mehrotra!, verbose = false)
visualize(mech, storage, vis = vis)

collect(mech.eqconstraints)[1]
collect(mech.eqconstraints)[1].isdamper
collect(mech.eqconstraints)[2]
collect(mech.eqconstraints)[2].isdamper
collect(mech.eqconstraints)[2].constraints[1].damper

# filename = "humanoid_mujoco"
# MeshCat.convert_frames_to_video(
#     "/home/simon/Downloads/$filename.tar",
#     "/home/simon/Documents/video/$filename.mp4", overwrite=true)
#
# using FFMPEG
# convert_video_to_gif(
#     "/home/simon/Documents/video/$filename.mp4",
#     "/home/simon/Documents/video/$filename.gif", overwrite=true)
