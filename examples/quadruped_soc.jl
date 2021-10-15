using ConstrainedDynamics
using ConstrainedDynamicsVis
using StaticArrays
using LinearAlgebra

# Utils
function module_dir()
    return joinpath(@__DIR__, "..")
end

# Activate package
using Pkg
Pkg.activate(module_dir())

# Load packages
# using ConstrainedDynamics
# using ConstrainedDynamicsVis
using Plots
using Random


# Include new files
include(joinpath(module_dir(), "examples", "loader.jl"))


begin
    path = "examples/examples_files/quadruped_simple.urdf"
    mech = Mechanism(joinpath(module_dir(), path), floating=false, g = 0)
    origin = mech.origin
    bodies = mech.bodies.values
    eqcs = mech.eqconstraints.values

    fricsandineqsFR = Friction(getbody(mech,"FR_calf"), [0;0;1.0], 0.2; p = [0.0;0;-0.1])
    fricsandineqsFL = Friction(getbody(mech,"FL_calf"), [0;0;1.0], 0.2; p = [0.0;0;-0.1])
    fricsandineqsRR = Friction(getbody(mech,"RR_calf"), [0;0;1.0], 0.2; p = [0.0;0;-0.1])
    fricsandineqsRL = Friction(getbody(mech,"RL_calf"), [0;0;1.0], 0.2; p = [0.0;0;-0.1])
    fricsFR = fricsandineqsFR[1]
    fricsFL = fricsandineqsFL[1]
    fricsRR = fricsandineqsRR[1]
    fricsRL = fricsandineqsRL[1]
    ineqcFR = fricsandineqsFR[2]
    ineqcFL = fricsandineqsFL[2]
    ineqcRR = fricsandineqsRR[2]
    ineqcRL = fricsandineqsRL[2]

    frics = [fricsFR;fricsFL;fricsRR;fricsRL]
    ineqcs = [ineqcFR;ineqcFL;ineqcRR;ineqcRL]
    imps = ineqcs[findall(x -> typeof(x.constraints[1]) <: Impact, ineqcs)]

    calves = ["FR_calf", "FL_calf", "RR_calf", "RL_calf"]
    coneineqcs = [ConeBound(getbody(mech, calf), [0;0;1.0], 0.2; p = [0.0;0;-0.1]) for calf in calves]
    conebounds = getindex.(coneineqcs,1)
    impactbounds = getindex.(coneineqcs,2)
    impactids = getfield.(impactbounds, :id)
    conebounds = [InequalityConstraint((cb, getbody(mech, calves[i]).id, impactids[i])) for (i,cb) in enumerate(conebounds)]
end

Δt_ = 0.001
# mech2 = Mechanism(origin, bodies, eqcs, ineqcs, frics, g = -9.81, Δt = Δt_)
# mech2 = Mechanism(origin, bodies, eqcs, imps, g = -9.81, Δt = Δt_)
mech2 = Mechanism(origin, bodies, eqcs, [impactbounds; conebounds], g = -9.81, Δt = Δt_)

begin
    setPosition!(mech2, geteqconstraint(mech2,"floating_base"),[0;0;0.23;0.;0.;0.])

    initangle = 0.95
    setPosition!(mech2, geteqconstraint(mech2,"FR_thigh_joint"),[initangle])
    setPosition!(mech2, geteqconstraint(mech2,"FR_calf_joint"),[-2*initangle])

    setPosition!(mech2, geteqconstraint(mech2,"FL_thigh_joint"),[initangle*0.9])
    setPosition!(mech2, geteqconstraint(mech2,"FL_calf_joint"),[-2*initangle])

    setPosition!(mech2, geteqconstraint(mech2,"RR_thigh_joint"),[initangle*0.9])
    setPosition!(mech2, geteqconstraint(mech2,"RR_calf_joint"),[-2*initangle])

    setPosition!(mech2, geteqconstraint(mech2,"RL_thigh_joint"),[initangle])
    setPosition!(mech2, geteqconstraint(mech2,"RL_calf_joint"),[-2*initangle])

    N = 100
    steps = Base.OneTo(N)
    storage = Storage(steps,length(bodies))

    traj21 = [initangle*(cos(i*0.001*2*pi)*0.1+0.9) for i=1:N]
    traj21 = [[initangle for i=1:500];traj21]
    traj21 = SA[traj21...]
    traj31 = [-2*initangle-sin(i*0.001*2*pi)*0.1 for i=1:N]
    traj31 = [[-2*initangle for i=1:500];traj31]
    traj31 = SA[traj31...]

    traj22 = [initangle*(cos(i*0.001*2*pi+pi)*0.1+0.9) for i=1:N]
    traj22 = [[initangle*0.9 for i=1:500];traj22]
    traj22 = SA[traj22...]
    traj32 = [-2*initangle-sin(i*0.001*2*pi+pi)*0.1 for i=1:N]
    traj32 = [[-2*initangle for i=1:500];traj32]
    traj32 = SA[traj32...]



    function singleleg(mechanism, leg, angles)
        j1 = geteqconstraint(mechanism, leg*"_hip_joint")
        j2 = geteqconstraint(mechanism, leg*"_thigh_joint")
        j3 = geteqconstraint(mechanism, leg*"_calf_joint")

        θ1 = minimalCoordinates(mechanism, j1)[1]
        θ2 = minimalCoordinates(mechanism, j2)[1]
        θ3 = minimalCoordinates(mechanism, j3)[1]
        dθ1 = minimalVelocities(mechanism, j1)[1]
        dθ2 = minimalVelocities(mechanism, j2)[1]
        dθ3 = minimalVelocities(mechanism, j3)[1]

        u1 = (100.0*(angles[1]-θ1) + 5.0*(0-dθ1)) * Δt_ #* 0.17
        u2 = (80.0*(angles[2]-θ2) + 4.0*(0-dθ2)) * Δt_ #* 0.17
        u3 = (60.0*(angles[3]-θ3) + 3.0*(0-dθ3)) * Δt_ #* 0.17

        setForce!(mechanism, j1, SA[u1])
        setForce!(mechanism, j2, SA[u2])
        setForce!(mechanism, j3, SA[u3])
    end

    function controller!(mechanism, k)
        singleleg(mechanism, "FR", SA[0.0;traj21[k];traj31[k]])
        singleleg(mechanism, "FL", SA[0.0;traj22[k];traj32[k]])
        singleleg(mechanism, "RR", SA[0.0;traj22[k];traj32[k]])
        singleleg(mechanism, "RL", SA[0.0;traj21[k];traj31[k]])
    end
end


# include(joinpath(module_dir(), "src", "solver", "mehrotra.jl"))
storage = simulate!(mech2, storage, controller!, record = true,
    debug=false, ε = 1e-3,
    solver = :mehrotra!)
visualize(mech2, storage, showframes = false)

plot(hcat(Vector.(storage.x[1])...)')
plot(hcat([[q.w, q.x, q.y, q.z] for q in storage.q[1]]...)')
plot(hcat(Vector.(storage.v[1])...)')
plot(hcat(Vector.(storage.ω[1])...)')


using FFMPEG
using MeshCat
filename = "atlas_rigid"
MeshCat.convert_frames_to_video(
    "/home/simon/Downloads/$filename.tar",
    "/home/simon/Documents/video/$filename.mp4", overwrite=true)

convert_video_to_gif(
    "/home/simon/Documents/video/$filename.mp4",
    "/home/simon/Documents/video/$filename.gif", overwrite=true)
