# Utils
function module_dir()
    return joinpath(@__DIR__, "..", "..")
end

# Activate package
using Pkg
Pkg.activate(module_dir())

using MeshCat

# Open visualizer
vis = Visualizer()
open(vis)

module_dir()
# Include new files
include(joinpath(module_dir(), "examples", "loader.jl"))


function getpendulum(; Δt::T = 0.01, g::T = -9.81, m::T = 1.0, l::T = 1.0,
        spring = 0.0, damper = 0.0, spring_offset = szeros(1), joint_limits = [-sones(1), sones(1)]) where T
    # Parameters
    joint_axis = [1.0; 0; 0]
    width, depth = 0.1, 0.1
    p2 = [0; 0; l/2] # joint connection point

    # Links
    origin = Origin{T}()
    link1 = Box(width, depth, l, m)

    # Constraints
    joint_between_origin_and_link1 = EqualityConstraint(Revolute(origin, link1,
        joint_axis; p2=p2, spring = spring, damper = damper, rot_spring_offset = spring_offset, rot_joint_limits = joint_limits))
    links = [link1]
    eqcs = [joint_between_origin_and_link1]

    mech = Mechanism(origin, links, eqcs, g = g, Δt = Δt, spring=spring, damper=damper)
    return mech
end




mech = getmechanism(:pendulum, Δt = 0.01, g = -9.81)
initialize!(mech, :pendulum, ϕ1 = 0.0)
storage = simulate!(mech, 3.1, record = true, verbose = true)
visualize(mech, storage, vis=vis)



eqcs = collect(mech.eqconstraints)
resetVars!.(eqcs)
eqcs[1].λsol

orig = mech.origin
body1 = collect(mech.bodies)[1]
eqc0 = EqualityConstraint(Revolute(orig, body1, [0,0,1.0], rot_joint_limits = [-sones(1), sones(1)]))
