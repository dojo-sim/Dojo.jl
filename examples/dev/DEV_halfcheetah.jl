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


mech = getmechanism(:halfcheetah, Nlink = 5, Δt = 0.01, g = -9.81, cf = 0.0, contact = true, conetype = :soc)

x = [0,-0.5,.1]
v = 0.1*[1,.3,4]
ω = 1.5*[.1,.8,0]
q1 = UnitQuaternion(RotX(π/1.5))
initialize!(mech, :snake, x = x, v = v, ω = ω, q1 = q1)

@elapsed storage = simulate!(mech, .28, record = true, solver = :mehrotra!)

visualize(mech, storage, vis = vis)



function gethalfcheetah(; Δt::T = 0.01, g::T = -9.81, cf::T = 0.8, spring::T = 0.0, damper::T = 0.0, contact::Bool = true) where {T}
    # TODO new feature: visualize capsule instead of cylinders
    # TODO new feature: visualize multiple shapes for a single body
    path = "examples/examples_files/halfcheetah.urdf"
    mech = Mechanism(joinpath(module_dir(), path), floating=false, g = g, Δt = Δt)

    # Adding springs and dampers
    for (i,eqc) in enumerate(collect(mech.eqconstraints[2:end]))
        eqc.isdamper = true
        eqc.isspring = true
        for joint in eqc.constraints
            joint.spring = spring
            joint.damper = damper
        end
    end

    if contact
        origin = Origin{T}()
        bodies = Vector{Body{T}}(collect(mech.bodies))
        eqs = Vector{EqualityConstraint{T}}(collect(mech.eqconstraints))

        # Foot contact
        contact1 = [0.0; 0.0; -0.140]
        contact2 = [0.0; 0.0; -0.188]
        normal = [0;0;1.0]

        contineqcs1 = contactconstraint(getbody(mech, "ffoot"), normal, cf, p = contact1)
        contineqcs2 = contactconstraint(getbody(mech, "bfoot"), normal, cf, p = contact2)

        setPosition!(mech, geteqconstraint(mech, "floating_joint"), [1.0,0,0])
        mech = Mechanism(origin, bodies, eqs, [contineqcs1; contineqcs2], g = g, Δt = Δt)
    end
    return mech
end

function initializehalfcheetah!(mechanism::Mechanism; x::T = 0.0, z::T = 1.2, θ::T = 0.1) where {T}
    setPosition!(mechanism,
                 geteqconstraint(mechanism, "floating_joint"),
                 [z, -x, -θ])
end

mech = getmechanism(:halfcheetah, g = -9.81, damper = 5.0, spring = 300.0)
initialize!(mech, :halfcheetah, x = 1.0, z = 1.0, θ = 0.2)

@elapsed storage = simulate!(mech, 1.5, record = true, solver = :mehrotra!, verbose = false)
visualize(mech, storage, vis = vis)


bodies = collect(mech.bodies)
eqcs = collect(mech.eqconstraints)

bodies
eqcs
