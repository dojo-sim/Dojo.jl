function gethalfcheetah(; Δt::T = 0.01, g::T = -9.81, cf::T = 0.8, spring::T = 0.0, damper::T = 0.0, contact::Bool = true) where {T}
    # TODO new feature: visualize capsule instead of cylinders
    # TODO new feature: visualize multiple shapes for a single body
    path = joinpath(@__DIR__, "../deps/halfcheetah.urdf")
    mech = Mechanism(path, floating=false, g = g, Δt = Δt)

    # Adding springs and dampers
    for (i,eqc) in enumerate(collect(mech.eqconstraints[2:end]))
    # @warn "damping on origin joint"
    # for (i,eqc) in enumerate(collect(mech.eqconstraints[1:end]))
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

        setPosition!(mech, geteqconstraint(mech, "floating_joint"), [0.530509,0,0.02792])
        mech = Mechanism(origin, bodies, eqs, [contineqcs1; contineqcs2], g = g, Δt = Δt)
    end
    return mech
end

function initializehalfcheetah!(mechanism::Mechanism; x::T = 0.0, z::T = 0.0, θ::T = 0.0) where {T}
    setPosition!(mechanism,
                 geteqconstraint(mechanism, "floating_joint"),
                 [z + 0.530509, -x, -θ + 0.02792])

end