function getquadruped(; Δt::T = 0.01, g::T = -9.81, cf::T = 0.8, spring = 0.0,
    damper = 0.0, contact::Bool = true) where {T}
    path = joinpath(@__DIR__, "../deps/quadruped_simon.urdf")
    mech = Mechanism(path, floating = true, g = g, Δt = Δt, spring=spring, damper=damper)

    # Adding springs and dampers
    for (i,eqc) in enumerate(collect(mech.eqconstraints)[2:end])
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
        contact = [0.0;0;-0.1]
        normal = [0;0;1.0]

        contineqcs1 = contactconstraint(getbody(mech,"FR_calf"), normal, cf; p = contact, name = "FR_contact")
        contineqcs2 = contactconstraint(getbody(mech,"FL_calf"), normal, cf; p = contact, name = "FL_contact")
        contineqcs3 = contactconstraint(getbody(mech,"RR_calf"), normal, cf; p = contact, name = "RR_contact")
        contineqcs4 = contactconstraint(getbody(mech,"RL_calf"), normal, cf; p = contact, name = "RL_contact")
        # contineqcs1 = impactconstraint(getbody(mech,"FR_calf"), normal; p = contact,)
        # contineqcs2 = impactconstraint(getbody(mech,"FL_calf"), normal; p = contact,)
        # contineqcs3 = impactconstraint(getbody(mech,"RR_calf"), normal; p = contact,)
        # contineqcs4 = impactconstraint(getbody(mech,"RL_calf"), normal; p = contact,)
        setPosition!(mech, geteqconstraint(mech, "auto_generated_floating_joint"), [0;0;0.23;0.;0.;0.])
        mech = Mechanism(origin, bodies, eqs, [contineqcs1; contineqcs2; contineqcs3; contineqcs4], g = g, Δt = Δt, spring=spring, damper=damper)
    end
    return mech
end

function initializequadruped!(mechanism::Mechanism; tran::AbstractVector{T} = [0,0,0.],
    rot::AbstractVector{T} = [0,0,0.0], v::AbstractVector{T} = [0,0,0.0], θ::T = 0.95) where {T}
    tran += [0,0,0.31]
    setPosition!(mechanism, geteqconstraint(mechanism, "auto_generated_floating_joint"), [tran; rot])

    setPosition!(mechanism, geteqconstraint(mechanism, "FR_thigh_joint"), [θ])
    setPosition!(mechanism, geteqconstraint(mechanism, "FR_calf_joint"), [-1.5*θ])

    setPosition!(mechanism, geteqconstraint(mechanism, "FL_thigh_joint"), [θ*0.9])
    setPosition!(mechanism, geteqconstraint(mechanism, "FL_calf_joint"), [-1.5*θ])

    setPosition!(mechanism, geteqconstraint(mechanism, "RR_thigh_joint"), [θ*0.9])
    setPosition!(mechanism, geteqconstraint(mechanism, "RR_calf_joint"), [-1.5*θ])

    setPosition!(mechanism, geteqconstraint(mechanism, "RL_thigh_joint"), [θ])
    setPosition!(mechanism, geteqconstraint(mechanism, "RL_calf_joint"), [-1.5*θ])

    setVelocity!(mechanism, geteqconstraint(mechanism, "auto_generated_floating_joint"), [v; zeros(3)])

    setVelocity!(mechanism, geteqconstraint(mechanism, "FR_thigh_joint"), [0.])
    setVelocity!(mechanism, geteqconstraint(mechanism, "FR_calf_joint"), [0.])

    setVelocity!(mechanism, geteqconstraint(mechanism, "FL_thigh_joint"), [0.])
    setVelocity!(mechanism, geteqconstraint(mechanism, "FL_calf_joint"), [0.])

    setVelocity!(mechanism, geteqconstraint(mechanism, "RR_thigh_joint"), [0.])
    setVelocity!(mechanism, geteqconstraint(mechanism, "RR_calf_joint"), [0.])

    setVelocity!(mechanism, geteqconstraint(mechanism, "RL_thigh_joint"), [0.])
    setVelocity!(mechanism, geteqconstraint(mechanism, "RL_calf_joint"), [0.])
end