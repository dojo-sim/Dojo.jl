function getquadruped(; Δt::T=0.01, g::T=-9.81, cf::T=0.8, spring=0.0,
    damper=0.0, contact::Bool=true, path=joinpath(@__DIR__, "../deps/quadruped.urdf")) where {T}
    mech=Mechanism(path, true, T, g=g, Δt=Δt, spring=spring, damper=damper)

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

        ineqcs1 = contactconstraint(getbody(mech,:FR_calf), normal; cf=cf, p=contact, name=:FR_contact)
        ineqcs2 = contactconstraint(getbody(mech,:FL_calf), normal; cf=cf, p=contact, name=:FL_contact)
        ineqcs3 = contactconstraint(getbody(mech,:RR_calf), normal; cf=cf, p=contact, name=:RR_contact)
        ineqcs4 = contactconstraint(getbody(mech,:RL_calf), normal; cf=cf, p=contact, name=:RL_contact)
        setPosition!(mech, geteqconstraint(mech, :auto_generated_floating_joint), [0;0;0.23;0.;0.;0.])
        mech = Mechanism(origin, bodies, eqs, [ineqcs1; ineqcs2; ineqcs3; ineqcs4], g=g, Δt=Δt, spring=spring, damper=damper)
    end
    return mech
end

function initializequadruped!(mechanism::Mechanism; tran::AbstractVector{T}=[0,0,0.],
    rot::AbstractVector{T}=[0,0,0.0], v::AbstractVector{T}=[0,0,0.0], θ::T=0.95) where {T}
    tran += [0,0,0.31]
    setPosition!(mechanism, geteqconstraint(mechanism, :auto_generated_floating_joint), [tran; rot])

    setPosition!(mechanism, geteqconstraint(mechanism, :FR_thigh_joint), [θ])
    setPosition!(mechanism, geteqconstraint(mechanism, :FR_calf_joint), [-1.5*θ])

    setPosition!(mechanism, geteqconstraint(mechanism, :FL_thigh_joint), [θ*0.9])
    setPosition!(mechanism, geteqconstraint(mechanism, :FL_calf_joint), [-1.5*θ])

    setPosition!(mechanism, geteqconstraint(mechanism, :RR_thigh_joint), [θ*0.9])
    setPosition!(mechanism, geteqconstraint(mechanism, :RR_calf_joint), [-1.5*θ])

    setPosition!(mechanism, geteqconstraint(mechanism, :RL_thigh_joint), [θ])
    setPosition!(mechanism, geteqconstraint(mechanism, :RL_calf_joint), [-1.5*θ])

    setVelocity!(mechanism, geteqconstraint(mechanism, :auto_generated_floating_joint), [v; zeros(3)])

    setVelocity!(mechanism, geteqconstraint(mechanism, :FR_thigh_joint), [0.])
    setVelocity!(mechanism, geteqconstraint(mechanism, :FR_calf_joint), [0.])

    setVelocity!(mechanism, geteqconstraint(mechanism, :FL_thigh_joint), [0.])
    setVelocity!(mechanism, geteqconstraint(mechanism, :FL_calf_joint), [0.])

    setVelocity!(mechanism, geteqconstraint(mechanism, :RR_thigh_joint), [0.])
    setVelocity!(mechanism, geteqconstraint(mechanism, :RR_calf_joint), [0.])

    setVelocity!(mechanism, geteqconstraint(mechanism, :RL_thigh_joint), [0.])
    setVelocity!(mechanism, geteqconstraint(mechanism, :RL_calf_joint), [0.])
end
