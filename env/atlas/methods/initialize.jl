
function getatlas(; Δt::T = 0.01, g::T = -9.81, cf::T = 0.8, spring::T = 0.0, damper::T = 0.0, contact::Bool = true, model_type::Symbol = :simple) where {T}
    path = joinpath(@__DIR__, "../deps/atlas_$(string(model_type)).urdf")
    mech = Mechanism(path, floating=true, g = g, Δt = Δt, spring=spring, damper=damper)

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
        contacts = [
            [-0.1; -0.05; -0.0095],
            [+0.1; -0.05; -0.0095],
            [-0.1; +0.05; -0.0095],
            [+0.1; +0.05; -0.0095],
            ]
        n = length(contacts)
        normal = [[0;0;1.0] for i = 1:n]
        offset = [[0.0; 0.0; 0.01] for i = 1:n]
        cf = cf * ones(T, n)
        names = ["RR", "FR", "RL", "RR"]

        contineqcs1 = contactconstraint(getbody(mech, "l_foot"), normal, cf, p = contacts, offset=offset, names = "l_" .* names)
        contineqcs2 = contactconstraint(getbody(mech, "r_foot"), normal, cf, p = contacts, offset=offset, names = "r_" .* names)

        setPosition!(mech, geteqconstraint(mech, "auto_generated_floating_joint"), [0;0;0.9385;0.;0.;0.])
        mech = Mechanism(origin, bodies, eqs, [contineqcs1; contineqcs2], g = g, Δt = Δt, spring=spring, damper=damper)
    end
    return mech
end

function initializeatlas!(mechanism::Mechanism; 
    tran::AbstractVector{T} = [0,0,0.2],
    rot::AbstractVector{T} = [0,0,0.], 
    v=[zeros(3) for i = 1:length(mechanism.bodies)], 
    ω=[zeros(3) for i = 1:length(mechanism.bodies)], 
    αhip::T = 0.0, αknee::T = 0.0) where {T}
    tran += [0,0,0.9385]

    # positions
    try
        setPosition!(mechanism,
                geteqconstraint(mechanism, "auto_generated_floating_joint"),
                [tran; rot])
        setPosition!(mechanism, geteqconstraint(mechanism, "l_leg_hpxyz"), [0.0, -αhip, 0.0])
        setPosition!(mechanism, geteqconstraint(mechanism, "r_leg_hpxyz"), [0.0, -αhip, 0.0])
        setPosition!(mechanism, geteqconstraint(mechanism, "l_leg_kny"), [αknee])
        setPosition!(mechanism, geteqconstraint(mechanism, "r_leg_kny"), [αknee])
        setPosition!(mechanism, geteqconstraint(mechanism, "l_leg_akxy"), [αhip-αknee, 0.0])
        setPosition!(mechanism, geteqconstraint(mechanism, "r_leg_akxy"), [αhip-αknee, 0.0])
    catch
        nothing
    end

    zeroVelocity!(mechanism)

    return nothing
end
