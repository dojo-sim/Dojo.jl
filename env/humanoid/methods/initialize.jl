function gethumanoid(; Δt::T=0.01, g::T=-9.81, cf=0.8, spring=0.0, damper=0.0, contact::Bool=true, contact_body::Bool=false) where {T}
    path = joinpath(@__DIR__, "../deps/humanoid.urdf")
    mech = Mechanism(path, floating=true, g=g, Δt=Δt, spring=spring, damper=damper)

    if contact
        origin = Origin{T}()
        bodies = Vector{Body{T}}(collect(mech.bodies))
        eqs = Vector{EqualityConstraint{T}}(collect(mech.eqconstraints))

        # Foot contact
        left_foot = getbody(mech, "left_foot") 

        pfl = [0.5 * left_foot.shape.shape[1].rh[2]; 0.0; 0.0]
        ofl = [0.0; 0.0; left_foot.shape.shape[1].rh[1]]
        pbl = [-0.5 * left_foot.shape.shape[1].rh[2]; 0.0; 0.0]
        obl = [0.0; 0.0; left_foot.shape.shape[1].rh[1]]
        contacts = [
                    pfl,
                    pbl,
                   ]
        offsets = [
                    ofl,
                    obl,
                  ]
        n = length(contacts)
        normal = [[0;0;1.0] for i = 1:n]
        cfs = cf * ones(T, n)

        contineqcs_left = contactconstraint(left_foot, normal, cfs, p = contacts, offset=offsets)

        right_foot = getbody(mech, "right_foot") 

        pfr = [0.5 * right_foot.shape.shape[1].rh[2]; 0.0; 0.0]
        ofr = [0.0; 0.0; right_foot.shape.shape[1].rh[1]]
        pbr = [-0.5 * right_foot.shape.shape[1].rh[2]; 0.0; 0.0]
        obr = [0.0; 0.0; right_foot.shape.shape[1].rh[1]]
        contacts = [
                    pfr,
                    pbr,
                   ]
        offsets = [
                    ofr,
                    obr,
                  ]
        n = length(contacts)
        normal = [[0;0;1.0] for i = 1:n]
        cfs = cf * ones(T, n)

        contineqcs_right = contactconstraint(right_foot, normal, cfs, p = contacts, offset=offsets)

        setPosition!(mech, geteqconstraint(mech, "auto_generated_floating_joint"), [0;0;1.2;0.1;0.;0.])
        mech = Mechanism(origin, bodies, eqs, [contineqcs_left; contineqcs_right], g = g, Δt = Δt, spring=spring, damper=damper)
    end
    return mech
end

function initializehumanoid!(mechanism::Mechanism; tran=[0,0,1.5], rot=[0.1,0,0]) where T
    setPosition!(mechanism, geteqconstraint(mechanism, "auto_generated_floating_joint"), [tran; rot])
    zeroVelocity!(mechanism)
end


