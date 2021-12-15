function gethopper(; Δt::T=0.01, g::T=-9.81, cf::T=2.0,
    contact::Bool=true,
    contact_body::Bool=true,
    limits::Bool = true,
    spring=0.0,
    damper=1.0,
    joint_limits=[[  0,   0, -45] * π/180,
                  [150, 150,  45] * π/180]) where T

    path = joinpath(@__DIR__, "../deps/hopper.urdf")
    mech = Mechanism(path, false, T, g=g, Δt=Δt, spring=spring, damper=damper)

    # joint limits
    eqcs = deepcopy(mech.eqconstraints)

    if limits
        thigh = geteqconstraint(mech, "thigh")
        eqcs[thigh.id] = add_limits(mech, thigh, rot_limits=[SVector{1}(joint_limits[1][1]), SVector{1}(joint_limits[2][1])])

        leg = geteqconstraint(mech, "leg")
        eqcs[leg.id] = add_limits(mech, leg, rot_limits=[SVector{1}(joint_limits[1][2]), SVector{1}(joint_limits[2][2])])

        foot = geteqconstraint(mech, "foot")
        eqcs[foot.id] = add_limits(mech, foot, rot_limits=[SVector{1}(joint_limits[1][3]), SVector{1}(joint_limits[2][3])])

        mech = Mechanism(Origin{T}(), [mech.bodies...], [eqcs...], g=g, Δt=Δt, spring=spring, damper=damper)
    end

    if contact
        origin = Origin{T}()
        bodies = mech.bodies.values
        eqcs = mech.eqconstraints.values

        normal = [0.0; 0.0; 1.0]
        names = contact_body ? getfield.(mech.bodies, :name) : ["ffoot", "foot"]
        bounds = []
        for name in names
            body = getbody(mech, name)
            if name == "foot" # need special case for foot
                # torso
                pf = [0,0, +0.5 * body.shape.rh[2]]
                pb = [0,0, -0.5 * body.shape.rh[2]]
                o = [0;0; body.shape.rh[1]]
                push!(bounds, contactconstraint(body, normal, cf, p=pf, offset=o))
                push!(bounds, contactconstraint(body, normal, cf, p=pb, offset=o))
            else
                p = [0;0; 0.5 * body.shape.rh[2]]
                o = [0;0; body.shape.rh[1]]
                push!(bounds, contactconstraint(body, normal, cf, p=p, offset=o))
            end
        end
        setPosition!(mech, geteqconstraint(mech, "floating_joint"), [0.576509, 0.0, 0.02792])
        mech = Mechanism(origin, bodies, eqcs, [bounds...], g=g, Δt=Δt, spring=spring, damper=damper)
    end
    return mech
end

function initializehopper!(mechanism::Mechanism; x::T=0.0, z::T=0.0, θ::T=0.0) where {T}
    setPosition!(mechanism,
                 geteqconstraint(mechanism, "floating_joint"),
                 [z + 1.25 , -x, -θ])
    for eqc in mechanism.eqconstraints
        (eqc.name != "floating_joint") && setPosition!(mechanism, eqc, zeros(controldim(eqc)))
    end
    zeroVelocity!(mechanism)
end
