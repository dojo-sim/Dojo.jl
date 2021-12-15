function gethalfcheetah(; Δt::T=0.01, g::T=-9.81, cf::T=0.8,
    contact::Bool=true,
    contact_body::Bool=true,
    limits::Bool = true,
    spring=[240, 180, 120, 180, 120, 60.],
    damper=[6., 4.5, 3., 4.5, 3., 1.5],
    joint_limits=[[-0.52, -0.785, -0.400, -1.0, -1.20, -0.5],
                  [ 1.05,  0.785,  0.785,  0.7,  0.87,  0.5]]) where T

    path = joinpath(@__DIR__, "../deps/halfcheetah.urdf")
    mech = Mechanism(path, false, T, g=g, Δt=Δt, spring=spring, damper=damper)

    # joint limits
    eqcs = deepcopy(mech.eqconstraints)

    if limits
        bthigh = geteqconstraint(mech, "bthigh")
        eqcs[bthigh.id] = add_limits(mech, bthigh, rot_limits=[SVector{1}(joint_limits[1][1]), SVector{1}(joint_limits[2][1])])

        bshin = geteqconstraint(mech, "bshin")
        eqcs[bshin.id] = add_limits(mech, bshin, rot_limits=[SVector{1}(joint_limits[1][2]), SVector{1}(joint_limits[2][2])])

        bfoot = geteqconstraint(mech, "bfoot")
        eqcs[bfoot.id] = add_limits(mech, bfoot, rot_limits=[SVector{1}(joint_limits[1][3]), SVector{1}(joint_limits[2][3])])

        fthigh = geteqconstraint(mech, "fthigh")
        eqcs[fthigh.id] = add_limits(mech, fthigh, rot_limits=[SVector{1}(joint_limits[1][4]), SVector{1}(joint_limits[2][4])])

        fshin = geteqconstraint(mech, "fshin")
        eqcs[fshin.id] = add_limits(mech, fshin, rot_limits=[SVector{1}(joint_limits[1][5]), SVector{1}(joint_limits[2][5])])

        ffoot = geteqconstraint(mech, "ffoot")
        eqcs[ffoot.id] = add_limits(mech, ffoot, rot_limits=[SVector{1}(joint_limits[1][6]), SVector{1}(joint_limits[2][6])])

        mech = Mechanism(Origin{T}(), [mech.bodies...], [eqcs...], g=g, Δt=Δt, spring=spring, damper=damper)
    end

    if contact
        origin = Origin{T}()
        bodies = mech.bodies.values
        eqcs = mech.eqconstraints.values

        normal = [0.0; 0.0; 1.0]
        names = contact_body ? getfield.(mech.bodies, :name) : ["ffoot", "bfoot"]
        bounds = []
        for name in names
            body = getbody(mech, name)
            if name == "torso" # need special case for torso
                # torso
                pf = [+0.5 * body.shape.shape[1].rh[2]; 0.0; 0.0]
                pb = [-0.5 * body.shape.shape[1].rh[2]; 0.0; 0.0]
                o = [0;0; body.shape.shape[1].rh[1]]
                push!(bounds, contactconstraint(body, normal, cf, p=pf, offset=o))
                push!(bounds, contactconstraint(body, normal, cf, p=pb, offset=o))

                # head
                pf = [+0.5 * body.shape.shape[1].rh[2] + 0.214; 0.0; 0.1935]
                o = [0;0; body.shape.shape[2].rh[1]]
                push!(bounds, contactconstraint(body, normal, cf, p=pf, offset=o))
            else
                p = [0;0; -0.5 * body.shape.rh[2]]
                o = [0;0; body.shape.rh[1]]
                push!(bounds, contactconstraint(body, normal, cf, p=p, offset=o))
            end
        end
        setPosition!(mech, geteqconstraint(mech, "floating_joint"), [0.576509, 0.0, 0.02792])
        mech = Mechanism(origin, bodies, eqcs, [bounds...], g=g, Δt=Δt, spring=spring, damper=damper)
    end
    return mech
end

function initializehalfcheetah!(mechanism::Mechanism; x::T=0.0, z::T=0.0, θ::T=0.0) where {T}
    setPosition!(mechanism,
                 geteqconstraint(mechanism, "floating_joint"),
                 [z + 0.576509, -x, -θ + 0.02792])
    for eqc in mechanism.eqconstraints
        (eqc.name != "floating_joint") && setPosition!(mechanism, eqc, zeros(controldim(eqc)))
    end
    zeroVelocity!(mechanism)
end
