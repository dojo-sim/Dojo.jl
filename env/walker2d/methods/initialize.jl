function getwalker2d(; Δt::T=0.01, g::T=-9.81, cf::T=1.9,
    contact::Bool=true,
    contact_body::Bool=true,
    limits::Bool=true,
    spring=0.0,
    damper=0.1,
    joint_limits=[[  0,   0, -45,   0,   0, -45] * π/180,
                  [150, 150,  45, 150, 150,  45] * π/180]) where T

    path = joinpath(@__DIR__, "../deps/walker2d.urdf")
    mech = Mechanism(path, false, T, g=g, Δt=Δt, spring=spring, damper=damper)

    # joint limits
    eqcs = deepcopy(mech.eqconstraints)

    if limits
        thigh = get_joint_constraint(mech, :thigh)
        eqcs[thigh.id] = add_limits(mech, thigh, rot_limits=[SVector{1}(joint_limits[1][1]), SVector{1}(joint_limits[2][1])])

        leg = get_joint_constraint(mech, :leg)
        eqcs[leg.id] = add_limits(mech, leg, rot_limits=[SVector{1}(joint_limits[1][2]), SVector{1}(joint_limits[2][2])])

        foot = get_joint_constraint(mech, :foot)
        eqcs[foot.id] = add_limits(mech, foot, rot_limits=[SVector{1}(joint_limits[1][3]), SVector{1}(joint_limits[2][3])])

        thigh_left = get_joint_constraint(mech, :thigh_left)
        eqcs[thigh_left.id] = add_limits(mech, thigh_left, rot_limits=[SVector{1}(joint_limits[1][4]), SVector{1}(joint_limits[2][4])])

        leg_left = get_joint_constraint(mech, :leg_left)
        eqcs[leg_left.id] = add_limits(mech, leg_left, rot_limits=[SVector{1}(joint_limits[1][5]), SVector{1}(joint_limits[2][5])])

        foot_left = get_joint_constraint(mech, :foot_left)
        eqcs[foot_left.id] = add_limits(mech, foot_left, rot_limits=[SVector{1}(joint_limits[1][6]), SVector{1}(joint_limits[2][6])])

        mech = Mechanism(Origin{T}(), [mech.bodies...], [eqcs...], g=g, Δt=Δt, spring=spring, damper=damper)
    end

    if contact
        origin = Origin{T}()
        bodies = mech.bodies.values
        eqcs = mech.eqconstraints.values

        normal = [0.0; 0.0; 1.0]
        names = contact_body ? getfield.(mech.bodies, :name) : [:ffoot, :foot]
        bounds = []
        for name in names
            body = get_body(mech, name)
            if name in [:foot, :foot_left] # need special case for foot
                # torso
                pf = [0,0, +0.5 * body.shape.rh[2]]
                pb = [0,0, -0.5 * body.shape.rh[2]]
                o = [0;0; body.shape.rh[1]]
                push!(bounds, contact_constraint(body, normal, cf=cf, p=pf, offset=o))
                push!(bounds, contact_constraint(body, normal, cf=cf, p=pb, offset=o))
            else
                p = [0;0; 0.5 * body.shape.rh[2]]
                o = [0;0; body.shape.rh[1]]
                push!(bounds, contact_constraint(body, normal, cf=cf, p=p, offset=o))
            end
        end
        set_position(mech, get_joint_constraint(mech, :floating_joint), [1.25, 0.0, 0.0])
        mech = Mechanism(origin, bodies, eqcs, [bounds...], g=g, Δt=Δt, spring=spring, damper=damper)
    end
    return mech
end

function initializewalker2d!(mechanism::Mechanism; x::T=0.0, z::T=0.0, θ::T=0.0) where {T}
    set_position(mechanism,
                 get_joint_constraint(mechanism, :floating_joint),
                 [z + 1.25 , -x, -θ])
    for eqc in mechanism.eqconstraints
        (eqc.name != :floating_joint) && set_position(mechanism, eqc, zeros(control_dimension(eqc)))
    end
    zeroVelocity!(mechanism)
end
