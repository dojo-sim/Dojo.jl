function get_halfcheetah(; timestep::T=0.01, gravity=[0.0; 0.0; -9.81], friction_coefficient::T=0.4,
    contact::Bool=true,
    contact_body::Bool=true,
    limits::Bool = true,
    spring=[240, 180, 120, 180, 120, 60.],
    damper=[6., 4.5, 3., 4.5, 3., 1.5],
    joint_limits=[[-0.52, -0.785, -0.400, -1.0, -1.20, -0.5],
                  [ 1.05,  0.785,  0.785,  0.7,  0.87,  0.5]]) where T

    path = joinpath(@__DIR__, "../deps/halfcheetah.urdf")
    mech = Mechanism(path, false, T, gravity=gravity, timestep=timestep, spring=spring, damper=damper)

    # joint limits
    joints = deepcopy(mech.joints)

    if limits
        bthigh = get_joint_constraint(mech, :bthigh)
        joints[bthigh.id] = add_limits(mech, bthigh, rot_limits=[SVector{1}(joint_limits[1][1]), SVector{1}(joint_limits[2][1])])

        bshin = get_joint_constraint(mech, :bshin)
        joints[bshin.id] = add_limits(mech, bshin, rot_limits=[SVector{1}(joint_limits[1][2]), SVector{1}(joint_limits[2][2])])

        bfoot = get_joint_constraint(mech, :bfoot)
        joints[bfoot.id] = add_limits(mech, bfoot, rot_limits=[SVector{1}(joint_limits[1][3]), SVector{1}(joint_limits[2][3])])

        fthigh = get_joint_constraint(mech, :fthigh)
        joints[fthigh.id] = add_limits(mech, fthigh, rot_limits=[SVector{1}(joint_limits[1][4]), SVector{1}(joint_limits[2][4])])

        fshin = get_joint_constraint(mech, :fshin)
        joints[fshin.id] = add_limits(mech, fshin, rot_limits=[SVector{1}(joint_limits[1][5]), SVector{1}(joint_limits[2][5])])

        ffoot = get_joint_constraint(mech, :ffoot)
        joints[ffoot.id] = add_limits(mech, ffoot, rot_limits=[SVector{1}(joint_limits[1][6]), SVector{1}(joint_limits[2][6])])

        mech = Mechanism(Origin{T}(), [mech.bodies...], [joints...], gravity=gravity, timestep=timestep, spring=spring, damper=damper)
    end

    if contact
        origin = Origin{T}()
        bodies = mech.bodies
        joints = mech.joints

        normal = [0.0; 0.0; 1.0]
        names = contact_body ? getfield.(mech.bodies, :name) : [:ffoot, :bfoot]
        bounds = []
        for name in names
            body = get_body(mech, name)
            if name == :torso # need special case for torso
                # torso
                pf = [+0.5 * body.shape.shape[1].rh[2]; 0.0; 0.0]
                pb = [-0.5 * body.shape.shape[1].rh[2]; 0.0; 0.0]
                o = [0;0; body.shape.shape[1].rh[1]]
                push!(bounds, contact_constraint(body, normal, friction_coefficient=friction_coefficient, contact_point=pf, offset=o))
                push!(bounds, contact_constraint(body, normal, friction_coefficient=friction_coefficient, contact_point=pb, offset=o))

                # head
                pf = [+0.5 * body.shape.shape[1].rh[2] + 0.214; 0.0; 0.1935]
                o = [0;0; body.shape.shape[2].rh[1]]
                push!(bounds, contact_constraint(body, normal, friction_coefficient=friction_coefficient, contact_point=pf, offset=o))
            else
                p = [0;0; -0.5 * body.shape.rh[2]]
                o = [0;0; body.shape.rh[1]]
                push!(bounds, contact_constraint(body, normal, friction_coefficient=friction_coefficient, contact_point=p, offset=o))
            end
        end
        set_position!(mech, get_joint_constraint(mech, :floating_joint), [0.576509, 0.0, 0.02792])
        mech = Mechanism(origin, bodies, joints, [bounds...], gravity=gravity, timestep=timestep, spring=spring, damper=damper)
    end
    return mech
end

function initialize_halfcheetah!(mechanism::Mechanism; x::T=0.0, z::T=0.0, θ::T=0.0) where T
    set_position!(mechanism,
                 get_joint_constraint(mechanism, :floating_joint),
                 [z + 0.576509, -x, -θ + 0.02792])
    for joint in mechanism.joints
        (joint.name != :floating_joint) && set_position!(mechanism, joint, zeros(control_dimension(joint)))
    end
    zero_velocity!(mechanism)
end
