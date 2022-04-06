function get_hopper(; 
    timestep=0.01, 
    gravity=[0.0; 0.0; -9.81], 
    friction_coefficient=2.0,
    contact_foot=true,
    contact_body=true,
    limits=false,
    spring=0.0,
    damper=1.0,
    joint_limits=[[  0.0,   0.0, -45.0] * π / 180.0,
                  [150.0, 150.0,  45.0] * π / 180.0],
    T=Float64)

    path = joinpath(@__DIR__, "../deps/hopper_good.urdf")
    mech = Mechanism(path, false, T, 
        gravity=gravity, 
        timestep=timestep, 
        spring=spring, 
        damper=damper)

    # joint limits
    joints = deepcopy(mech.joints)

    if limits
        thigh = get_joint(mech, :thigh)
        joints[thigh.id] = add_limits(mech, thigh, 
            rot_limits=[SVector{1}(joint_limits[1][1]), SVector{1}(joint_limits[2][1])])

        @warn "uncomment limits"
        # leg = get_joint(mech, "leg")
        # joints[leg.id] = add_limits(mech, leg, 
            # rot_limits=[SVector{1}(joint_limits[1][2]), SVector{1}(joint_limits[2][2])])
        #
        # foot = get_joint(mech, "foot")
        # joints[foot.id] = add_limits(mech, foot, 
            # rot_limits=[SVector{1}(joint_limits[1][3]), SVector{1}(joint_limits[2][3])])

        mech = Mechanism(Origin{T}(), [mech.bodies...], [joints...], 
            gravity=gravity, 
            timestep=timestep, 
            spring=spring, 
            damper=damper)
    end

    if contact_foot
        origin = Origin{T}()
        bodies = mech.bodies
        joints = mech.joints

        normal = [0.0; 0.0; 1.0]
        names = contact_body ? getfield.(mech.bodies, :name) : [:ffoot, :foot]
        models = []
        for name in names
            body = get_body(mech, name)
            if name == :foot # need special case for foot
                # torso
                pf = [0.0, 0.0, +0.5 * body.shape.rh[2]]
                pb = [0.0, 0.0, -0.5 * body.shape.rh[2]]
                o = body.shape.rh[1]
                push!(models, contact_constraint(body, normal, 
                    friction_coefficient=friction_coefficient, 
                    contact_origin=pf, 
                    contact_radius=o))
                push!(models, contact_constraint(body, normal, 
                    friction_coefficient=friction_coefficient, 
                    contact_origin=pb, 
                    contact_radius=o))
            else
                p = [0.0; 0.0; 0.5 * body.shape.rh[2]]
                o = body.shape.rh[1]
                push!(models, contact_constraint(body, normal, 
                    friction_coefficient=friction_coefficient, 
                    contact_origin=p, 
                    contact_radius=o))
            end
        end
        set_minimal_coordinates!(mech, get_joint(mech, :floating_joint), [1.25, 0.0, 0.0])
        mech = Mechanism(origin, bodies, joints, [models...], 
            gravity=gravity, 
            timestep=timestep, 
            spring=spring, 
            damper=damper)
    end
    return mech
end

function initialize_hopper!(mechanism::Mechanism{T}; 
    body_position=[0.0, 0.0],
    body_orientation=0.0) where T
    #TODO add leg length

    set_minimal_coordinates!(mechanism,
                 get_joint(mechanism, :floating_joint),
                 [body_position[2] + 1.25 , -body_position[1], -body_orientation])
    for joint in mechanism.joints
        (joint.name != :floating_joint) && set_minimal_coordinates!(mechanism, joint, zeros(input_dimension(joint)))
    end
    zero_velocity!(mechanism)
end
