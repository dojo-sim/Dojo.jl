function get_humanoid(; 
    timestep=0.01, 
    gravity=[0.0; 0.0; -9.81], 
    friction_coefficient=0.8, 
    spring=0.0, 
    damper=0.0,
	contact_feet=true, 
    contact_body=false,
    T=Float64)

    path = joinpath(@__DIR__, "../deps/humanoid.urdf")
    mech = Mechanism(path, true, T, 
        gravity=gravity, 
        timestep=timestep, 
        spring=spring, 
        damper=damper)

    if contact_feet
        origin = Origin{T}()
        bodies = mech.bodies
        eqs = mech.joints

        # Foot contact
        left_foot = get_body(mech, :left_foot)

		aa = -0.43000 * [-0.44721, 0.00000, 0.89442]
		ql = axis_angle_to_quaternion(aa)
        qll = ql * UnitQuaternion(RotXYZ(roll=-1.57080, pitch=1.47585, yaw=-1.47585)) # roll pitch yaw
        qlr = ql * UnitQuaternion(RotXYZ(roll=+1.57080, pitch=1.47585, yaw=+1.47585)) # roll pitch yaw

        pfll = vector_rotate([ 0.5 * left_foot.shape.shape[1].rh[2] + 0.03500; -0.03; 0.0], qll)
        pbll = vector_rotate([-0.5 * left_foot.shape.shape[1].rh[2] + 0.03500; -0.03; 0.0], qll)
        pflr = vector_rotate([ 0.5 * left_foot.shape.shape[1].rh[2] + 0.03500; +0.01; 0.0], qlr)
		pblr = vector_rotate([-0.5 * left_foot.shape.shape[1].rh[2] + 0.03500; +0.01; 0.0], qlr)
        p = [0.0,0.054,0.]
        o = [0.0; 0.0; left_foot.shape.shape[1].rh[1]]
        contacts = [
					p,
                    # pfll,
                    # pbll,
                    # pflr,
                    # pblr,
                   ]
        offsets = [
                    o,
                    # o,
                    # o,
                    # o,
                  ]
        n = length(contacts)
        normal = [[0;0;1.0] for i = 1:n]
        friction_coefficients = friction_coefficient * ones(T, n)

        contacts_left = contact_constraint(left_foot, normal, 
            friction_coefficient=friction_coefficients, 
            contact_points=contacts, offset=offsets)

        right_foot = get_body(mech, :right_foot)

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
        friction_coefficients = friction_coefficient * ones(T, n)

        contacts_right = contact_constraint(right_foot, normal, 
            friction_coefficient=friction_coefficients, 
            contact_points=contacts, 
            offset=offsets)

        set_minimal_coordinates!(mech, get_joint_constraint(mech, :auto_generated_floating_joint), [0;0;1.2;0.1;0.;0.])
        # mech = Mechanism(origin, bodies, eqs, [contacts_left; contacts_right], 
            # gravity=gravity, 
            # timestep=timestep, 
            # spring=spring, 
            # damper=damper)
        mech = Mechanism(origin, bodies, eqs, [contacts_left; ], 
            gravity=gravity, 
            timestep=timestep, 
            spring=spring, 
            damper=damper)
    end
    return mech
end

function initialize_humanoid!(mechanism::Mechanism{T}; 
    body_position=[0.0, 0.0, 1.5], 
    body_orientation=[0.1, 0.0, 0.0]) where T
    set_minimal_coordinates!(mechanism, 
        get_joint_constraint(mechanism, :auto_generated_floating_joint), 
        [body_position; body_orientation])
    zero_velocity!(mechanism)
end
