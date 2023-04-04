function get_humanoid(; 
    timestep=0.01, 
    input_scaling=timestep, 
    gravity=-9.81, 
    urdf=:humanoid, 
    springs=0.0, 
    dampers=0.0,
    parse_springs=true, 
    parse_dampers=true,
    limits=false,
    joint_limits=Dict(),
    keep_fixed_joints=false, 
    friction_coefficient=0.8, 
	contact_feet=true, 
    contact_body=false,
    T=Float64)

    # mechanism
    path = joinpath(@__DIR__, "../dependencies/$(string(urdf)).urdf")
    mech = Mechanism(path; floating=true, T,
        gravity, timestep, input_scaling, 
        parse_damper, keep_fixed_joints)

    # springs and dampers
    !parse_springs && set_springs!(mech.joints, springs)
    !parse_dampers && set_dampers!(mech.joints, dampers)

    # joint limits
    if limits
        joints = set_limits(mech, joint_limits)

        mech = Mechanism(Origin{T}(), mech.bodies, joints;
            gravity, timestep, input_scaling)
    end

    # contacts
    origin = Origin{T}()
    bodies = mech.bodies
    joints = mech.joints
    contacts = ContactConstraint{T}[]

    if contact_feet
        # Foot contact
        left_foot = get_body(mech, :left_foot)

		aa = -0.43000 * [-0.44721, 0.00000, 0.89442]
		ql = Dojo.axis_angle_to_quaternion(aa)
        qll = ql * RotX(-1.57080)*RotY(1.47585)*RotZ(-1.47585) # Quaternion(RotXYZ(roll=-1.57080, pitch=1.47585, yaw=-1.47585)) # roll pitch yaw
        qlr = ql * RotX(+1.57080)*RotY(1.47585)*RotZ(+1.47585) # Quaternion(RotXYZ(roll=+1.57080, pitch=1.47585, yaw=+1.47585)) # roll pitch yaw

        pfll = vector_rotate([ 0.5 * left_foot.shape.shapes[1].shapes[1].rh[2] + 0.03500; -0.03; 0.0], qll)
        pbll = vector_rotate([-0.5 * left_foot.shape.shapes[1].shapes[1].rh[2] + 0.03500; -0.03; 0.0], qll)
        pflr = vector_rotate([ 0.5 * left_foot.shape.shapes[1].shapes[1].rh[2] + 0.03500; +0.01; 0.0], qlr)
		pblr = vector_rotate([-0.5 * left_foot.shape.shapes[1].shapes[1].rh[2] + 0.03500; +0.01; 0.0], qlr)
        p = [0.0,0.054,0.]
        o = left_foot.shape.shapes[1].shapes[1].rh[1]
        contacts_points = [
					p,
                    # pfll,
                    # pbll,
                    # pflr,
                    # pblr,
                   ]
        contact_radius = [
                    o,
                    # o,
                    # o,
                    # o,
                  ]
        n = length(contacts_points)
        normal = [Z_AXIS for i = 1:n]
        friction_coefficients = friction_coefficient * ones(T, n)

        contacts_left = contact_constraint(left_foot, normal; 
            friction_coefficient=friction_coefficients, 
            contact_origins=contacts_points, 
            contact_radius)

        right_foot = get_body(mech, :right_foot)

        pfr = [0.5 * right_foot.shape.shapes[1].shapes[1].rh[2]; 0.0; 0.0]
        ofr = right_foot.shape.shapes[1].shapes[1].rh[1]
        pbr = [-0.5 * right_foot.shape.shapes[1].shapes[1].rh[2]; 0.0; 0.0]
        obr = right_foot.shape.shapes[1].shapes[1].rh[1]

        contact_points = [
                    pfr,
                    pbr,
                   ]

        contact_radius = [
                            ofr,
                            obr,
                         ]
                  
        n = length(contact_points)
        normal = [Z_AXIS for i = 1:n]
        friction_coefficients = friction_coefficient * ones(T, n)

        contacts_right = contact_constraint(right_foot, normal; 
            friction_coefficient=friction_coefficients, 
            contact_origins=contact_points, 
            contact_radius)

        contacts = [contacts_left, contacts_right]

        # set_minimal_coordinates!(mech, get_joint(mech, :floating_base), [0.0; 0.0; 1.2; 0.1; 0.0; 0.0])
    end

    mech = Mechanism(origin, bodies, joints, contacts;
        gravity, timestep, input_scaling)

    # construction finished
    return mech
end

function initialize_humanoid!(mechanism::Mechanism{T}; 
    body_position=[0.0, 0.0, 1.5], 
    body_orientation=[0.1, 0.0, 0.0]) where T
    set_minimal_coordinates!(mechanism, 
        get_joint(mechanism, :floating_base), 
        [body_position; body_orientation])
    zero_velocity!(mechanism)
end
