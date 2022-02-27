function get_ant(; 
    timestep=0.05, 
    gravity=[0.0; 0.0; -9.81], 
    friction_coefficient=0.5,
    spring=0.0, 
    damper=1.0, 
    contact_feet=true, 
    contact_body=true,
    limits::Bool=true,
    joint_limits=[[-30.0, 30.0, -30.0, -70.0, -30.0, -70.0, -30.0, 30.0] * π / 180.0,
                  [ 30.0, 70.0,  30.0, -30.0,  30.0, -30.0,  30.0, 70.0] * π / 180.0],
    T=Float64)

    path = joinpath(@__DIR__, "../deps/ant.urdf")
    mech = Mechanism(path, true, T, 
        gravity=gravity, 
        timestep=timestep, 
        spring=spring, 
        damper=damper)

    # joint limits
    joints = deepcopy(mech.joints)

    if limits
        hiparent_vertex = get_joint_constraint(mech, :hip_1)
        joints[hiparent_vertex.id] = add_limits(mech, hiparent_vertex, 
            rot_limits=[SVector{1}(joint_limits[1][1]), SVector{1}(joint_limits[2][1])])

        ankle1 = get_joint_constraint(mech, :ankle_1)
        joints[ankle1.id] = add_limits(mech, ankle1, 
            rot_limits=[SVector{1}(joint_limits[1][2]), SVector{1}(joint_limits[2][2])])

        hip2 = get_joint_constraint(mech, :hip_2)
        joints[hip2.id] = add_limits(mech, hip2, 
            rot_limits=[SVector{1}(joint_limits[1][3]), SVector{1}(joint_limits[2][3])])

        ankle2 = get_joint_constraint(mech, :ankle_2)
        joints[ankle2.id] = add_limits(mech, ankle2, 
            rot_limits=[SVector{1}(joint_limits[1][4]), SVector{1}(joint_limits[2][4])])

        hip3 = get_joint_constraint(mech, :hip_3)
        joints[hip3.id] = add_limits(mech, hip3, 
            rot_limits=[SVector{1}(joint_limits[1][5]), SVector{1}(joint_limits[2][5])])

        ankle3 = get_joint_constraint(mech, :ankle_3)
        joints[ankle3.id] = add_limits(mech, ankle3, 
            rot_limits=[SVector{1}(joint_limits[1][6]), SVector{1}(joint_limits[2][6])])

        hip4 = get_joint_constraint(mech, :hip_4)
        joints[hip4.id] = add_limits(mech, hip4, 
            rot_limits=[SVector{1}(joint_limits[1][7]), SVector{1}(joint_limits[2][7])])

        ankle4 = get_joint_constraint(mech, :ankle_4)
        joints[ankle4.id] = add_limits(mech, ankle4, 
            rot_limits=[SVector{1}(joint_limits[1][8]), SVector{1}(joint_limits[2][8])])

        mech = Mechanism(Origin{T}(), mech.bodies, [joints...], 
            gravity=gravity, 
            timestep=timestep, 
            spring=spring, 
            damper=damper)
    end

    if contact_feet
        origin = Origin{T}()
        bodies = mech.bodies
        joints = mech.joints

        # foot contact
        normal = [0.0; 0.0; 1.0]
        foot_names = [:front_left_foot, :front_right_foot, :left_back_foot, :right_back_foot]
        foot = [get_body(mech, name) for name in foot_names]
        p = [[0.2; 0.2; 0.0], [-0.2; 0.2; 0.0], [-0.2; -0.2; 0.0], [0.2; -0.2; 0.0]]
        o = [[0.0; 0.0; f.shape.rh[1]] for f in foot]
        contacts = [contact_constraint(foot[i], normal, friction_coefficient=friction_coefficient, contact_point=p[i], offset=o[i]) for i = 1:length(foot_names)]

        if contact_body
            # torso contact
            torso = get_body(mech, :torso)
            p = [0.0; 0.0; 0.0]
            o = [0.0; 0.0; torso.shape.r]
            torso_contacts = contact_constraint(torso, normal, friction_coefficient=friction_coefficient, contact_point=p, offset=o)

            # elbow contact
            elbow_names = [:aux_1, :aux_2, :aux_3, :aux_4]
            elbow = [get_body(mech, e) for e in elbow_names]
            p = [-[0.1; 0.1; 0.0], -[-0.1; 0.1; 0.0], -[-0.1;-0.1; 0.0], -[0.1; -0.1; 0.0]]
            o = [[0.0; 0.0; e.shape.rh[1]] for e in elbow]
            elbow_contacts = [contact_constraint(elbow[i], normal, friction_coefficient=friction_coefficient, contact_point=p[i], offset=o[i]) for i = 1:length(elbow_names)]

            contacts = [contacts..., torso_contacts, elbow_contacts...]
        end
        mech = Mechanism(origin, bodies, joints, contacts, 
        gravity=gravity, 
        timestep=timestep, 
        spring=spring, 
        damper=damper)
    end

    return mech
end

function initialize_ant!(mechanism::Mechanism; 
        ankle=0.25, 
        alt=0.15, 
        pos=[0.0; 0.0; 0.48 + alt], 
        rot=[0.0; 0.0; 0.00 * π]) where T
        
    set_minimal_coordinates!(mechanism, get_joint_constraint(mechanism, :auto_generated_floating_joint), [pos; rot])

    for i in [1,4]
        set_minimal_coordinates!(mechanism, get_joint_constraint(mechanism, Symbol("hip_$i")), [0.0 * π])
        set_minimal_coordinates!(mechanism, get_joint_constraint(mechanism, Symbol("ankle_$i")), [ankle * π])
    end

    for i in [2,3]
        set_minimal_coordinates!(mechanism, get_joint_constraint(mechanism, Symbol("hip_$i")), [0.0 * π])
        set_minimal_coordinates!(mechanism, get_joint_constraint(mechanism, Symbol("ankle_$i")), [-ankle * π])
    end

    zero_velocity!(mechanism)
end
