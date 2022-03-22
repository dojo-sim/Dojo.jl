function get_panda(;
        timestep=0.01,
        gravity=-9.81,
        friction_coefficient=0.8,
        spring=0.0,
        damper=0.01, # this value comes from the official URDF https://github.com/frankaemika/franka_ros/blob/develop/franka_gazebo/test/launch/panda-gazebo.urdf
        contact=false,
        limits=true,
        joint_limits=[[-2.9671, -1.8326, -2.9671, -0.25, -0.25, -0.25, -0.25],
                      [ 2.9671,  1.8326,  2.9671,  0.25,  0.25,  0.25,  0.25]],
        T=Float64)

    path = joinpath(@__DIR__, "../deps/panda.urdf")
    mech = Mechanism(path, false, T,
        gravity=gravity,
        timestep=timestep,
        spring=spring,
        damper=damper)

    # Adding springs and dampers
    for joint in mech.joints
        joint.damper = true
        joint.spring = true
        joint.translational.spring=spring
        joint.translational.damper=damper
        joint.rotational.spring=spring
        joint.rotational.damper=damper
    end

    # joint limits
    joints = deepcopy(mech.joints)
    if limits
        for (i,joint) in enumerate(joints)
            id = joint.id
            if input_dimension(joint) > 0
                joints[i] = add_limits(mech, joint,
                    rot_limits=[SVector{1}(joint_limits[1][id-1]), SVector{1}(joint_limits[2][id-1])])
            end
        end
        mech = Mechanism(Origin{T}(), [mech.bodies...], [joints...],
            gravity=gravity,
            timestep=timestep,
            spring=spring,
            damper=damper)
    end


    origin = Origin{T}()
    bodies = mech.bodies
    joints = mech.joints
    contacts = ContactConstraint{T}[]

    if contact
        # spherical end-effector contact
        location1 = [-0.01; 0.004; 0.01]
        normal = [0.0; 0.0; 1.0]
        offset = [0.0; 0.0; 0.05]
        contact = contact_constraint(
            get_body(mech, :link7),
            normal,
            friction_coefficient=friction_coefficient,
            contact_point=location1,
            offset=offset,
            name=:end_effector)
        push!(contacts, contact)

        # spherical end-effector contact
        location = location1 + [0.00,0.,0.03]
        normal = [0.0; 0.0; 1.0]
        offset = [0.0; 0.0; 0.01]
        contact = contact_constraint(
            get_body(mech, :link7),
            normal,
            friction_coefficient=friction_coefficient,
            contact_point=location,
            offset=offset,
            name=:end_effector)
        push!(contacts, contact)

        # spherical end-effector contact
        location = location1 + [0,0,0]
        normal = [0.0; 0.0; 1.0]
        offset = [0.0; 0.0; 0.01]
        contact = contact_constraint(
            get_body(mech, :link7),
            normal,
            friction_coefficient=friction_coefficient,
            contact_point=location,
            offset=offset,
            name=:end_effector)
        push!(contacts, contact)

        # spherical end-effector contact
        location = location1 + [0,0,0]
        normal = [0.0; 0.0; 1.0]
        offset = [0.0; 0.0; 0.01]
        contact = contact_constraint(
            get_body(mech, :link7),
            normal,
            friction_coefficient=friction_coefficient,
            contact_point=location,
            offset=offset,
            name=:end_effector)
        push!(contacts, contact)
    end


    mech = Mechanism(origin, bodies, joints, contacts,
        gravity=gravity,
        timestep=timestep,
        spring=spring,
        damper=damper)

    set_minimal_state!(mech, szeros(minimal_dimension(mech)))
    return mech
end

function initialize_panda!(mechanism::Mechanism{T};
    joint_angles=zeros(7),
    joint_velocities=zeros(7)) where T

    zero_velocity!(mechanism)
    y = vcat([[joint_angles[i], joint_velocities[i]] for i=1:7]...)
    set_minimal_state!(mechanism, y)
    return nothing
end
