function get_panda(;
        timestep=0.01,
        gravity=-9.81,
        friction_coefficient=0.8,
        spring=0.0,
        damper=0.0,
        parse_damper=true,
        contact=false,
        limits=true,
        model_type=:end_effector,
        # joint_limits=[[-0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25],
        #               [ 0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25]],
        # joint_limits=[[-2.9671, -1.8326, -2.9671, -0.4000, -2.9671, -3.8223, -2.9671],
        #               [ 2.9671,  1.8326,  2.9671,  3.1416,  2.9671,  0.0873,  2.9671]],
        joint_limits=[[-2.8973, -1.7628, -2.8973, -0.0698, -2.8973, -3.7525, -2.8973, -0.00],
                      [ 2.8973,  1.7628,  2.8973,  3.0718,  2.8973,  0.0175,  2.8973,  0.04]],
        T=Float64)

    path = joinpath(@__DIR__, "../deps/panda_$(String(model_type)).urdf")
    mech = Mechanism(path; floating=false, T,
        gravity,
        timestep,
        parse_damper)

    # Adding springs and dampers
    set_springs!(mech.joints, spring)
    set_dampers!(mech.joints, damper)

    # joint limits
    joints = deepcopy(mech.joints)

    if limits
        for (i,joint) in enumerate(joints)
            id = joint.id
            if input_dimension(joint.translational) == 0 && input_dimension(joint.rotational) == 1
                joints[i] = add_limits(mech, joint,
                    rot_limits=[SVector{1}(joint_limits[1][id-1]), SVector{1}(joint_limits[2][id-1])])
            end
            if input_dimension(joint.translational) == 1 && input_dimension(joint.rotational) == 0
                joints[i] = add_limits(mech, joint,
                    tra_limits=[SVector{1}(0.00), SVector{1}(0.04)])
            end
        end
        mech = Mechanism(Origin{T}(), [mech.bodies...], [joints...];
            gravity,
            timestep)
    end


    origin = Origin{T}()
    bodies = mech.bodies
    joints = mech.joints
    contacts = ContactConstraint{T}[]

    if contact
        if model_type == :end_effector
            # TODO place the contact points for each finger of the end-effector
        elseif model_type == :no_end_effector
            # spherical end-effector contact
            location = [-0.01; 0.004; 0.01]
            normal = [0.0; 0.0; 1.0]
            offset = [0.0; 0.0; 0.05]
            contact = contact_constraint(
                get_body(mech, :link7),
                normal;
                friction_coefficient,
                contact_point=location,
                offset,
                name=:end_effector)
            push!(contacts, contact)
        end
    end


    mech = Mechanism(origin, bodies, joints, contacts;
        gravity,
        timestep)

    set_minimal_state!(mech, szeros(minimal_dimension(mech)))
    return mech
end

function initialize_panda!(mechanism::Mechanism{T};
    joint_angles=[[0,-0.8,0,1.6,0,-2.4,0]; zeros(input_dimension(mechanism)-7)],
    joint_velocities=zeros(input_dimension(mechanism))) where T

    nu = input_dimension(mechanism)
    zero_velocity!(mechanism)
    y = vcat([[joint_angles[i], joint_velocities[i]] for i=1:nu]...)
    set_minimal_state!(mechanism, y)
    return nothing
end
