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
        object_type=:none,
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
        # if model_type == :end_effector
            # TODO place the contact points for each finger of the end-effector
        # elseif model_type == :no_end_effector
            # spherical end-effector contact
        # location = [-0.01; 0.004; 0.1125]
        location = [-0.01; 0.004; 0.0625]
        normal = [0.0; 0.0; 1.0]
        contact = contact_constraint(
            get_body(mech, :link7),
            normal,
            friction_coefficient=friction_coefficient,
            contact_origin=location,
            contact_radius=0.00,
            # contact_radius=0.05,
            name=:end_effector)
        push!(contacts, contact)
        # end
    end

    if object_type == :box
        sidex = 0.1
        sidey = 0.1
        sidez = 0.25
        box = Box(sidex, sidey, sidez, 10.0,
            name=:box,
            color=RGBA(1.0, 0.0, 0.0))
        ee = get_body(mech, :link6)


        corners = [
                    [[ sidex / 2.0;  sidey / 2.0; -sidez / 2.0]]
                    [[ sidex / 2.0; -sidey / 2.0; -sidez / 2.0]]
                    [[-sidex / 2.0;  sidey / 2.0; -sidez / 2.0]]
                    [[-sidex / 2.0; -sidey / 2.0; -sidez / 2.0]]
                    [[ sidex / 2.0;  sidey / 2.0;  sidez / 2.0]]
                    [[ sidex / 2.0; -sidey / 2.0;  sidez / 2.0]]
                    [[-sidex / 2.0;  sidey / 2.0;  sidez / 2.0]]
                    [[-sidex / 2.0; -sidey / 2.0;  sidez / 2.0]]
                ]

        normal = [[0.0, 0.0, 1.0] for i = 1:8]
        contact_radius = [0.0 for i = 1:8]
        friction_coefficient = 0.5 * ones(8)

        box_contacts = contact_constraint(box, normal,
            friction_coefficient=friction_coefficient,
            contact_origins=corners,
            contact_radius=contact_radius,
            contact_type=:nonlinear)

        ee_location = SA[0.025; -0.075; 0.01]
        ee_radius = 0.06

        collision = Dojo.SphereBoxCollision{Float64,2,3,6}(
            ee_location,
            sidex,
            sidey,
            sidez,
            ee_radius,
        )

        # collision = Dojo.SphereSphereCollision{Float64,2,3,6}(
        #     SA[0.0; 0.0; 0.0],
        #     SA[0.0; 0.0; 0.3],
        #     0.075,
        #     # SA[0.5 * sidex; 0.0; 0.0],
        #     # SA[0.0; 0.5 * sidey; 0.0],
        #     # SA[0.0; 0.0; 0.5 * sidez],
        #     0.075,
        # )

        # collision = SphereCapsuleCollision{Float64,2,3,6}(
        #     szeros(3),
        #     SA[0.0; 0.0; 0.5],
        #     SA[0.0; 0.0; -0.5],
        #     0.05,
        #     0.05,
        # )

        friction_parameterization = SA{Float64}[
            1.0  0.0
            0.0  1.0
        ]
        body_body_contact = NonlinearContact{Float64,8}(0.1, friction_parameterization, collision)

        object_contacts = [
                ContactConstraint((body_body_contact, ee.id, box.id), name=:body_body1),
                box_contacts...,
        ]

        bodies = [bodies..., box]
        contacts = [contacts..., object_contacts...]

        joints = [joints..., JointConstraint(Floating(origin, box))]
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
