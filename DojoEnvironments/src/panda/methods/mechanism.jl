function get_panda(;
    timestep=0.01,
    input_scaling=timestep, 
    gravity=-9.81,
    urdf=:panda_end_effector,
    springs=0,
    dampers=5,
    parse_springs=true, 
    parse_dampers=false,
    limits=true,
    joint_limits=Dict([
        (:joint1, [-2.8973, 2.8973]),
        (:joint2, [-1.7628, 1.7628]),
        (:joint3, [-2.8973, 2.8973]),
        (:joint4, [-3.0718,-0.0698]),
        (:joint5, [-2.8973, 2.8973]),
        (:joint6, [-0.0175, 3.7525]),
        (:joint7, [-2.8973, 2.8973]),
        (:jointf1, [-0.00, 0.04]),
        (:jointf2, [-0.00, 0.04])]),
    keep_fixed_joints=false, 
    friction_coefficient=0.8,
    contact=false,
    object_type=:none,
    T=Float64)

    # mechanism
    path = joinpath(@__DIR__, "../dependencies/$(string(urdf)).urdf")
    mechanism = Mechanism(path; floating=false, T,
        gravity, timestep, input_scaling, 
        parse_dampers, keep_fixed_joints)

    # springs and dampers
    !parse_springs && set_springs!(mechanism.joints, springs)
    !parse_dampers && set_dampers!(mechanism.joints, dampers)

    # joint limits    
    if limits
        joints = set_limits(mechanism, joint_limits)

        mechanism = Mechanism(mechanism.origin, mechanism.bodies, joints;
            gravity, timestep, input_scaling)
    end

    # contacts
    origin = mechanism.origin
    bodies = mechanism.bodies
    joints = mechanism.joints
    contacts = ContactConstraint{T}[]

    if contact
        # if urdf == :end_effector
            # TODO place the contact points for each finger of the end-effector
        # elseif urdf == :no_end_effector
            # spherical end-effector contact
        location = [-0.01; 0.004; 0.1125]
        normal = Z_AXIS
        contact = contact_constraint(
            get_body(mechanism, :link7),
            normal,
            friction_coefficient=friction_coefficient,
            contact_origin=location,
            contact_radius=0.05,
            name=:end_effector)
        push!(contacts, contact)
        # end
    end

    if object_type == :box 
        sidex = 0.1
        sidey = 0.1 
        sidez = 0.25
        box = Box(sidex, sidey, sidez, 10, 
            name=:box,
            color=RGBA(1, 0, 0))
        ee = get_body(mechanism, :link6)

        
        corners = [
                    [[ sidex / 2;  sidey / 2; -sidez / 2]]
                    [[ sidex / 2; -sidey / 2; -sidez / 2]]
                    [[-sidex / 2;  sidey / 2; -sidez / 2]]
                    [[-sidex / 2; -sidey / 2; -sidez / 2]]
                    [[ sidex / 2;  sidey / 2;  sidez / 2]]
                    [[ sidex / 2; -sidey / 2;  sidez / 2]]
                    [[-sidex / 2;  sidey / 2;  sidez / 2]]
                    [[-sidex / 2; -sidey / 2;  sidez / 2]]
                ]
               
        normal = [Z_AXIS for i = 1:8]
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
        #     SA[0; 0; 0],
        #     SA[0; 0; 0.3],
        #     0.075,
        #     # SA[0.5 * sidex; 0; 0],
        #     # SA[0; 0.5 * sidey; 0],
        #     # SA[0; 0; 0.5 * sidez],
        #     0.075,
        # )

        # collision = SphereCapsuleCollision{Float64,2,3,6}(
        #     szeros(3),
        #     SA[0; 0; 0.5],
        #     SA[0; 0; -0.5],
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

    mechanism = Mechanism(origin, bodies, joints, contacts;
        gravity, timestep, input_scaling)

    # zero configuration
    initialize_panda!(mechanism)

    # construction finished
    return mechanism
end

function initialize_panda!(mechanism::Mechanism;
    joint_angles=[0;0.5;0;-0.5;0;0.5;0; zeros(input_dimension(mechanism)-7)],
    joint_velocities=zeros(input_dimension(mechanism)))

    zero_velocity!(mechanism)
    zero_coordinates!(mechanism)

    nu = input_dimension(mechanism)
    x = vcat([[joint_angles[i], joint_velocities[i]] for i=1:nu]...)
    set_minimal_state!(mechanism, x)

    return
end
