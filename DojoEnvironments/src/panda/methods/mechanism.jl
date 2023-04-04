function get_panda(;
    timestep=0.01,
    input_scaling=timestep, 
    gravity=-9.81,
    urdf=:panda_end_effector,
    springs=0.0,
    dampers=0.0,
    parse_springs=true, 
    parse_dampers=true,
    limits=true,
    joint_limits=Dict([
        (:joint1, [-2.8973, 2.8973]),
        (:joint2, [-1.7628, 1.7628]),
        (:joint3, [-2.8973, 2.8973]),
        (:joint4, [-0.0698, 3.0718]),
        (:joint5, [-2.8973, 2.8973]),
        (:joint6, [-3.7525, 0.0175]),
        (:joint7, [-2.8973, 2.8973]),
        (:jointf1, [-0.00, 0.04]),
        (:jointf2, [-0.00, 0.04])])
    keep_fixed_joints=false, 
    friction_coefficient=0.8,
    contact=false,
    object_type=:none,
    T=Float64)

    # mechanism
    path = joinpath(@__DIR__, "../dependencies/$(string(urdf)).urdf")
    mech = Mechanism(path; floating=false, T,
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

    if contact
        # if urdf == :end_effector
            # TODO place the contact points for each finger of the end-effector
        # elseif urdf == :no_end_effector
            # spherical end-effector contact
        location = [-0.01; 0.004; 0.1125]
        normal = Z_AXIS
        contact = contact_constraint(
            get_body(mech, :link7),
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
        gravity, timestep, input_scaling)

    # set_minimal_state!(mech, szeros(minimal_dimension(mech)))

    # construction finished
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
