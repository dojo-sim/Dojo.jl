function get_raiberthopper(; 
    timestep=0.05, 
    input_scaling=timestep, 
    gravity=-9.81, 
    body_mass=4.18,
    foot_mass=0.52,
    body_radius=0.1,
    foot_radius=0.05,
    color=RGBA(0.9, 0.9, 0.9),
    springs=[0;0], 
    dampers=[0;0.1],
    limits=false,
    joint_limits=Dict(), 
    contact_foot=true, 
    contact_body=true,
    T=Float64)

    # mechanism
    origin = Origin{Float64}()

    body = Sphere(body_radius, body_mass; color)
    foot = Sphere(foot_radius, foot_mass; color)
    bodies = [body, foot]

    joint_origin_body = JointConstraint(Floating(origin, body))
    joint_body_foot = JointConstraint(Prismatic(body, foot, Z_AXIS;
        parent_vertex=szeros(Float64, 3), child_vertex=szeros(Float64, 3)))
    joints = [joint_origin_body, joint_body_foot]

    mechanism = Mechanism(origin, bodies, joints;
        gravity, timestep, input_scaling)

    # springs and dampers
    set_springs!(mechanism.joints, springs)
    set_dampers!(mechanism.joints, dampers)

    # joint limits    
    if limits
        joints = set_limits(mechanism, joint_limits)

        mechanism = Mechanism(Origin{T}(), mechanism.bodies, joints;
            gravity, timestep, input_scaling)
    end

    # contacts
    origin = Origin{T}()
    bodies = mechanism.bodies
    joints = mechanism.joints
    contacts = ContactConstraint{T}[]

    if contact_foot
         # Contact
        contact_normal = Z_AXIS
        friction_coefficient = 0.5

        # foot
        foot_contacts = contact_constraint(foot, contact_normal; 
            friction_coefficient,
            contact_origin=[0; 0; 0], 
            contact_radius=foot_radius)

        contacts = [foot_contacts]

        # body
        if contact_body
            body_contacts = contact_constraint(body, contact_normal; 
                friction_coefficient,
                contact_origin=[0; 0; 0], 
                contact_radius=body_radius)
            push!(contacts, body_contacts)
        end
    end

    mechanism = Mechanism(origin, bodies, joints, contacts; 
        gravity, timestep, input_scaling)

    # zero configuration
    zero_coordinates!(mechanism)
    set_minimal_coordinates!(mechanism, joint_origin_body, [0; 0; body_radius; 0; 0; 0])
    set_minimal_coordinates!(mechanism, joint_body_foot, [2*body_radius+foot_radius])
    
    # construction finished
    return mechanism
end

function initialize_raiberthopper!(mech::Mechanism{T,Nn,Ne,Nb}; 
    body_position=[0, 0, 0.05],
    leg_length=0.5, 
    body_linear_velocity=zeros(3), 
    body_angular_velocity=zeros(3)) where {T,Nn,Ne,Nb}

    pbody = mechanism.bodies[1]
    cbody = mechanism.bodies[2]
    joint2 = mechanism.joints[2]
    tra2 = joint2.translational

    # origin to body
    set_maximal_configurations!(mechanism.origin, pbody, 
        Δx=[body_position[1:2]; leg_length + body_position[3]])
    set_maximal_velocities!(pbody, 
        v=body_linear_velocity, 
        ω=body_angular_velocity)

    # body to foot
    set_maximal_configurations!(pbody, cbody, 
        Δx=[0; 0; -leg_length], 
        Δq=RotX(0))
    set_maximal_velocities!(pbody, cbody, 
        parent_vertex=tra2.vertices[1], 
        child_vertex=tra2.vertices[2], 
        Δv=zeros(3), 
        Δω=zeros(3))
end
