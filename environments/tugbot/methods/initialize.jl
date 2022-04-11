function get_tugbot(;
    timestep=0.01,
    gravity=[0.0; 0.0; -9.81],
    friction_coefficient=0.2,
    radius=0.2,
    object_dimension=[0.4, 0.4, 0.1],
    contact=true,
    contact_type=:nonlinear,
    T=Float64)

    origin = Origin{T}(name=:origin)
    drone_mass = 1.0
    object_mass = 2.0
    bodies = [
        Sphere(radius, drone_mass, name=:drone),
        Box(object_dimension..., object_mass, name=:object)]
    joints = [
        JointConstraint(Floating(origin, bodies[1]), name=:drone_joint)
        # JointConstraint(Floating(origin, bodies[2]), name=:object_joint)
        ]
    mechanism = Mechanism(origin, bodies, joints,
        timestep=timestep,
        gravity=gravity)

    if contact
        contacts = [[0.2, -0.2, -0.05], [-0.2, 0.2, -0.05], [0.2, 0.2, -0.05], [-0.2, -0.2, -0.05]]
        normal = [0.0, 0.0, 1.0]
        contacts = [contact_constraint(get_body(mechanism, :object), normal,
            friction_coefficient=friction_coefficient,
            contact_origin=contacts[i],
            contact_radius=0.0,
            contact_type=contact_type,
            name=Symbol(:contact,i)) for i=1:length(contacts)]

        collision = StringCollision{Float64,0,3,0}(
                0*[0,0,-0.2],
                1*[0.2, 0.2, 0.05],
                2.0)
        parameterization = szeros(T, 0, 2)
        body_body_contact = ImpactContact{Float64,2}(parameterization, collision)

        contacts = [contacts; ContactConstraint((body_body_contact,
            get_body(mechanism, :drone).id,
            get_body(mechanism, :object).id,), name=:body_body)]
        mechanism = Mechanism(origin, bodies, joints, contacts,
            gravity=gravity,
            timestep=timestep)
    end
    return mechanism
end

function initialize_tugbot!(mechanism::Mechanism{T};
    drone_position=[0,1,1.0],
    drone_orientation=one(Quaternion),
    drone_velocity=zeros(3),
    drone_angular_velocity=zeros(3),
    object_position=zeros(3),
    object_orientation=one(Quaternion),
    object_velocity=zeros(3),
    object_angular_velocity=zeros(3),
    ) where T

    r_drone = get_body(mechanism, :drone).shape.r
    r_object = get_body(mechanism, :object).shape.xyz[3]/2
    joint = get_joint(mechanism, :drone_joint)
    zero_velocity!(mechanism)
    z_drone = [drone_position + [0,0,r_drone]; drone_velocity; vector(drone_orientation); drone_angular_velocity]
    z_object = [object_position + [0,0,r_object]; object_velocity; vector(object_orientation); object_angular_velocity]
    set_maximal_state!(mechanism, [z_drone; z_object])
end
