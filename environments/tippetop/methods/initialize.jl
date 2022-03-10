function get_tippetop(;
    timestep=0.01,
    gravity=[0.0; 0.0; -9.81],
    friction_coefficient=0.4,
    contact=true,
    contact_type=:nonlinear,
    T=Float64)

    origin = Origin{T}(name=:origin)
    radius = 0.5
    mass = 1.0
    α = 0.2
    bodies = [
        Sphere(radius, mass, name=:sphere1, color=gray_light),
        Sphere(radius*α, mass*α^3, name=:sphere2, color=gray_light)]

    joints = [
        JointConstraint(Floating(origin, bodies[1]),
                name=:floating_joint),
        JointConstraint(Fixed(bodies[1], bodies[2],
                parent_vertex=[0,0,radius],
                child_vertex=zeros(3)),
                name = :fixed_joint)]
    mechanism = Mechanism(origin, bodies, joints,
        timestep=timestep,
        gravity=gravity)

    # modify inertias
    mechanism.bodies[1].inertia = Diagonal([1.9, 2.1, 2.0])

    if contact
        contact_point = [0.0, 0.0, 0.0]
        normal = [0.0, 0.0, 1.0]
        contacts = [
            contact_constraint(get_body(mechanism, :sphere1), normal,
                friction_coefficient=friction_coefficient,
                contact_point=contact_point, contact_radius=radius,
                contact_type=contact_type),
            contact_constraint(get_body(mechanism, :sphere2), normal,
                friction_coefficient=friction_coefficient,
                contact_point=contact_point,
                contact_radius=radius * α,
                contact_type=contact_type)
            ]
        set_minimal_coordinates!(mechanism, get_joint(mechanism, :floating_joint), [0.0; 0.0; radius; zeros(3)])
        mechanism = Mechanism(origin, bodies, joints, contacts,
            gravity=gravity,
            timestep=timestep)
    end

    return mechanism
end

function initialize_tippetop!(mechanism::Mechanism{T};
    body_position=zeros(3),
    body_orientation=zeros(3),
    body_linear_velocity=zeros(3),
    body_angular_velocity=zeros(3)) where T

    floating_joint = get_joint(mechanism, :floating_joint)
    fixed_joint = get_joint(mechanism, :fixed_joint)
    radius = fixed_joint.translational.vertices[1][3]

    zero_velocity!(mechanism)
    set_minimal_coordinates_velocities!(mechanism, floating_joint,
        xmin=[body_position; body_orientation; body_linear_velocity; body_angular_velocity])
    set_minimal_coordinates_velocities!(mechanism, fixed_joint)
    return nothing
end
