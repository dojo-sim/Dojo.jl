function get_bunny(collider::Collider;
    timestep=0.01,
    gravity=[0.0; 0.0; -9.81],
    friction_coefficient=0.8,
    radius=0.5,
    T=Float64)

    origin = Origin{T}(name=:origin)
    mass = 1.0
    bodies = [Sphere(radius, mass, name=:bunny, color=RGBA(1,1,1,0.1))]
    shape0 = bodies[1].shape
    shape1 = Mesh(joinpath(module_dir(), "environments/bunny/deps/bunny_inner_mesh.obj"),
        position_offset = -collider.center_of_mass,
        color=RGBA(0.2,0.2,0.2,1.0))
    shape2 = Mesh(joinpath(module_dir(), "environments/bunny/deps/bunny_outer_mesh.obj"),
        position_offset = -collider.center_of_mass,
        color=RGBA(0.9,0.9,0.9,0.3))
    shape_vec = [shape1, shape2, shape0]
    shapes = Shapes(shape_vec)
    bodies[1].shape = shapes

    joints = [JointConstraint(Floating(origin, bodies[1]), name=:floating_joint)]
    mechanism = Mechanism(origin, bodies, joints,
        timestep=timestep,
        gravity=gravity)

    contact = [0.0, 0.0, 0.0]
    normal = [0.0, 0.0, 1.0]
    contacts = [contact_constraint(get_body(mechanism, :bunny), normal,
        friction_coefficient=friction_coefficient,
        contact_origin=contact,
        contact_radius=radius,
        contact_type=:soft)]
    set_minimal_coordinates!(mechanism, get_joint(mechanism, :floating_joint), [0.0; 0.0; radius; zeros(3)])
    mechanism = Mechanism(origin, bodies, joints, contacts,
        gravity=gravity,
        timestep=timestep)
    return mechanism
end

function initialize_bunny!(mechanism::Mechanism{T};
    position=zeros(3),
    orientation=one(Quaternion),
    velocity=zeros(3),
    angular_velocity=zeros(3)) where T

    r = 0.50
    joint = get_joint(mechanism, :floating_joint)
    zero_velocity!(mechanism)
    set_minimal_coordinates!(mechanism, joint, [position + [0.0, 0.0, r] rotation_vector(orientation)])
    set_minimal_velocities!(mechanism, joint, [velocity; angular_velocity])
end
