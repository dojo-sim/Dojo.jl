function get_bunny(;
    timestep=0.01,
    gravity=[0.0; 0.0; -9.81],
    friction_coefficient=0.8,
    radius=0.5,
    T=Float64)

    origin = Origin{T}(name=:origin)
    mass = 1.0
    bodies = [Sphere(radius, mass, name=:bunny, color=RGBA(1,1,1,0.1))]
    shape0 = bodies[1].shape
    shape1 = Mesh(joinpath(module_dir(), "environments/bunny/deps/bunny_inner_mesh.obj"), color=RGBA(0.2,0.2,0.2,1.0))
    shape2 = Mesh(joinpath(module_dir(), "environments/bunny/deps/bunny_outer_mesh.obj"), color=RGBA(0.9,0.9,0.9,0.3))
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
    x=zeros(3),
    q=one(Quaternion),
    v=zeros(3),
    ω=zeros(3)) where T

    r = mechanism.bodies[1].shape.r
    joint = get_joint(mechanism, :floating_joint)
    zero_velocity!(mechanism)
    set_minimal_coordinates!(mechanism, joint, [x + [0.0, 0.0, r] rotation_vector(q)])
    set_minimal_velocities!(mechanism, joint, [v; ω])
end
