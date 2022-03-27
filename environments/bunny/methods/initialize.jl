function get_bunny(;
    collider="soft_collider.jld2",
    timestep=0.01,
    gravity=[0.0; 0.0; -9.81],
    friction_coefficient=0.8,
    radius=0.5,
    T=Float64)

    deps_folder = joinpath(module_dir(), "environments/bunny/deps")
    collider = jldopen(joinpath(deps_folder, collider))["soft"]
    inner_mesh_path = joinpath(deps_folder, "bunny_inner_mesh.obj")
    outer_mesh_path = joinpath(deps_folder, "bunny_outer_mesh.obj")
    soft_body = SoftBody(collider, inner_mesh_path, outer_mesh_path, name=:bunny)

    origin = Origin{T}(name=:origin)
    bodies = [soft_body]

    joints = [JointConstraint(Floating(origin, bodies[1]), name=:floating_joint)]
    mechanism = Mechanism(origin, bodies, joints,
        timestep=timestep,
        gravity=gravity)

    contact = [0.0, 0.0, 0.0]
    normal = [0.0, 0.0, 1.0]
    contacts = [soft_contact_constraint(get_body(mechanism, :bunny), normal, collider,
        friction_coefficient=friction_coefficient,
        contact_origin=contact,
        contact_radius=radius,
        name=:contact_1)]
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
