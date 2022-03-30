function get_bunny(;
    collider="soft_collider.jld2",
    timestep=0.01,
    gravity=[0.0; 0.0; -9.81],
    friction_coefficient=0.8,
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
        name=:contact_1)]
    set_minimal_coordinates!(mechanism, get_joint(mechanism, :floating_joint), [0.0; 0.0; 0.5; zeros(3)])
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

function get_bunny_sphere(;
    collider="soft_collider.jld2",
    timestep=0.01,
    gravity=[0.0; 0.0; -9.81],
    friction_coefficient=0.8,
    radius=0.5,
    mass=1.0,
    T=Float64)

    deps_folder = joinpath(module_dir(), "environments/bunny/deps")
    collider = jldopen(joinpath(deps_folder, collider))["soft"]
    inner_mesh_path = joinpath(deps_folder, "bunny_inner_mesh.obj")
    outer_mesh_path = joinpath(deps_folder, "bunny_outer_mesh.obj")

    origin = Origin{T}(name=:origin)
    bunny = SoftBody(collider, inner_mesh_path, outer_mesh_path, name=:bunny)
    sphere = Sphere(radius, mass, name=:sphere)
    bodies = [bunny, sphere]

    # joints
    joints = [
        JointConstraint(Floating(origin, bodies[1]), name=:joint_origin_bunny),
        # JointConstraint(Floating(origin, bodies[2]), name=:joint_origin_sphere),
        # JointConstraint(Floating(bodies[1], bodies[2]), name=:joint_bunny_sphere),
        ]
    mechanism = Mechanism(origin, bodies, joints,
        timestep=timestep,
        gravity=gravity)

    # Contacts
    normal = [0,0,1.0]
    contact_bunny_halfspace = soft_contact_constraint(get_body(mechanism, :bunny), normal, collider,
        friction_coefficient=friction_coefficient,
        contact_origin=[0,0,0.0],
        name=:contact_bunny_halfspace)
    contact_sphere_halfspace = contact_constraint(get_body(mechanism, :sphere), normal,
        friction_coefficient=friction_coefficient,
        contact_origin=[0,0,0.0],
        contact_radius=radius,
        contact_type=:nonlinear,
        name=:contact_sphere_halfspace)
    model = SoftContact(get_body(mechanism, :bunny), normal, collider)
    contact_bunny_sphere = SoftContactConstraint((
        model,
        get_body(mechanism, :bunny).id,
        get_body(mechanism, :sphere).id); name=:contact_bunny_sphere)
    contacts = [contact_bunny_halfspace; contact_sphere_halfspace; contact_bunny_sphere]

    mechanism = Mechanism(origin, bodies, joints, contacts,
        gravity=gravity,
        timestep=timestep)
    return mechanism
end

function initialize_bunny_sphere!(mechanism::Mechanism{T};
    bunny_position=[0,0,0],
    sphere_position=[2,0,0],
    bunny_orientation=one(Quaternion),
    sphere_orientation=one(Quaternion),
    bunny_velocity=zeros(3),
    sphere_velocity=zeros(3),
    bunny_angular_velocity=zeros(3),
    sphere_angular_velocity=zeros(3)) where T

    r = 0.50
    z_bunny = [bunny_position + [0,0,r]; bunny_velocity; vector(bunny_orientation); bunny_angular_velocity]
    z_sphere = [sphere_position + [0,0,r]; sphere_velocity; vector(sphere_orientation); sphere_angular_velocity]
    z = [z_bunny; z_sphere]
    set_maximal_state!(mechanism, z)
end


#
# function get_sphere(;
#     timestep=0.01,
#     gravity=[0.0; 0.0; -9.81],
#     friction_coefficient=0.8,
#     radius=0.5,
#     contact=true,
#     contact_type=:nonlinear,
#     T=Float64)
#
#     origin = Origin{T}(name=:origin)
#     mass = 1.0
#     bodies = [Sphere(radius, mass, name=:sphere)]
#     joints = [JointConstraint(Floating(origin, bodies[1]), name=:floating_joint)]
#     mechanism = Mechanism(origin, bodies, joints,
#         timestep=timestep,
#         gravity=gravity)
#
#     if contact
#         contact = [0.0, 0.0, 0.0]
#         normal = [0.0, 0.0, 1.0]
#         contacts = [contact_constraint(get_body(mechanism, :sphere), normal,
#             friction_coefficient=friction_coefficient,
#             contact_origin=contact,
#             contact_radius=radius,
#             contact_type=contact_type)]
#         set_minimal_coordinates!(mechanism, get_joint(mechanism, :floating_joint), [0.0; 0.0; radius; zeros(3)])
#         mechanism = Mechanism(origin, bodies, joints, contacts,
#             gravity=gravity,
#             timestep=timestep)
#     end
#     return mechanism
# end
