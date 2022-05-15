function get_nerf(;
    nerf::Symbol,
    collider_options=ColliderOptions(),
    timestep=0.01,
    gravity=[0.0; 0.0; -9.81],
    friction_coefficient=0.8,
    T=Float64)

    collider = SoftCollider(nerf=nerf, opts=collider_options)

    assets_folder = joinpath(OSFLoader.osf_loader_dir(), "assets/mesh")
    inner_mesh_path = joinpath(assets_folder, String(nerf) * "_low.obj")
    outer_mesh_path = joinpath(assets_folder, String(nerf) * "_high.obj")
    soft_body = SoftBody(collider, inner_mesh_path, outer_mesh_path, name=nerf)

    origin = Origin{T}(name=:origin)
    bodies = [soft_body]

    joints = [JointConstraint(Floating(origin, bodies[1]), name=:floating_joint)]
    mechanism = Mechanism(origin, bodies, joints,
        timestep=timestep,
        gravity=gravity)

    contact = [0.0, 0.0, 0.0]
    normal = [0.0, 0.0, 1.0]
    contacts = [soft_contact_constraint(get_body(mechanism, nerf), normal, collider,
        friction_coefficient=friction_coefficient,
        collider_origin=-collider.center_of_mass,
        name=:contact_1)]
    set_minimal_coordinates!(mechanism, get_joint(mechanism, :floating_joint), [0.0; 0.0; 0.5; zeros(3)])
    mechanism = Mechanism(origin, bodies, joints, contacts,
        gravity=gravity,
        timestep=timestep)
    return mechanism
end

function initialize_nerf!(mechanism::Mechanism{T};
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

function get_nerf_sphere(;
    nerf::Symbol,
    collider_options=ColliderOptions(),
    timestep=0.01,
    gravity=[0.0; 0.0; -9.81],
    friction_coefficient=0.8,
    radius=0.5,
    mass=10.0,
    T=Float64)

    collider = SoftCollider(nerf=nerf, opts=collider_options)

    assets_folder = joinpath(OSFLoader.osf_loader_dir(), "assets/mesh")
    inner_mesh_path = joinpath(assets_folder, String(nerf) * "_low.obj")
    outer_mesh_path = joinpath(assets_folder, String(nerf) * "_high.obj")
    soft_body = SoftBody(collider, inner_mesh_path, outer_mesh_path, name=nerf)

    origin = Origin{T}(name=:origin)
    sphere_origin = szeros(3)
    sphere = Sphere(radius, mass, name=:sphere)
    bodies = [soft_body, sphere]

    # joints
    joints = [
        JointConstraint(Floating(origin, bodies[1]), name=:joint_origin_nerf),
        # JointConstraint(Floating(origin, bodies[2]), name=:joint_origin_sphere),
        # JointConstraint(Floating(bodies[1], bodies[2]), name=:joint_nerf_sphere),
        ]
    mechanism = Mechanism(origin, bodies, joints,
        timestep=timestep,
        gravity=gravity)

    # Contacts
    normal = [0,0,1.0]
    contact_nerf_halfspace = soft_contact_constraint(get_body(mechanism, nerf), normal, collider,
        friction_coefficient=friction_coefficient,
        collider_origin=[0,0,0.0],
        name=:contact_nerf_halfspace)
    contact_sphere_halfspace = contact_constraint(get_body(mechanism, :sphere), normal,
        friction_coefficient=friction_coefficient,
        contact_origin=sphere_origin,
        contact_radius=radius,
        contact_type=:nonlinear,
        name=:contact_sphere_halfspace)
    model = SoftContact(get_body(mechanism, nerf), normal, collider, child_origin=sphere_origin,
        radius=radius, collision_type=:soft_sphere)
    contact_nerf_sphere = SoftContactConstraint((
        model,
        get_body(mechanism, nerf).id,
        get_body(mechanism, :sphere).id); name=:contact_nerf_sphere)
    contacts = [contact_nerf_halfspace; contact_sphere_halfspace; contact_nerf_sphere]
    # contacts = [contact_nerf_sphere]

    mechanism = Mechanism(origin, bodies, joints, contacts,
        gravity=gravity,
        timestep=timestep)
    return mechanism
end

function initialize_nerf_sphere!(mechanism::Mechanism{T};
    nerf_position=[0,0,0],
    sphere_position=[2,0,0],
    nerf_orientation=one(Quaternion),
    sphere_orientation=one(Quaternion),
    nerf_velocity=zeros(3),
    sphere_velocity=zeros(3),
    nerf_angular_velocity=zeros(3),
    sphere_angular_velocity=zeros(3)) where T

    r = 0.50
    z_nerf = [nerf_position + [0,0,r]; nerf_velocity; vector(nerf_orientation); nerf_angular_velocity]
    z_sphere = [sphere_position + [0,0,r]; sphere_velocity; vector(sphere_orientation); sphere_angular_velocity]
    z = [z_nerf; z_sphere]
    set_maximal_state!(mechanism, z)
end
#
# function get_bunny_triumvirate(;
#     collider="bunny_collider.jld2",
#     timestep=0.01,
#     gravity=[0.0; 0.0; -9.81],
#     friction_coefficient=0.8,
#     radius=0.5,
#     mass=10.0,
#     T=Float64)
#
#     assets_folder = joinpath(module_dir(), "environments/bunny/deps")
#     collider1 = jldopen(joinpath(assets_folder, collider))["collider"]
#     collider2 = deepcopy(collider1)
#     inner_mesh_path = joinpath(assets_folder, "bunny_inner_mesh.obj")
#     outer_mesh_path = joinpath(assets_folder, "bunny_outer_mesh.obj")
#
#     origin = Origin{T}(name=:origin)
#     bunny1 = SoftBody(collider1, inner_mesh_path, outer_mesh_path, name=:bunny1)
#     bunny2 = SoftBody(collider2, inner_mesh_path, outer_mesh_path, name=:bunny2)
#     sphere_origin = szeros(3)
#     sphere = Sphere(radius, mass, name=:sphere)
#     bodies = [bunny1, bunny2, sphere]
#
#     # joints
#     joints = [JointConstraint(Floating(origin, bodies[1]), name=:joint_origin_bunny)]
#
#     mechanism = Mechanism(origin, bodies, joints,
#         timestep=timestep,
#         gravity=gravity)
#
#     # Contacts
#     normal = [0,0,1.0]
#     contact_bunny1_halfspace = soft_contact_constraint(get_body(mechanism, :bunny1), normal, collider1,
#         friction_coefficient=friction_coefficient,
#         collider_origin=[0,0,0.0],
#         name=:contact_bunny1_halfspace)
#     contact_bunny2_halfspace = soft_contact_constraint(get_body(mechanism, :bunny2), normal, collider2,
#         friction_coefficient=friction_coefficient,
#         collider_origin=[0,0,0.0],
#         name=:contact_bunny2_halfspace)
#     contact_sphere_halfspace = contact_constraint(get_body(mechanism, :sphere), normal,
#         friction_coefficient=friction_coefficient,
#         contact_origin=sphere_origin,
#         contact_radius=radius,
#         contact_type=:nonlinear,
#         name=:contact_sphere_halfspace)
#     model1 = SoftContact(get_body(mechanism, :bunny1), normal, collider1, child_origin=sphere_origin,
#         radius=radius, collision_type=:soft_sphere)
#     model2 = SoftContact(get_body(mechanism, :bunny2), normal, collider2, child_origin=sphere_origin,
#         radius=radius, collision_type=:soft_sphere)
#     model3 = SoftContact(get_body(mechanism, :bunny1), normal, collider1, child_collider=collider2,
#         radius=radius, collision_type=:soft_soft)
#     contact_bunny1_sphere = SoftContactConstraint((
#         model1,
#         get_body(mechanism, :bunny1).id,
#         get_body(mechanism, :sphere).id); name=:contact_bunny1_sphere)
#     contact_bunny2_sphere = SoftContactConstraint((
#         model2,
#         get_body(mechanism, :bunny2).id,
#         get_body(mechanism, :sphere).id); name=:contact_bunny2_sphere)
#     contact_bunny1_bunny2 = SoftContactConstraint((
#         model3,
#         get_body(mechanism, :bunny1).id,
#         get_body(mechanism, :bunny2).id); name=:contact_bunny1_bunny2)
#     contacts = [contact_bunny1_halfspace;
#         contact_bunny2_halfspace;
#         contact_sphere_halfspace;
#         contact_bunny1_sphere;
#         contact_bunny2_sphere;
#         contact_bunny1_bunny2]
#
#     mechanism = Mechanism(origin, bodies, joints, contacts,
#         gravity=gravity,
#         timestep=timestep)
#     return mechanism
# end
#
# function initialize_bunny_triumvirate!(mechanism::Mechanism{T};
#     positions=[[0,0,0] for i=1:3],
#     orientations=[one(Quaternion) for i=1:3],
#     velocities=[zeros(3) for i=1:3],
#     angular_velocities=[zeros(3) for i=1:3]) where T
#
#     r = 0.50
#     z_bunny1 = [positions[1] + [0,0,r]; velocities[1]; vector(orientations[1]); angular_velocities[1]]
#     z_bunny2 = [positions[2] + [0,0,r]; velocities[2]; vector(orientations[2]); angular_velocities[2]]
#     z_sphere = [positions[3] + [0,0,r]; velocities[3]; vector(orientations[3]); angular_velocities[3]]
#     z = [z_bunny1; z_bunny2; z_sphere]
#     set_maximal_state!(mechanism, z)
# end
