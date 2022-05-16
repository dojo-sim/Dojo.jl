
function get_nerf_sphere(;
    nerf::Symbol=:bunny,
    collider_options=ColliderOptions(),
    timestep=0.01,
    gravity=[0.0; 0.0; -9.81],
    friction_coefficient=0.8,
    radius=0.5,
    mass=10.0,
    T=Float64)

    collider = SoftCollider(nerf=nerf, opts=collider_options)

    assets_folder = joinpath(OSFLoader.osf_loader_dir(), "assets/mesh")
    inner_mesh_path = joinpath(assets_folder, String(nerf) * "_high.obj")
    outer_mesh_path = joinpath(assets_folder, String(nerf) * "_low.obj")
    soft_body = SoftBody(collider, inner_mesh_path, outer_mesh_path, name=nerf)

    origin = Origin{T}(name=:origin)
    sphere_origin = szeros(3)
    sphere = Sphere(radius, mass, name=:sphere)
    bodies = [soft_body, sphere]

    # joints
    joints = [
        JointConstraint(Floating(origin, bodies[1]), name=:joint_origin_nerf),
        JointConstraint(Floating(origin, bodies[2]), name=:joint_origin_sphere),
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



################################################################################
# Simulate nerf & sphere
################################################################################
mech = get_nerf_sphere(nerf=:bunny,
    collider_options=ColliderOptions(impact_spring=1e4),
    timestep=0.01,
    gravity=-9.81,
    friction_coefficient=0.1)

initialize!(mech, :nerf_sphere,
    nerf_position=[0,0,0],
    nerf_velocity=[0,0,0],
    sphere_position=[0,4.0,0.4],
    sphere_velocity=[0,-5.0,0],
    )
# Main.@profiler
@elapsed storage = simulate!(mech, 5.0,
    opts=SolverOptions(verbose=false, rtol=1e-4))
visualize(mech, storage, vis=vis)
mech.contacts[1]
mech.contacts[2]
mech.contacts[3]

similar(mech.system.matrix_entries, Float64)
mech.joints
mech.bodies
mech.contacts

mech.system.cyclic_children
mech.system.acyclic_children
mech.system.graph.fadjlist
mech.system.dfs_graph.fadjlist
