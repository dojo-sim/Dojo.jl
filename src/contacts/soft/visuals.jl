function build_collider!(collider::SoftCollider{T,N}, vis::Visualizer;
        name::Symbol=:collider, visualize_particle::Bool=false,
        mesh_color=RGBA(1,1,1,0.3),
        particle_color=RGBA(1,1,1,1.0),
        gradient_color=RGBA(0.0,0.5,1,0.8),
        com_color=RGBA(0,0,0,1.0),
        ) where {T,N}
    # setobject!(vis[name][:mesh], collider.mesh, MeshPhongMaterial(color=mesh_color))
    setobject!(vis[name][:center_of_mass], HyperSphere(Point(collider.center_of_mass...), 0.03), MeshPhongMaterial(color=com_color))
    if visualize_particle
        for i = 1:1000
            particle = collider.particles[i]
            gradient = collider.density_gradients[i]
            setobject!(vis[name][:particles]["$i"], HyperSphere(Point(particle...), 0.005), MeshPhongMaterial(color=particle_color))
            arrow_vis = ArrowVisualizer(vis[name][:gradients]["$i"])
            setobject!(arrow_vis, MeshPhongMaterial(color=gradient_color))
            settransform!(arrow_vis, Point(particle...), Vec(-1e-5gradient...), shaft_radius=0.0025, max_head_radius=0.0025)
        end
    end
    return nothing
end

# vis = Visualizer()
# open(vis)

set_light!(vis)
set_background!(vis)
nerf_object = get_nerf_object()
# collider = SoftCollider(nerf_object, N=1000)
# mech = get_bunny()
# collider = mech.contacts[1].model.collision.collider
build_collider!(collider, vis, visualize_particle=true)
build_robot(mech, vis=vis, color=RGBA(1,0.5,0,0.4))
set_floor!(vis, origin=[0,0,-1.0])
set_robot(vis, mech, [-com; zeros(3); 1; zeros(3); zeros(3);], show_contact=false)
com = mech.contacts[1].model.collision.collider_origin
