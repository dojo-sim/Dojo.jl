function build_mesh!(mesh, vis::Visualizer; name::Symbol=:robot, color=RGBA(1,1,1,1.0))
    setobject!(vis[name], mesh, MeshPhongMaterial(color=color))
    return nothing
end

function build_sphere!(sphere, vis::Visualizer; name::Symbol=:robot, color=RGBA(1,1,1,1.0))
    setobject!(vis[name], HyperSphere(Point(sphere.origin...), sphere.radius),
        MeshPhongMaterial(color=color))
    return nothing
end

function set_robot!(x, q, vis::Visualizer; name::Symbol=:robot)
    settransform!(vis[name], MeshCat.compose(
        MeshCat.Translation(x),
        MeshCat.LinearMap(q),)
        )
    return nothing
end

function build_collider!(collider::SoftCollider110{N,T}, vis::Visualizer;
        name::Symbol=:collider, visualize_particle::Bool=false,
        mesh_color=RGBA(1,1,1,0.3),
        particle_color=RGBA(0.8,0.8,0.8,1.0),
        gradient_color=RGBA(0.9,0.2,0.2,1.0),
        com_color=RGBA(0,0,0,1.0),
        ) where {T,N}
    setobject!(vis[name][:mesh], collider.mesh, MeshPhongMaterial(color=mesh_color))
    setobject!(vis[name][:center_of_mass], HyperSphere(Point(collider.center_of_mass...), 0.03), MeshPhongMaterial(color=com_color))
    if visualize_particle
        for i = 1:N
            particle = collider.particles[i]
            gradient = collider.gradients[i]
            setobject!(vis[name][:particles]["$i"], HyperSphere(Point(particle...), 0.005), MeshPhongMaterial(color=particle_color))
            arrow_vis = ArrowVisualizer(vis[name][:gradients]["$i"])
            setobject!(arrow_vis, MeshPhongMaterial(color=gradient_color))
            settransform!(arrow_vis, Point(particle...), Vec(-2e-5gradient...), shaft_radius=0.0025)
        end
    end
    return nothing
end
