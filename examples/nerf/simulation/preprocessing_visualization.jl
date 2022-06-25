using Dojo
using Plots
using OSFLoader
using DojoEnvironments

vis = Visualizer()
render(vis)
# open(vis)
render_static(vis)

################################################################################
# visuaize collider
################################################################################

function build_collider!(collider::SoftCollider{T,N}, vis::Visualizer;
        name::Symbol=:collider,
        visualize_particle::Bool=false,
        mesh_color=RGBA(1,1,1,0.3),
        particle_color=RGBA(1,1,1,1.0),
        gradient_color=RGBA(0.0,0.5,1,0.8),
        com_color=RGBA(0,0,0,1.0),
        particle_radius=0.012,
        ) where {T,N}
    # setobject!(vis[name][:mesh], collider.mesh, MeshPhongMaterial(color=mesh_color))
    setobject!(vis[name][:center_of_mass], HyperSphere(Point(collider.center_of_mass...), 0.05), MeshPhongMaterial(color=com_color))
    if visualize_particle
        for i = 1:10:length(collider.particles)
            particle = collider.particles[i]
            gradient = collider.density_gradients[i]
            setobject!(vis[name][:particles]["$i"], HyperSphere(Point(particle...), particle_radius), MeshPhongMaterial(color=particle_color))
            arrow_vis = ArrowVisualizer(vis[name][:gradients]["$i"])
            setobject!(arrow_vis, MeshPhongMaterial(color=gradient_color))
            settransform!(arrow_vis, Point(particle...), Vec(-1e-5gradient...), shaft_radius=0.005, max_head_radius=0.005)
        end
    end
    return nothing
end


nerf_model = :bunny
collider = SoftCollider(nerf=nerf_model)
mech = get_mechanism(:nerf, nerf=nerf_model)
p_offset = [-0.1,-0.5-0.2,0.38]
q_offset = Quaternion(normalize([1,+0.1,0.2,0.])...)
z = [p_offset + collider.center_of_mass; zeros(3); vector(q_offset); zeros(3)]

build_robot(mech, vis=vis, show_contact=false, color=RGBA(1, 165/255, 0.0, 0.15))
set_robot(vis, mech, z, show_contact=false)

build_collider!(collider, vis, visualize_particle=true)
settransform!(vis[:collider], compose(
    Translation(p_offset...),
    LinearMap(q_offset),
    ))
set_floor!(vis, x=2, y=2.0, z=0.001, color=RGBA(0.4, 0.4, 0.4, 0.4))
set_light!(vis, fill=0.8)
set_background!(vis)
# add centroid
centroid = Point(0.34,-0.1,-0.32)
setobject!(vis[:collider][:centroid], HyperSphere(centroid, 0.05),
    MeshPhongMaterial(color=RGBA(1,0,0,1)))

arrow_vis = ArrowVisualizer(vis[:collider][:force])
setobject!(arrow_vis, MeshPhongMaterial(color=RGBA(1,165/255,0,1)))
settransform!(arrow_vis, centroid, Vec(0,0,1.0),
    shaft_radius=0.025, max_head_radius=0.05)
# open(vis)
################################################################################
# animation
################################################################################
anim = MeshCat.Animation()
for i = 1:100
    atframe(anim, i) do
        settransform!(vis[:collider], compose(LinearMap(RotZ(0.01i*2π)), Translation(0,0,0.6)))
        # settransform!(vis, compose(LinearMap(RotZ(0.01i*2π)), Translation(0,0,0.6)))
        # settransform!(vis, LinearMap(RotZ(0.01i*2π)))
    end
end
setanimation!(vis, anim)



nerf_model = :bunny
collider = SoftCollider(nerf=nerf_model)
mech = get_mechanism(:nerf, nerf=nerf_model)
p_offset = [0,0,0.6]
q_offset = Quaternion(normalize([1,0,0,0.])...)
z = [p_offset + collider.center_of_mass; zeros(3); vector(q_offset); zeros(3)]

build_robot(mech, vis=vis, show_contact=false, color=RGBA(1, 165/255, 0.0, 0.25))
set_robot(vis, mech, z, show_contact=false)

build_collider!(collider, vis, visualize_particle=true)
settransform!(vis[:collider], compose(
    Translation(p_offset...),
    LinearMap(q_offset),
    ))
set_floor!(vis, x=2, y=2.0, z=0.001, color=RGBA(0.4, 0.4, 0.4, 0.4))
set_light!(vis, fill=0.8)
set_background!(vis)
# open(vis)


# render_static(vis)
# open(joinpath(@__DIR__, "preprocessing_visualization.html"), "w") do file
#     write(file, static_html(vis))
# end
