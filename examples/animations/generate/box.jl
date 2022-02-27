# x goal
shift = [0.25; 0.0; 0.0]
u = 0.65 * [1.0; 0.0; 0.0] ./ norm(ones(3))

# z goal
shift = [0.0; 0.0; 0.25]
u = 0.65 * [0.0; 0.0; 1.0] ./ norm(ones(3))

anim = MeshCat.Animation(convert(Int, floor(1.0 / dt)))
v = env.vis
z = z_sol[1]
u_max = maximum([norm(u) for u in u_sol])
force_vis = ArrowVisualizer(v[:force])
setobject!(force_vis, MeshPhongMaterial(color=orange))
settransform!(force_vis,
    Point(z[1] - u[1] - shift[1], z[2] - u[2] - shift[2], z[3] - u[3] - shift[3]),
    Vec(u[1], u[2], u[3]),
    shaft_radius=0.05,
    max_head_radius=0.1)

z_vis = [[z_sol[1] for t = 1:15]..., z_sol..., [z_sol[end] for t = 1:15]...]
u_vis = [[u_sol[1] for t = 1:15]..., u_sol..., [u_sol[end] for t = 1:15]...]
for t = 1:length(z_vis)
    z = z_vis[t]
    u = (t == length(z_vis) ? 0.5 * u_vis[end] ./ u_max : 0.5 * u_vis[t] ./ u_max)
    MeshCat.atframe(anim, t) do
        
        settransform!(v[:robot], MeshCat.compose(MeshCat.Translation(z[2], z[1], z[3]), MeshCat.LinearMap(UnitQuaternion(z[6 .+ (1:4)]...))))
        settransform!(force_vis,
            Point(z[2] - u[2] - shift[2], z[1] - u[1] - shift[1], z[3] - u[3] - shift[3]),
            Vec(u[2], u[1], u[3]),
            shaft_radius=0.05,
            max_head_radius=0.1)
        # set_floor!(env.vis, x=0.0, y=2, z=0.02, color=RGBA(0.7,0.7,0.7,1))
    end 
end
MeshCat.setanimation!(v, anim)


