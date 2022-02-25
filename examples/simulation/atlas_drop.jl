using Dojo

# ## Open visualizer
vis = Visualizer()
open(vis)

# ## Mechanism
# mech = get_mechanism(:atlas, timestep=0.01, gravity=-9.81, friction_coefficient=0.5, damper=100.0, spring=1.0, contact=true)
mech = get_mechanism(:atlas, timestep=0.1, gravity=-9.81, friction_coefficient=0.5, damper=25.0, spring=1.0, contact=true, body_contact=true, model_type=:v2)

# ## Simulate
initialize_atlasstance!(mech, tran=[0,0,0.5], rot=[0.0,0.2,0.1])
storage = simulate!(mech, 2.5, record=true, opts=SolverOptions(rtol=1.0e-6, btol=1e-6))

A = full_matrix(mech.system)
Ne = sum([length(joint) for joint in mech.joints])
bi = 2
maximum([cond(A[Ne + (i - 1) * 6 .+ (1:6), Ne + (i - 1) * 6  .+ (1:6)]) for i = 1:length(mech.bodies)])

# ## Visualize
render(vis)
visualize(mech, storage, vis=vis, show_contact=true)

# ## Contact innerpenetration 
res = get_sdf(mech, storage)
minimum(minimum([min.(0.0, r) for r in res]))

# ## Animation
z = get_maximal_state(storage)
z = [[z[1] for t = 1:100]..., z..., [z[end] for t = 1:100]...]
T = length(z)
anim = MeshCat.Animation(convert(Int, floor(1.0 / 0.01)))
build_robot(mech, vis=vis)
for t = 1:T
    MeshCat.atframe(anim, t) do
        set_robot(vis, mech, z[t])
    end
end
MeshCat.setanimation!(vis, anim)
set_camera!(vis, cam_pos=[3,-0,0], zoom=0.84)

set_floor!(vis, x=0.2, y=10, z=0.02, color=RGBA(0,0,0,1))

# ## HMTL scene
render_static(vis)
open(joinpath(@__DIR__, "atlas_drop.html"), "w") do file
    write(file, static_html(vis))
end
