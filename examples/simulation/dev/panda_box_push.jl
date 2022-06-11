using MeshCat
using GeometryBasics

using Dojo
using DojoEnvironments

vis = Visualizer()
render(vis)

include(joinpath(@__DIR__, "..", "..", "..", "DojoEnvironments/src/utilities.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "DojoEnvironments/src/panda/methods/initialize.jl"))

mech = get_panda(
    timestep=0.01,
    gravity=0.0 * 9.81,
    spring=1.0,
    damper=100.0,
    # contact=true,
    model_type=:no_end_effector,
    limits=false,
    object_type=:box,
)

nu = input_dimension(mech)
function ctrl!(m, k; kp=10.0, kd=5.0)
    pref = [0,-0.8,0,1.6,0,-2.4]#,0.5 * π]
    y = get_minimal_state(m)
    p = y[12 .+ (1:2:12)]
    v = y[12 .+ (2:2:12)]
    u = [
            5.0 * kp * [1,1.0e-1,1.0e-1,1.0e-1,1.0e-1,1.0e-1] .* (pref - p) - kd * v;   
            [0.0; 0.0; -1.0 * get_body(m, :box).mass * 9.81];
            zeros(3);         
    ]
    set_input!(m, u*m.timestep)
    return nothing
end

joint_angles = [1.3,-0.85,0.3,1.3,0.3,-2.0]#,0.5 * π]
joint_velocities = [0,0,0,0,0,0]#,0]

y = zeros(length(get_minimal_state(mech)))
y[1:3] = [0.5; 0.5; 0.125]
for i = 1:6
    y[12 + (i - 1) * 2 .+ (1:2)] = [joint_angles[i]; joint_velocities[i]]
end
set_minimal_state!(mech, y)
# set_minimal_state!(mech, copy(y_start[12 .+ (1:18)]))
# set_minimal_state!(mech, copy(y_start))

box_body = get_body(mech, :box)
box_body.state.q2 =  RotZ(0.25 * π)

storage = simulate!(mech, 5.0, 
    ctrl!, 
    verbose=true,
    opts=SolverOptions(verbose=true, rtol=1.0e-6, btol=1.0e-6));

vis, anim = visualize(mech, storage, vis=vis, show_contact=false)

setobject!(vis["ee"], GeometryBasics.Sphere(Point3f0(0),
        convert(Float32, mech.contacts[1].model.collision.radius_sphere)),
        MeshPhongMaterial(color = RGBA(0.0, 1.0, 0.0, 1.0)))


for t in 1:length(storage)
    MeshCat.atframe(anim, t) do 
        MeshCat.settransform!(vis["ee"], 
        MeshCat.Translation(storage.x[5][t] + vector_rotate(mech.contacts[1].model.collision.origin_sphere, storage.q[5][t])))
    end
end
MeshCat.setanimation!(vis, anim)

