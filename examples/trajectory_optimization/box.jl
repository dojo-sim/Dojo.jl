using Dojo
using IterativeLQR
using LinearAlgebra 

# ## system
gravity = -9.81
dt = 0.1
env = make("block", 
    mode=:max, 
    dt=dt,
    friction_coefficient=0.5,
    gravity=gravity)

@show env.mechanism.bodies[1].m
@show env.mechanism.bodies[1].J

# ## visualizer 
open(env.vis) 

# ## dimensions
n = env.nx
m = env.nu
d = 0

# ## states
z1 = [0.0; 0.0; 0.25; 0.0; 0.0; 0.0; 1.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0]
zT = [1.0; 0.0; 0.25; 0.0; 0.0; 0.0; 1.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0] # right goal
# zT = [0.0; 0.0; 0.25 + 1.0; 0.0; 0.0; 0.0; 1.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0] # up goal

# ## horizon
T = 11

# ## model
dyn = IterativeLQR.Dynamics(
    (y, x, u, w) -> f(y, env, x, u, w), 
    (dx, x, u, w) -> fx(dx, env, x, u, w),
    (du, x, u, w) -> fu(du, env, x, u, w),
    n, n, m, d)
model = [dyn for t = 1:T-1]

# ## rollout
ū = [[0.0; 0.0; 0.0] for t = 1:T-1]
w = [zeros(d) for t = 1:T-1]
x̄ = rollout(model, z1, ū, w)
visualize(env, x̄)

# ## objective
ot = (x, u, w) -> transpose(x - zT) * Diagonal(1.0 * ones(n)) * (x - zT) + transpose(u) * Diagonal(1.0e-2 * ones(m)) * u
oT = (x, u, w) -> transpose(x - zT) * Diagonal(1.0 * ones(n)) * (x - zT)

ct = Cost(ot, n, m, d)
cT = Cost(oT, n, 0, 0)
obj = [[ct for t = 1:T-1]..., cT]

# ## constraints
goal(x, u, w) = x - zT

cont = IterativeLQR.Constraint()
conT = IterativeLQR.Constraint(goal, n, 0)
cons = [[cont for t = 1:T-1]..., conT]

# ## problem 
prob = IterativeLQR.problem_data(model, obj, cons)
IterativeLQR.initialize_controls!(prob, ū)
IterativeLQR.initialize_states!(prob, x̄)

# ## solve
@time IterativeLQR.solve!(prob,
    linesearch=:armijo,
    α_min=1.0e-5,
    obj_tol=1.0e-3,
    grad_tol=1.0e-3,
    con_tol=0.005,
    max_iter=100,
    max_al_iter=10,
    ρ_init=1.0,
    ρ_scale=10.0,
    verbose=false)

# ## solution
z_sol, u_sol = IterativeLQR.get_trajectory(prob)
@show IterativeLQR.eval_obj(prob.m_data.obj.costs, prob.m_data.x, prob.m_data.u, prob.m_data.w)
@show prob.s_data.iter[1]
@show norm(goal(prob.m_data.x[T], zeros(0), zeros(0)), Inf)

# ## visualize
v, anim = visualize(env, [[z_sol[1] for t = 1:10]..., z_sol..., [z_sol[end] for t = 1:10]...])

set_camera!(env.vis, zoom=50.0, cam_pos=[100,0,0])
set_floor!(env.vis, x=0.0, y=4, z=0.02, color=RGBA(0.7,0.7,0.7,1))

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


