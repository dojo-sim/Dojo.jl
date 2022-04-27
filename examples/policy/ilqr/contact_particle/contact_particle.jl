using Pkg
Pkg.activate(joinpath(module_dir(), "examples"))

using Dojo
using IterativeLQR
const iLQR = IterativeLQR
using LinearAlgebra
using Plots
using Symbolics
using BenchmarkTools


# ## contact particle
vis = Visualizer()
open(vis)


include("model_contact_particle.jl")

timestep = 0.05
gravity = -9.81
rtol = 1e-7
btol = 1e-5
model = ContactParticle19(gravity=gravity, timestep=timestep)
add_symbolics!(model)
env = particle(timestep=timestep)


nx = model.nx
nu = model.nu

x_hist = [[0,0,1,2,1,0.0]]
for i = 1:100
    y = zeros(nx)
    dynamics(model, y, x_hist[end], [0,0,0.0], btol=1e-3, rtol=1e-5)
    push!(x_hist, y)
end
plot(hcat(x_hist...)'[:,1:3])

z_hist = [minimal_to_maximal(env.mechanism, x) for x in x_hist]
storage_hist = generate_storage(env.mechanism, z_hist)
visualize(env.mechanism, storage_hist, vis=vis)

x = [1,1,1,1,1,1.0]
u = [0,0,0.0]
dx = zeros(nx,nx)
du = zeros(nx,nu)
dynamics_jacobian_state(model, dx, x, u; btol=1e-3, rtol=1e-5)
dx
dynamics_jacobian_input(model, du, x, u; btol=1e-3, rtol=1e-5)
du


# ## initialization
x1 = [0.0; 0.0; 0.25; 0.0; 0.0; 0.0]
xT = [1.0; 2.0; 0.25; 0.0; 0.0; 0.0]
u_hover = [0.0, 0.0, 0.0]

# ## (1-layer) multi-layer perceptron policy
l_input = nx
l1 = 6
l2 = nu
nθ = l1 * l_input + l2 * l1

function policy(θ, x, goal)
    shift = 0
    # input
    input = x - goal

    # layer 1
    W1 = reshape(θ[shift .+ (1:(l1 * l_input))], l1, l_input)
    z1 = W1 * input
    o1 = tanh.(z1)
    shift += l1 * l_input

    # layer 2
    W2 = reshape(θ[shift .+ (1:(l2 * l1))], l2, l1)
    z2 = W2 * o1

    o2 = z2
    return o2
end

# ## horizon
T = 31

# ## model
h = timestep

function f1(y, x, u, w)
    u_ctrl = u[1:nu]
    x_di = x[1:nx]
    θ = u[nu .+ (1:nθ)]
    # dynamics(drone, h, x_di, u_ctrl, w);
    # dynamics(view(y, 1:nx), env, x_di, u_ctrl, w)
    dynamics(model, view(y, 1:nx), x_di, u_ctrl)
    y[nx .+ (1:nθ)] .= θ
    # [
    #     y;
    #     θ;
    # ]
end

function f1x(dx, x, u, w)
    u_ctrl = u[1:nu]
    x_di = x[1:nx]
    θ = u[nu .+ (1:nθ)]
    # dynamics(drone, h, x_di, u_ctrl, w);
    dx .= 0.0
    # dynamics_jacobian_state(view(dx, 1:nx, 1:nx), env, x_di, u_ctrl, w)
    dynamics_jacobian_state(model, view(dx, 1:nx, 1:nx), x_di, u_ctrl)
    # [
    #     dx;
    #     zeros(nθ,nx);
    # ]
end

function f1u(du, x, u, w)
    u_ctrl = u[1:nu]
    x_di = x[1:nx]
    θ = u[nu .+ (1:nθ)]
    # dynamics(drone, h, x_di, u_ctrl, w);
    du .= 0.0
    # dynamics_jacobian_input(view(du, 1:nx, 1:nu), env, x_di, u_ctrl, w)
    dynamics_jacobian_input(model, view(du, 1:nx, 1:nu), x_di, u_ctrl)
    du[nx .+ (1:nθ), nu .+ (1:nθ)] .= I(nθ)
    # [
    #     du zeros(nx,nθ);
    #     zeros(nθ,nu) I(nθ);
    # ]
end


function ft(y, x, u, w)
    u_ctrl = u[1:nu]
    x_di = x[1:nx]
    θ = x[nx .+ (1:nθ)]
    # dynamics(drone, h, x_di, u_ctrl, w);
    # dynamics(view(y, 1:nx), env, x_di, u_ctrl, w)
    dynamics(model, view(y, 1:nx), x_di, u_ctrl)
    y[nx .+ (1:nθ)] .= θ
    # [
    #     y
    #     θ;
    # ]
end

function ftx(dx, x, u, w)
    u_ctrl = u[1:nu]
    x_di = x[1:nx]
    θ = x[nx .+ (1:nθ)]
    # dynamics(drone, h, x_di, u_ctrl, w);
    dx .= 0.0
    # dynamics_jacobian_state(view(dx, 1:nx, 1:nx), env, x_di, u_ctrl, w)
    dynamics_jacobian_state(model, view(dx, 1:nx, 1:nx), x_di, u_ctrl)
    dx[nx .+ (1:nθ), nx .+ (1:nθ)] .= I(nθ)
    # [
    #     dx;
    #     zeros(nθ,nx);
    # ]
end

function ftu(du, x, u, w)
    u_ctrl = u[1:nu]
    x_di = x[1:nx]
    θ = x[nx .+ (1:nθ)]
    # dynamics(drone, h, x_di, u_ctrl, w);
    # dynamics_jacobian_input(view(du, 1:nx, 1:nu), env, x_di, u_ctrl, w)
    dynamics_jacobian_input(model, view(du, 1:nx, 1:nu), x_di, u_ctrl)
end



# user-provided dynamics and gradients
dyn1 = iLQR.Dynamics(f1, f1x, f1u, nx + nθ, nx, nu + nθ)
dynt = iLQR.Dynamics(ft, ftx, ftu, nx + nθ, nx + nθ, nu)

dyn = [dyn1, [dynt for t = 2:T-1]...]

# ## objective
function o1(x, u, w)
    J = 0.0
    q = [1.0e-1, 1.0e-1, 1.0e-1, 1.0e-1, 1.0e-1, 1.0e-1]
    r = [1.0e-2, 1.0e-2, 1.0e-0]
    ex = x - xT
    J += 0.5 * transpose(ex) * Diagonal(q) * ex
    J += 1.0e-1 * transpose(u[1:nu] - u_hover) * Diagonal(r) * (u[1:nu] - u_hover)
    J += 1.0e-1 * dot(u[nu .+ (1:nθ)], u[nu .+ (1:nθ)])
    return J
end

function ot(x, u, w)
    J = 0.0
    q = [1.0, 1.0, 1.0, 1.0e-1, 1.0e-1, 1.0e-1]
    r = [1.0e-2, 1.0e-2, 1.0e-0]
    ex = x[1:nx] - xT
    J += 0.5 * transpose(ex) * Diagonal(q) * ex
    J += 1.0e-1 * transpose(u[1:nu] - u_hover) * Diagonal(r) * (u[1:nu] - u_hover)
    J += 1.0e-1 * dot(x[nx .+ (1:nθ)], x[nx .+ (1:nθ)])
    return J
end

function ott(x, u, w)
    J = 0.0
    q = [10.0, 10.0, 10.0, 1.0e-1, 1.0e-1, 1.0e-1]
    r = [1.0e-2, 1.0e-2, 1.0e-0]
    ex = x[1:nx] - xT
    J += 0.5 * transpose(ex) * Diagonal(q) * ex
    J += 1.0e-1 * transpose(u[1:nu] - u_hover) * Diagonal(r) * (u[1:nu] - u_hover)
    J += 1.0e-1 * dot(x[nx .+ (1:nθ)], x[nx .+ (1:nθ)])
    return J
end

function oT(x, u, w)
    J = 0.0
    # q = [1000.0, 1000.0, 1000.0, 100.0, 100.0, 100.0]
    # ex = x[1:nx] - xT
    # J += 0.5 * transpose(ex) * Diagonal(q) * ex
    # J += 1.0e-1 * dot(x[nx .+ (1:nθ)], x[nx .+ (1:nθ)])
    return J
end

c1 = iLQR.Cost(o1, nx, nu + nθ)
ct = iLQR.Cost(ot, nx + nθ, nu)
ctt = iLQR.Cost(ott, nx + nθ, nu)
cT = iLQR.Cost(oT, nx + nθ, 0)
obj = [c1, [ct for t = 2:16]..., [ctt for t = 17:(T - 1)]..., cT]

# ## constraints
ul = -10.0 * ones(nu)
uu = +10.0 * ones(nu)
# bnd1 = iLQR.Bound(nx, nu + nθ, xl=x1, xu=x1, ul=[ul; -Inf * ones(nθ)], uu=[uu; Inf * ones(nθ)])
# bndt = iLQR.Bound(nx + nθ, nu, ul=ul, uu=uu)
# bndT = iLQR.Bound(nx + nθ, 0)# xl=[xT; -Inf * ones(nθ)], xu=[xT; Inf * ones(nθ)])
# bnds = [bnd1, [bndt for t = 2:T-1]..., bndT]

function con1(x, u, w)
    θ = u[nu .+ (1:nθ)]
    [
        ul - u[1:nu];
        u[1:nu] - uu;
        1.0 * (u[1:nu] - policy(θ, x[1:nx], xT));
    ]
end

function cont(x, u, w)
    θ = x[nx .+ (1:nθ)]
    [
        ul - u[1:nu];
        u[1:nu] - uu;
        1.0 * (u[1:nu] - policy(θ, x[1:nx], xT))
    ]
end

function goal(x, u, w)
    [
        x[1:nx] - xT[1:nx];
        # x[nx .+ (1:nθ)]
    ]
end
con_policy1 = iLQR.Constraint(con1, nx, nu + nθ, indices_inequality=collect(1:2nu))
con_policyt = iLQR.Constraint(cont, nx + nθ, nu, indices_inequality=collect(1:2nu))

cons = [con_policy1, [con_policyt for t = 2:T-1]..., iLQR.Constraint(goal, nx + nθ, 0)]

# ## problem
opts = iLQR.Options(line_search=:armijo,
    max_iterations=250,
    max_dual_updates=10,
    objective_tolerance=1e-3,
    lagrangian_gradient_tolerance=1e-3,
    constraint_tolerance=1e-3,
    # initial_constraint_penalty,
    # scaling_penalty,
    # max_penalty,
    # reset_cache,
    verbose=true)

# p = iLQR.problem_data(dyn, obj, cons)
p = iLQR.Solver(dyn, obj, cons, options=opts)

# ## initialize
θ0 = 1.0 * randn(nθ)
u_guess = [t == 1 ? [u_hover; θ0] : u_hover for t = 1:T-1]
x_guess = iLQR.rollout(dyn, x1, u_guess)

iLQR.initialize_controls!(p, u_guess)
iLQR.initialize_states!(p, x_guess)

# ## solve
@time iLQR.solve!(p)

# ## solution
x_sol, u_sol = iLQR.get_trajectory(p)
@show u_sol[1]
@show x_sol[1]
@show x_sol[T]
θ_sol = u_sol[1][nu .+ (1:nθ)]

# ## state
plot(hcat([x[1:nx] for x in x_sol]...)', label="", color=:orange, width=2.0)

# ## control
plot(hcat([u[1:nu] for u in u_sol]..., u_sol[end])', linetype = :steppost)

# ## plot xy
plot([x[1] for x in x_sol], [x[2] for x in x_sol], label="", color=:black, width=2.0)

# ## visualization
z_sol = [minimal_to_maximal(env.mechanism, x) for x in x_sol]
storage_sol = generate_storage(env.mechanism, z_sol)
visualize(env.mechanism, storage_sol, vis=vis)

# ## simulate policy
x_hist = [x1]
u_hist = [u_hover]

for t = 1:(3 * T)
    push!(u_hist, policy(θ_sol, x_hist[end], xT))
    y = zeros(nx)
    dynamics(model, y, x_hist[end], u_hist[end], btol=1e-3, rtol=1e-5)
    # dynamics(y, env, x_hist[end], u_hist[end], zeros(0))
    push!(x_hist, y)
end

z_hist = [minimal_to_maximal(env.mechanism, x) for x in x_hist]
storage_hist = generate_storage(env.mechanism, z_hist)
visualize(env.mechanism, storage_hist, vis=vis)

# convert_frames_to_video_and_gif("particle_slide_po_open_loop")
# convert_frames_to_video_and_gif("particle_slide_po_feedback")
