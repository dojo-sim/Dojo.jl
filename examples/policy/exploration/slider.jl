using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using Dojo
using Plots
using ForwardDiff

vis = Visualizer()
open(vis)

include("utils.jl")
include("solver.jl")
include("methods.jl")

################################################################################
# build environment
################################################################################
H0 = 20
N0 = 2
timestep = 0.05
gravity = [0, 0, 0.0]
spring = 20.0
damper = 1.0
env = slider(timestep=timestep, gravity=gravity, spring=spring, damper=damper, vis=vis)
nx = minimal_dimension(env.mechanism)
nu = input_dimension(env.mechanism)
nθ = (nx+1) * nu
# nva = nθ + nx*(N0-1) + nu*(H0-1)
nva = nθ + nu*(H0-1)


function dummy_dynamics(y, dx, du, x, u)
	A = [1 timestep
		 0 1]
	B = 1/timestep * [0.5*timestep^2 timestep]'
	y .= A * x + B * u
	dx .= A
	du .= B
	return nothing
end

################################################################################
# reference trajectory
################################################################################
initialize!(env.mechanism, :slider, position=0.50)
ref_traj = simulate!(env.mechanism, H0*timestep)
visualize(env.mechanism, ref_traj, vis=vis)
Zref = get_maximal_state(ref_traj)
Xref = [maximal_to_minimal(env.mechanism, Zref[i]) for i=1:H0]
x0 = Xref[1]
xF = Xref[H0]
plot([x[1] for x in Xref])


################################################################################
# multiple shooting
################################################################################
splits0 = split_trajectory(H0, N0)
# vars0 = [0;0;0; vcat([Xref[s[end]] for s in splits0[1:N0-1]]...); zeros(nu*(H0-1))]
vars0 = [0;0;0; zeros(nu*(H0-1))]
X0 = policy_rollout_only(env, x0, splits0, vars0)
plot_rollout(X0, Xref, splits0)


Q0 = 1*1e-1*Diagonal(ones(nx))
Qθ0 = 1*1e-2*Diagonal(ones(nx+1))
Qstart0 = 0*1e-5*Diagonal(ones(nx))
R0 = 0*1e-5*Diagonal(ones(nu))
λx0 = [zeros(nx) for i=1:N0-1]
λu0 = [zeros(nu) for i=1:H0-1]
ρx0 = 1*1e0
ρu0 = 1e-0
u_scale0 = 1e-0# for N=9 ->  1e-3
obj0 = AugmentedObjective120(Q0, Qθ0, Qstart0, R0, λx0, λu0, ρx0, ρu0, u_scale0)

vars1 = augmented_lagrangian_solver(env, obj0, Xref, splits0, vars0)
Xpolicy = pure_policy_rollout(env, x0, N0, H0, vars1, Hsim=H0, mode=:policy)
Xcontrol = pure_policy_rollout(env, x0, N0, H0, vars1, Hsim=H0, mode=:control)
plt = plot()
scatter!(plt, [x[1] for x in Xpolicy], color=:green, markersize=10, label="policy")
scatter!(plt, [x[2] for x in Xpolicy], color=:green, markersize=10, markershape=:square, label=nothing)
scatter!(plt, [x[1] for x in Xcontrol], color=:yellow, markersize=7, label="control")
scatter!(plt, [x[2] for x in Xcontrol], color=:yellow, markersize=7, markershape=:square, label=nothing)
scatter!(plt, [x[1] for x in Xref], color=:black, markersize=3, label="reference")
scatter!(plt, [x[2] for x in Xref], color=:black, markersize=3, markershape=:square, label=nothing)

θ1, us1 = unpack_vars(vars1, N=N0, H=H0)
con_u = [zeros(nu) for i=1:H0-1]
for i = 1:H0-1
	x = Xcontrol[i]
	up = policy(env, x, θ1)
	uc = us1[i]
	con_u[i] = up - uc
end
norm([con_u...;], Inf)

X4 = policy_rollout_only(env, x0, splits0, vars1)
con_x = [zeros(nx) for i=1:N0-1]
for i = 1:N0-1
	xend = X4[i][end]
	xstart = X4[i+1][1]
	con_x[i] = xstart - xend
end
con_x
norm([con_x...;], Inf)

obj0


function continuation_solve(env, obj, Xref, vars, H, N)
	x0 = Xref[1]
	for Ni in [20, 10, 5, 2]
		splits = split_trajectory(H, Ni)
		obj.ρx = 1e0
		obj.ρu = 1e0
		obj.λx = [zeros(nx) for i=1:Ni-1]
		obj.λu = [zeros(nu) for i=1:H-1]
		obj.Q *= 3

		θ, us = unpack_vars(vars, N=Ni, H=H)
		# vars = [θ; vcat([Xref[s[end]] for s in splits[1:Ni-1]]...); us...]
		vars = [θ; us...]
		vars = augmented_lagrangian_solver(env, obj, Xref, splits, vars)
	end
	return vars
end


vars0 = [0;0;0; vcat([Xref[s[end]] for s in splits0[1:N0-1]]...); zeros(nu*(H0-1))]
vars1 = continuation_solve(env, obj0, Xref, vars0, H0, N0)
vars1
splits1 = split_trajectory(H0, 2)
X2 = pure_policy_rollout(env, x0, splits1, vars1, mode=:control)
X3 = pure_policy_rollout(env, x0, splits1, vars1, mode=:policy)
plt = plot()
scatter!(plt, [x[1] for x in X2], color=:green, markersize=10, label="policy")
scatter!(plt, [x[2] for x in X2], color=:green, markersize=10, markershape=:square, label=nothing)
scatter!(plt, [x[1] for x in X2], color=:yellow, markersize=7, label="control")
scatter!(plt, [x[2] for x in X2], color=:yellow, markersize=7, markershape=:square, label=nothing)
scatter!(plt, [x[1] for x in Xref], color=:black, markersize=3, label="reference")
scatter!(plt, [x[2] for x in Xref], color=:black, markersize=3, markershape=:square, label=nothing)
