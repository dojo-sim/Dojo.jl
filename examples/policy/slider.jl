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
N0 = 20
timestep = 0.05
gravity = [0, 0, 0.0]
spring = 20.0
damper = 1.0
env = slider(timestep=timestep, gravity=gravity, spring=spring, damper=damper, vis=vis)
nx = minimal_dimension(env.mechanism)
nu = input_dimension(env.mechanism)
nθ = (nx+1) * nu
nva = nθ + nx*(N0-1) + nu*(H0-1)


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
vars0 = [0;0;0; vcat([Xref[s[end]] for s in splits0[1:N0-1]]...); zeros(nu*(H0-1))]
X0 = policy_rollout_only(env, x0, splits0, vars0)
plot_rollout(X0, Xref, splits0)


Q0 = 1*1e2*Diagonal(ones(nx))
Qθ0 = 1*1e1*Diagonal(ones(nx+1))
Qstart0 = 1*1e2*Diagonal(ones(nx))
R0 = 1*1e-3*Diagonal(ones(nu))
λx0 = [zeros(nx) for i=1:N0-1]
λu0 = [zeros(nu) for i=1:H0-1]
ρx0 = 1e-3
ρu0 = 1e-3
obj0 = AugmentedObjective119(Q0, Qθ0, Qstart0, R0, λx0, λu0, ρx0, ρu0)

vars1 = augmented_lagrangian_solver(env, obj0, Xref, splits0, vars0)
X2 = pure_policy_rollout(env, x0, splits0, vars1, mode=:control)
X3 = pure_policy_rollout(env, x0, splits0, vars1, mode=:policy)
plt = plot()
scatter!(plt, [x[1] for x in X2], color=:green, markersize=10, label="policy")
scatter!(plt, [x[2] for x in X2], color=:green, markersize=10, markershape=:square, label=nothing)
scatter!(plt, [x[1] for x in X2], color=:yellow, markersize=7, label="control")
scatter!(plt, [x[2] for x in X2], color=:yellow, markersize=7, markershape=:square, label=nothing)
scatter!(plt, [x[1] for x in Xref], color=:black, markersize=3, label="reference")
scatter!(plt, [x[2] for x in Xref], color=:black, markersize=3, markershape=:square, label=nothing)




function continuation_solve(env, obj, Xref, vars, H, N)
	x0 = Xref[1]
	for Ni in [20, 10, 5, 2]
		splits = split_trajectory(H, Ni)
		obj.ρx = 1e0
		obj.ρu = 1e0
		obj.λx = [zeros(nx) for i=1:Ni-1]
		obj.λu = [zeros(nu) for i=1:H-1]
		obj.Q *= 3

		θ, xstarts, us = unpack_vars(vars, N=Ni, H=H)
		vars = [θ; vcat([Xref[s[end]] for s in splits[1:Ni-1]]...); us...]
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
