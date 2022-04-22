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
H0 = 40
N0 = 20
timestep = 0.05
gravity = [0, 0, 0.0]
spring = 20.0
damper = 1.0
env = slider(timestep=timestep, gravity=gravity, spring=spring, damper=damper, vis=vis)
nx = minimal_dimension(env.mechanism)
nu = input_dimension(env.mechanism)
nθ = (nx+1) * nu
nva = nθ + nx*(N0-1)


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
vars0 = [0;0;0; vcat([Xref[s[end]] for s in splits0[1:N0-1]]...)]
X0 = policy_rollout_only(env, x0, splits0, vars0)
plot_rollout(X0, Xref, splits0)


Q0 = 1*1e2*Diagonal(ones(nx))
Qθ0 = 1*1e1*Diagonal(ones(nx+1))
Qstart0 = 1*1e3*Diagonal(ones(nx))
R0 = 1*1e-1*Diagonal(ones(nu))
λ0 = [zeros(nx) for i=1:N0-1]
ρ0 = 1e-3
obj0 = AugmentedObjective115(Q0, Qθ0, Qstart0, R0, λ0, ρ0)

vars1 = augmented_lagrangian_solver(env, obj0, Xref, splits0, vars0)
X1 = policy_rollout_only(env, x0, splits0, vars1)
plot_rollout(X1, Xref, splits0)

X2 = pure_policy_rollout(env, x0, splits0, vars1)
plt = plot(legend=false)
scatter!(plt, [x[1] for x in X2], color=:yellow, markersize=10)
scatter!(plt, [x[2] for x in X2], color=:yellow, markersize=10, markershape=:square)
scatter!(plt, [x[1] for x in Xref], color=:black)
scatter!(plt, [x[2] for x in Xref], color=:black, markershape=:square)


continuation_solve(env, obj0, Xref, vars0, H0, N0)
