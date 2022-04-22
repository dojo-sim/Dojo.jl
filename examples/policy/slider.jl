using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using Dojo
using Plots

vis = Visualizer()
open(vis)

include("utils.jl")

################################################################################
# build environment
################################################################################
H = 25
N = 5
timestep = 0.05
gravity = [0, 0, 1.0]
spring = 0.0
damper = 0.0
env = slider(timestep=timestep, gravity=gravity, spring=spring, damper=damper, vis=vis)
nx = minimal_dimension(env.mechanism)
nu = input_dimension(env.mechanism)
nθ = (nx+1) * nu
nva = nθ + nx*(N-1)


################################################################################
# reference trajectory
################################################################################
initialize!(env.mechanism, :slider, position=0.50)
ref_traj = simulate!(env.mechanism, H*timestep)
visualize(env.mechanism, ref_traj, vis=vis)
Zref = get_maximal_state(ref_traj)
Xref = [maximal_to_minimal(env.mechanism, Zref[i]) for i=1:H]
x0 = Xref[1]
xF = Xref[H]


################################################################################
# multiple shooting
################################################################################
x = get_minimal_state(env.mechanism)
u = rand(nu)
y = zeros(nx)
dx = zeros(nx,nx)
du = zeros(nx,nu)
dynamics_and_jacobian(y, dx, du, env, x, u, nothing;
	attitude_decompress=false)

function policy_rollout(env::Environment, x0, vars::Vector)
	θ, xstarts = unpack_vars(vars)
	splits = split_trajectory(H, N)

	X = [[zeros(nx) for j=1:length(splits[i])] for i=1:N]
	dX = [[zeros(nx,nx) for j=1:length(splits[i])] for i=1:N]
	dU = [[zeros(nx,nu) for j=1:length(splits[i])] for i=1:N]

	# Initial states
	X[1][1] .= deepcopy(x0)
	for i = 1:N-1
		X[i+1][1] .= deepcopy(xstarts[i])
	end
	for i = 1:N
		for (j,ind) in enumerate(splits[i][1:end-1])
			x = X[i][j]
			y = X[i][j+1]
			dx = dX[i][j+1]
			du = dU[i][j+1]
			xref = Xref[ind]
			u = policy(env, x, θ)
			dynamics_and_jacobian(y, dx, du, env, x, u, nothing;
				attitude_decompress=false)
		end
	end
    return X, dX, dU
end

x0
splits0 = split_trajectory(H, N)
vars0 = [0;-10;-1;[1,-1];[2,-2];[3,-3];[4,-4]]
X0, dX0, dU0 = policy_rollout(env, [0,0], vars0)
plt = plot()
for i = 1:N
	scatter!(plt, splits0[i], [x[1] for x in X0[i]])
end
display(plt)

function cost(env::Environment, obj, Xref, vars)
	X, dX, dU = policy_rollout(env, Xref[1], vars)
	cost(env, obj, X, dX, dU, Xref, vars)
end

function cost_only(env::Environment, obj, Xref, vars)
	X, dX, dU = policy_rollout(env, Xref[1], vars)
	cost_only(env, obj, X, dX, dU, Xref, vars)
end

function cost(env::Environment, obj, X, dX, dU, Xref, vars)
	splits = split_trajectory(H, N)
	θ, xstarts = unpack_vars(vars)

	# Initialization
	f = 0.0
	g_θ = zeros(nθ)
	g_xstarts = [zeros(nx) for i=1:N-1]

	## Evaluation
	# Objective
	for i = 1:N
		for (j,ind) in enumerate(splits[i])
			x = X[i][j]
			xref = Xref[ind]
			f += 0.5 * (x - xref)' * obj.Q * (x - xref)
		end
	end
	# Constraints
	for i = 1:N-1
		λ = obj.λ[i]
		xend = X[i][end]
		xstart = X[i+1][1]
		f += λ' * (xstart - xend)
		f += 0.5 * obj.ρ * (xstart - xend)' * (xstart - xend)
	end

	## Gradient
	for i = 1:N-1
		λ = obj.λ[i]
		ind = splits[i+1][i] # starting index
		xref = Xref[ind]
		xend = X[i][end]
		xstart = X[i+1][1]
		# Objective
		g_xstarts[i] += obj.Q * (xstart - xref)
		# Constraints
		g_xstarts[i] += λ + obj.ρ * (xstart - xend)
	end
	for i = 1:N
		∂x∂θ = zeros(nx, nθ)
		for (j,ind) in enumerate(splits[i][1:end-1])
			x0 = X[i][j] # x0
			x1 = X[i][j+1] # x1
			x1ref = Xref[ind+1]

			u = policy(env, x0, θ)
			∂u∂θ = policy_jacobian_parameters(env, x0, θ)
			∂u∂x0 = policy_jacobian_state(env, x0, θ)

			∂x1∂u = dU[i][j+1]
			∂x1∂x0 = dX[i][j+1]

			∂x1∂θ = ∂x1∂u * ∂u∂θ
			∂x1∂θ += (∂x1∂x0 + ∂x1∂u * ∂u∂x0) * ∂x∂θ

			∂o∂x1 = (x1 - x1ref)' * obj.Q
			∂o∂u = u' * obj.R * u
			g_θ += (∂o∂x1 * ∂x1∂θ + ∂o∂u * ∂u∂θ)'

			# For the last element we have a stitching constraint, except for the last split
			if j == length(splits[i])-1 && i < N
				xend = x1
				xstart = X[i+1][1]
				λ = obj.λ[i]
				∂xend∂θ = ∂x1∂θ
				∂c∂xend = -λ' + obj.ρ * (xend - xstart)' # xend = x1
				g_θ += (∂c∂xend * ∂xend∂θ)'
			end
			# update
			∂x∂θ = ∂x1∂θ
		end
	end
	g_θ = g_θ[:,1]
	g = pack_vars(g_θ, g_xstarts)
	return f, g
end

function cost_only(env::Environment, obj, X, dX, dU, Xref, vars)
	splits = split_trajectory(H, N)
	θ, xstarts = unpack_vars(vars)

	# Initialization
	f = 0.0

	## Evaluation
	# Objective
	for i = 1:N
		for (j,ind) in enumerate(splits[i])
			x = X[i][j]
			xref = Xref[ind]
			f += 0.5 * (x - xref)' * obj.Q * (x - xref)
		end
	end
	# Constraints
	for i = 1:N-1
		λ = obj.λ[i]
		xend = X[i][end]
		xstart = X[i+1][1]
		f += λ' * (xstart - xend)
		f += 0.5 * obj.ρ * (xstart - xend)' * (xstart - xend)
	end
	return f
end

mutable struct AugmentedObjective114{T}
    Q::AbstractMatrix{T}
    R::AbstractMatrix{T}
	λ::Vector{Vector{T}}
	ρ::T
end


vars0 = zeros(nva)
θ, xs = unpack_vars(vars0)
X0, dX0, dU0 = policy_rollout(env, x0, vars0)

Q0 = 0*1e-0*Diagonal(ones(nx))
R0 = 0*1e-2*Diagonal(ones(nu))
λ0 = [(i in (1,2))*ones(nx) for i=1:N-1]
ρ0 = 0.0
obj0 = AugmentedObjective114(Q0, R0, λ0, ρ0)

cost(env, obj0, X0, dX0, dU0, Xref, vars0)
cost(env, obj0, Xref, vars0)

g1 = FiniteDiff.finite_difference_gradient(vars -> cost_only(env, obj0, Xref, vars), vars0)
f0, g0 = cost(env, obj0, Xref, vars0)
g1
g0
g0 - g1
norm(g0 - g1)

function constraints(X)
	con = [zeros(nx) for i=1:N-1]
	for i = 1:N-1
		xend = X[i][end]
		xstart = X[i+1][1]
		con[i] = xstart - xend
	end
	return con
end

function violation(X)
	con = constraints(X)
	return norm([con...], Inf)
end


function augmented_lagrangian_solver(env, obj, Xref, vars)
	x0 = Xref[1]
	X, _, _ = policy_rollout(env, x0, vars)
	for i = 1:6
		(violation(X) < 1e-3) && break
		for i = 1:20
			f = cost_only(env, obj, Xref, vars)
			g = FiniteDiff.finite_difference_gradient(vars -> cost_only(env, obj, Xref, vars), vars)
			(norm(g, Inf) < 1e-3) && break
			H = FiniteDiff.finite_difference_hessian(vars -> cost_only(env, obj, Xref, vars), vars)
			Δ = - H \ g
			α = linesearch(vars, Δ, g)
			vars += α * Δ
		end

		# Dual ascent and penalty update
		X, _, _ = policy_rollout(env, x0, vars)
		con = constraints(X)
		for i = 1:N-1
			obj.λ += ρ * con[i]
		end
		obj.ρ = min(obj.ρ * 10, 1e6)
	end

	return vars
end
