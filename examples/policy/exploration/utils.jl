#################################################################################\
# variables
#################################################################################\
function unpack_vars(vars; N=N, H=H)
	nθ = (nx+1) * nu
	θ = vars[1:nθ]
	# xstarts = [vars[nθ + nx*(i-1) .+ (1:nx)] for i=1:N-1]
	us = [vars[nθ + (i-1)*nu .+ (1:nu)] for i=1:H-1]
	# return θ, xstarts, us
	return θ, us
end

function pack_vars(θ, us; N=N, H=H)
	# vars = vcat(θ, xstarts..., us...)
	vars = vcat(θ, us...)
	return vars
end
# vars = rand(nva)
# norm(vars - pack_vars(unpack_vars(vars, N=N0, H=H0)..., N=N0, H=H0))

function split_trajectory(H::Int, N::Int)
	C = H + N - 1
	S = Int(floor(C/N))
	splits = [1:S]
	C -= S
	for i = 2:N
		S = Int(floor(C / (N-i+1)))
		push!(splits, splits[end][end]-1 .+ (1:S))
		C -= S
	end
    return splits
end


################################################################################
# Policy
################################################################################
function policy(env::Environment, x, θ)
	mechanism = env.mechanism
	timestep = mechanism.timestep
	nx = minimal_dimension(mechanism)
	nu = input_dimension(mechanism)
	θmat = reshape(θ, (nu,nx+1))
    u = θmat * [1; x]
    return u
end

function policy_jacobian_parameters(env::Environment, x, θ)
	# FiniteDiff.finite_difference_jacobian(θ -> policy(env, x, θ), θ)
	ForwardDiff.jacobian(θ -> policy(env, x, θ), θ)
end

function policy_jacobian_state(env::Environment, x, θ)
	# FiniteDiff.finite_difference_jacobian(x -> policy(env, x, θ), x)
	ForwardDiff.jacobian(x -> policy(env, x, θ), x)
end

function ctrl!(mechanism::Mechanism, k::Int; θ=nothing)
    timestep = mechanism.timestep
    nu = input_dimension(mechanism)
    nx = minimal_dimension(mechanism)

    x = get_minimal_state(mechanism)
    if θ != nothing
        u = θ * x
    else
        u = szeros(nu)
    end
    set_input!(mechanism, timestep * u)

    return nothing
end


################################################################################
# Plot
################################################################################
function plot_rollout(X, Xref, splits)
	N = length(splits)
	colors = [:red, :green, :orange, :blue, :yellow]
	plt = plot(legend=false)
	for i = 1:N
		scatter!(plt, splits[i], [x[1] for x in X[i]], color=colors[1+i%5], markersize=10)
		scatter!(plt, splits[i], [x[2] for x in X[i]], color=colors[1+i%5], markersize=10, markershape=:square)
	end
	scatter!(plt, [x[1] for x in Xref], color=:black)
	scatter!(plt, [x[2] for x in Xref], color=:black, markershape=:square)
	display(plt)
	return nothing
end

mutable struct AugmentedObjective120{T}
	Q::AbstractMatrix{T}
	Qθ::AbstractMatrix{T}
    Qstart::AbstractMatrix{T}
	R::AbstractMatrix{T}
	λx::Vector{Vector{T}}
	λu::Vector{Vector{T}}
	ρx::T
	ρu::T
	u_scale::T
end

function constraints(X, obj, splits, vars)
	N = length(splits)
	H = splits[end][end]
	θ, us = unpack_vars(vars, N=N, H=H)
	con_x = [zeros(nx) for i=1:N-1]
	for i = 1:N-1
		xend = X[i][end]
		xstart = X[i+1][1]
		con_x[i] = xstart - xend
	end
	con_u = [zeros(nu) for i=1:H-1]
	for i = 1:N
		for (j,ind) in enumerate(splits[i])
			if ind < H
				x = X[i][j]
				up = policy(env, x, θ)
				uc = us[ind]
				con_u[ind] = obj.u_scale * (up - uc)
			end
		end
	end
	return con_x, con_u
end

function violation(X, obj, splits, vars)
	con_x, con_u = constraints(X, obj, splits, vars)
	return norm([con_x...; con_u...], Inf)
	# @warn "ignore con_u"
	# return norm([con_x...], Inf)
end
