function unpack_vars(vars; N=N)
	nθ = (nx+1) * nu
	θ = vars[1:nθ]
	xstarts = [vars[nθ + nx*(i-1) .+ (1:nx)] for i=1:N-1]
	return θ, xstarts
end

function pack_vars(θ, xstarts; N=N)
	nθ = nu * nx
	vars = vcat(θ, xstarts...)
	return vars
end
# vars = rand(nva)
# vars - pack_vars(unpack_vars(vars)...)

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
	FiniteDiff.finite_difference_jacobian(θ -> policy(env, x, θ), θ)
end

function policy_jacobian_state(env::Environment, x, θ)
	FiniteDiff.finite_difference_jacobian(x -> policy(env, x, θ), x)
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




mutable struct AugmentedObjective115{T}
	Q::AbstractMatrix{T}
	Qθ::AbstractMatrix{T}
    Qstart::AbstractMatrix{T}
    R::AbstractMatrix{T}
	λ::Vector{Vector{T}}
	ρ::T
end

function constraints(X, N)
	con = [zeros(nx) for i=1:N-1]
	for i = 1:N-1
		xend = X[i][end]
		xstart = X[i+1][1]
		con[i] = xstart - xend
	end
	return con
end

function violation(X, N)
	con = constraints(X, N)
	(N == 1) && return 1e-3
	return norm([con...], Inf)
end
