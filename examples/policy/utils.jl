function unpack_vars(vars)
	nθ = (nx+1) * nu
	θ = vars[1:nθ]
	xstarts = [vars[nθ + nx*(i-1) .+ (1:nx)] for i=1:N-1]
	return θ, xstarts
end

function pack_vars(θ, xstarts)
	nθ = nu * nx
	vars = vcat(θ, xstarts...)
	return vars
end
# vars = rand(nva)
# vars - pack_vars(unpack_vars(vars)...)

function split_trajectory(H::Int, N::Int)
    S = Int(floor(H / N))
    rest = H % N
    splits = [S*(i-1) .+ (1:S+1) for i=1:N-1]
    push!(splits, splits[end][end]:H)
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
