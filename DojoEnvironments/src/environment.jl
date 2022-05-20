"""
    Environment{T} 

    simulation object containing a mechanism along with additional information useful for 
    reinforcement learning and trajectory optimization 

    mechanism: Mechanism
    representation: :minimal or :maximal state representation
    input_space: Space, limits on inputs
    observation_space: Space, limits on observations
    state: contains current minimal or maximal states
    dynamics_jacobian_state: dynamics Jacobian wrt to state
    dynamics_jacobian_input: dynamics Jacobian wrt to input
    input_previous: input applied to mechanism at previous time step
	control_map: mapping for inputs to mechanism translational and rotational dynamics
    num_states: dimension of minimal or maximal state
    num_inputs: dimension of inputs
    num_observations: dimension of observation
    info: object containing environment specific information
    rng: random number generator
    vis: Visualizer
    opts_step: SolverOptions
    opts_grad: SolverOptions
"""
mutable struct Environment{X,T,M,A,O,I}
    mechanism::M
    representation::Symbol
    input_space::A
    observation_space::O
    state::Vector{T}
    dynamics_jacobian_state::Matrix{T}
    dynamics_jacobian_input::Matrix{T}
    input_previous::Vector{T}
	control_map::Matrix{T}
    num_states::Int
    num_inputs::Int
    num_observations::Int
    info::I
    rng::Vector{MersenneTwister}
    vis::Visualizer
    opts_step::SolverOptions{T}
    opts_grad::SolverOptions{T}
end

# function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, environment::Environment)
#     summary(io, environment)
#     println(io, "")
#     println(io, " mechanism:               "*string(environment.mechanism))
#     println(io, " representation:          "*string(environment.representation))
#     println(io, " input_space:             "*string(environment.input_space))
#     println(io, " observation_space:       "*string(environment.observation_space))
#     println(io, " state:                   "*string(environment.state))
#     println(io, " dynamics_jacobian_state: "*string(environment.dynamics_jacobian_state))
#     println(io, " dynamics_jacobian_input: "*string(environment.dynamics_jacobian_input))
#     println(io, " input_previous:          "*string(environment.input_previous))
# 	println(io, " control_map:             "*string(environment.control_map))
#     println(io, " num_states:              "*string(environment.num_states))
#     println(io, " num_inputs:              "*string(environment.num_inputs))
#     println(io, " num_observations:        "*string(environment.num_observations))
#     println(io, " info:                    "*string(environment.info))
#     println(io, " rng:                     "*string(environment.rng))
#     println(io, " vis:                     "*string(environment.vis))
#     println(io, " opts_step:               "*string(environment.opts_step))
#     println(io, " opts_grad:               "*string(environment.opts_grad))
# end

"""
    get_environment(model; kwargs...)

    construct existing environment 

    model: name of of environment 
    kwargs: environment specific parameters
"""
function get_environment(model; kwargs...)
    return eval(string_to_symbol(model))(; string_to_symbol(kwargs)...)
end

"""
    step(env, x, u; gradients, attitude_decompress)

    simulates environment one time step 

    env: Environment 
    x: state 
    u: input 
    gradients: flag for computing gradients of dynamics
    attitude_decompress: flag for pre- and post-concatenating Jacobians with attitude Jacobians
"""
function Base.step(env::Environment, x, u;
    gradients=false,
    attitude_decompress=false)

    mechanism = env.mechanism
    timestep= mechanism.timestep

    x0 = x
    # u = clip(env.input_space, u) # control limits
    env.input_previous .= u  # for rendering in Gym
	u_scaled = env.control_map * u

    z0 = env.representation == :minimal ? minimal_to_maximal(mechanism, x0) : x0
    z1 = step!(mechanism, z0, u_scaled; opts=env.opts_step)
    env.state .= env.representation == :minimal ? maximal_to_minimal(mechanism, z1) : z1

    # Compute cost
    costs = cost(env, x, u)

	# Check termination
	done = is_done(env, x)

    # Gradients
    if gradients
        if env.representation == :minimal
            fx, fu = get_minimal_gradients!(env.mechanism, z0, u_scaled, opts=env.opts_grad)
        elseif env.representation == :maximal
            fx, fu = get_maximal_gradients!(env.mechanism, z0, u_scaled, opts=env.opts_grad)
            if attitude_decompress
                A0 = attitude_jacobian(z0, length(env.mechanism.bodies))
                A1 = attitude_jacobian(z1, length(env.mechanism.bodies))
                fx = A1 * fx * A0'
                fu = A1 * fu
            end
        end
        env.dynamics_jacobian_state .= fx
        env.dynamics_jacobian_input .= fu * env.control_map
    end

    info = Dict()
    return get_observation(env), -costs, done, info
end

function Base.step(env::Environment, u;
    gradients=false, 
    attitude_decompress=false) 
    step(env, env.state, u; 
        gradients=gradients, 
        attitude_decompress=attitude_decompress)
end

"""
    get_observation(env) 

    return observation for current state 

    env: Environment
"""
function get_observation(env::Environment)
    return env.state
end

"""
    cost(env, x, u) 

    return cost (-reward) for current state and input

    env: Environment
    x: state 
    u: input
"""
cost(env::Environment, x, u) = 0.0

"""
    is_done(env) 

    check for termination of simulation

    env: Environment
"""
is_done(env::Environment, x) = false

"""
    reset(env; x) 

    returns environment to nominal state

    env: Environment
    x: state
"""
function Base.reset(env::Environment{X};
    x=nothing) where X

    initialize!(env.mechanism, type2symbol(X))
    if x != nothing
        env.state = x
    else
        if env.representation == :minimal
            env.state .= get_minimal_state(env.mechanism)
        elseif env.representation == :maximal
            env.state .= get_maximal_state(env.mechanism)
        end
        env.input_previous .= 0.0
    end
    return get_observation(env)
end

function MeshCat.render(env::Environment, 
    mode="human")
    z = env.representation == :minimal ? minimal_to_maximal(env.mechanism, env.state) : env.state
    set_robot(env.vis, env.mechanism, z, name=:robot)
    return nothing
end

function seed(env::Environment; s=0)
    env.rng[1] = MersenneTwister(seed)
    return nothing
end

function Base.close(env::Environment; kwargs...)
    return nothing
end

"""
    Space{T,N} 

    Abstract type for domains
"""
abstract type Space{T,N} end

""" 
    BoxSpace{T,N} <: Environment{T,N}

    domain with lower and upper limits 

    n: dimension of domain 
    low:: lower limit 
    high: upper limit 
    shape: tuple (n,)
    dtype: type for domain
"""
mutable struct BoxSpace{T,N} <: Space{T,N}
    n::Int # box dimension
    low::AbstractVector{T} # minimum value
    high::AbstractVector{T} # maximum value
    shape::Tuple{Int} # this is always (n,), it's needed to interface with Stable-Baselines
    dtype::DataType # this is always T, it's needed to interface with Stable-Baselines
end

# function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, space::BoxSpace)
#     summary(io, space)
#     println(io, "")
#     println(io, " n:     "*string(space.n))
#     println(io, " low:   "*string(space.law))
#     println(io, " high:  "*string(space.high))
#     println(io, " shape: "*string(space.shape))
#     println(io, " dtype: "*string(space.dtype))
# end

function BoxSpace(n::Int; low::AbstractVector = -ones(n), high::AbstractVector = ones(n))
    T = promote_type(eltype.((low, high))...)
    return BoxSpace{T,n}(n, low, high, (n,), T)
end

function sample(s::BoxSpace{T,N}) where {T,N}
    return rand(T,N) .* (s.high .- s.low) .+ s.low
end

function contains(s::BoxSpace{T,N}, v::AbstractVector{T}) where {T,N}
    all(v .>= s.low) && all(v .<= s.high)
end

function clip(s::BoxSpace, u)
    clamp.(u, s.low, s.high)
end



