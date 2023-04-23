# Running from within Dojo repo does not work somehow; create temp environment

# Copied from https://github.com/JuliaReinforcementLearning/ReinforcementLearning.jl/blob/main/src/ReinforcementLearningEnvironments/src/environments/examples/CartPoleEnv.jl

# ### Setup
# PKG_SETUP
using ReinforcementLearning # v0.10.2
using Random
using ClosedIntervals
using StableRNGs: StableRNG
using Flux
using Flux: params, ignore
using Flux.Losses: huber_loss
using Dojo
using DojoEnvironments

# ### BEGIN Redefinition because of breaking changes in ReinforcementLearning.jl repo
function RLBase.update!(learner::BasicDQNLearner, batch::NamedTuple{SARTS})

    Q = learner.approximator
    γ = learner.γ
    loss_func = learner.loss_func

    s, a, r, t, s′ = send_to_device(device(Q), batch)
    a = CartesianIndex.(a, 1:length(a))

    gs = gradient(params(Q)) do
        q = Q(s)[a]
        q′ = vec(maximum(Q(s′); dims = 1))
        G = @. r + γ * (1 - t) * q′
        loss = loss_func(G, q)
        ignore() do
            learner.loss = loss
        end
        loss
    end

    update!(Q, gs)
end

RLBase.update!(app::NeuralNetworkApproximator, gs) =
    Flux.Optimise.update!(app.optimizer, params(app), gs)
# ### END Redefinition because of breaking changes in ReinforcementLearning.jl repo


struct CartPoleEnvParams{T}
    gravity::T
    masscart::T
    masspole::T
    totalmass::T
    halflength::T
    polemasslength::T
    forcemag::T
    dt::T
    thetathreshold::T
    xthreshold::T
    max_steps::Int
end

Base.show(io::IO, params::CartPoleEnvParams) = print(
    io,
    join(["$p=$(getfield(params, p))" for p in fieldnames(CartPoleEnvParams)], ","),
)

function CartPoleEnvParams{T}(;
    gravity=9.8,
    masscart=1.0,
    masspole=0.1,
    halflength=0.5,
    forcemag=10.0,
    max_steps=200,
    dt=0.02,
    thetathreshold=12.0,
    xthreshold=2.4
) where {T}
    CartPoleEnvParams{T}(
        gravity,
        masscart,
        masspole,
        masscart + masspole,
        halflength,
        masspole * halflength,
        forcemag,
        dt,
        thetathreshold * π / 180,
        xthreshold,
        max_steps,
    )
end

mutable struct CartPoleEnv{T,ACT} <: AbstractEnv
    dojo_env::DojoEnvironments.CartpoleDQN # Dojo Environment
    record::Bool # Select if trajectory should be recorded in Dojo Storage
    params::CartPoleEnvParams{T}
    state::Vector{T}
    action::ACT
    done::Bool
    t::Int
    rng::AbstractRNG
end

function CartPoleEnv(; T=Float64, record=false, continuous=false, rng=Random.GLOBAL_RNG, kwargs...)
    params = CartPoleEnvParams{T}(; kwargs...)

    ## create Dojo Environment
    dojo_env = get_environment(:cartpole_dqn;
        horizon=params.max_steps+1,
        timestep = params.dt,
        gravity = -params.gravity,
        slider_mass = params.masscart,
        pendulum_mass = params.masspole,
        link_length = params.halflength*2
    )

    env = CartPoleEnv(dojo_env, record, params, zeros(T, 4), continuous ? zero(T) : zero(Int), false, 0, rng)
    reset!(env)
    env
end

CartPoleEnv{T}(; kwargs...) where {T} = CartPoleEnv(T=T, kwargs...)

Random.seed!(env::CartPoleEnv, seed) = Random.seed!(env.rng, seed)
RLBase.reward(env::CartPoleEnv{T}) where {T} = env.done ? zero(T) : one(T)
RLBase.is_terminated(env::CartPoleEnv) = env.done
RLBase.state(env::CartPoleEnv) = env.state

function RLBase.state_space(env::CartPoleEnv{T}) where {T}
    ((-2 * env.params.xthreshold) .. (2 * env.params.xthreshold)) ×
    (typemin(T) .. typemax(T)) ×
    ((-2 * env.params.thetathreshold) .. (2 * env.params.thetathreshold)) ×
    (typemin(T) .. typemax(T))
end

RLBase.action_space(env::CartPoleEnv{<:AbstractFloat,Int}, player) = Base.OneTo(2)
RLBase.action_space(env::CartPoleEnv{<:AbstractFloat,<:AbstractFloat}, player) = -1.0 .. 1.0

function RLBase.reset!(env::CartPoleEnv{T}) where {T}
    env.state[:] = T(0.1) * rand(env.rng, T, 4) .- T(0.05)
    env.t = 0
    env.action = rand(env.rng, action_space(env))
    env.done = false
    nothing
end

function (env::CartPoleEnv)(a::AbstractFloat)
    @assert a in action_space(env)
    env.action = a
    _step!(env, a)
end

function (env::CartPoleEnv)(a::Int)
    @assert a in action_space(env)
    env.action = a
    _step!(env, a == 2 ? 1 : -1)
end

function _step!(env::CartPoleEnv, a)
    env.t += 1
    force = a * env.params.forcemag
    #= 
    Remove ReinforcementLearningEnvironments dynamics
    x, xdot, theta, thetadot = env.state
    costheta = cos(theta)
    sintheta = sin(theta)
    tmp = (force + env.params.polemasslength * thetadot^2 * sintheta) / env.params.totalmass
    thetaacc =
        (env.params.gravity * sintheta - costheta * tmp) / (
            env.params.halflength *
            (4 / 3 - env.params.masspole * costheta^2 / env.params.totalmass)
        )
    xacc = tmp - env.params.polemasslength * thetaacc * costheta / env.params.totalmass
    env.state[1] += env.params.dt * xdot
    env.state[2] += env.params.dt * xacc
    env.state[3] += env.params.dt * thetadot
    env.state[4] += env.params.dt * thetaacc
    =#

    ## Use Dojo dynamics
    DojoEnvironments.step!(env.dojo_env, env.state, force; k=env.t, env.record)
    env.state = DojoEnvironments.get_state(env.dojo_env)
    
    env.done =
        abs(env.state[1]) > env.params.xthreshold ||
        abs(env.state[3]) > env.params.thetathreshold ||
        env.t > env.params.max_steps
    nothing
end

# Copied from https://juliareinforcementlearning.org/docs/experiments/experiments/DQN/JuliaRL_BasicDQN_CartPole/#JuliaRL\\_BasicDQN\\_CartPole

env = CartPoleEnv(; record=true, rng = MersenneTwister(123))

seed = 123
rng = StableRNG(seed)
ns, na = length(state(env)), length(action_space(env))

policy = Agent(
    policy = QBasedPolicy(
        learner = BasicDQNLearner(
            approximator = NeuralNetworkApproximator(
                model = Chain(
                    Dense(ns, 128, relu; init = glorot_uniform(rng)),
                    Dense(128, 128, relu; init = glorot_uniform(rng)),
                    Dense(128, na; init = glorot_uniform(rng)),
                ) |> gpu,
                optimizer = ADAM(),
            ),
            batch_size = 32,
            min_replay_history = 100,
            loss_func = huber_loss,
            rng = rng,
        ),
        explorer = EpsilonGreedyExplorer(
            kind = :exp,
            ϵ_stable = 0.01,
            decay_steps = 500,
            rng = rng,
        ),
    ),
    trajectory = CircularArraySARTTrajectory(
        capacity = 1000,
        state = Vector{Float32} => (ns,),
    ),
);

# ## Train policy
run(
        policy,
        CartPoleEnv(),
        StopAfterEpisode(1000),
        TotalRewardPerEpisode()
    )

# ## Rollout policy once and record (record=true in env)
run(
        policy,
        env,
        StopAfterEpisode(1),
        TotalRewardPerEpisode()
    )

# ### Visualize trained cartpole
vis = visualize(env.dojo_env)
render(vis)