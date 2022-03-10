using ReinforcementLearningBase: RLBase

mutable struct DojoRLEnv{T} <: RLBase.AbstractEnv
    dojoenv::Environment
    state::Vector{T}
    reward::T
    done::Bool
    info::Dict
end

function DojoRLEnv(dojoenv::Environment{X,T}) where {X,T}
    state = reset(dojoenv)
    return DojoRLEnv{T}(dojoenv, state, convert(T, 0.0), false, Dict())
end

function DojoRLEnv(name::String; kwargs...)
    DojoRLEnv(Dojo.get_environment(name; kwargs...))
end

RLBase.action_space(env::DojoRLEnv) = env.dojoenv.input_space
RLBase.state_space(env::DojoRLEnv) = env.dojoenv.observation_space
RLBase.is_terminated(env::DojoRLEnv) = env.done

RLBase.reset!(env::DojoRLEnv) = reset(env.dojoenv)

RLBase.reward(env::DojoRLEnv) = env.reward
RLBase.state(env::DojoRLEnv) = env.state

Random.seed!(env::DojoRLEnv, seed) = Dojo.seed(env.dojoenv, seed)

# TODO:
# RLBase.ChanceStyle(env::DojoRLEnv) = RLBase.DETERMINISTIC

function (env::DojoRLEnv)(a)
    s, r, d, i = step(env.dojoenv, a)
    env.state .= s
    env.reward = r
    env.done = d
    env.info = i
    return nothing
end
(env::DojoRLEnv)(a::AbstractFloat) = env([a])
