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

function Base.convert(::Type{RLBase.Space}, s::BoxSpace)
    RLBase.Space([BoxSpace(1; low = s.low[i:i], high = s.high[i:i]) for i in 1:s.n])
end

RLBase.action_space(env::DojoRLEnv) = convert(RLBase.Space, env.dojoenv.input_space)
RLBase.state_space(env::DojoRLEnv) = convert(RLBase.Space, env.dojoenv.observation_space)
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
(env::DojoRLEnv)(a::Number) = env([a])
