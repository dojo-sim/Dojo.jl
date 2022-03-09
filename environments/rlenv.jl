using ReinforcementLearningBase: RLBase

mutable struct DojoRLEnv <: RLBase.AbstractEnv
    dojoenv
    action_space
    observation_space
    state
    reward
    done::Bool
    info::Dict
end

function DojoRLEnv(dojoenv::Environment)
    action_space = convert(RLBase.Space, dojoenv.input_space)
    observation_space = convert(RLBase.Space, dojoenv.observation_space)
    state = reset(dojoenv)
    return DojoRLEnv(dojoenv, action_space, observation_space, state, 0.0, false, Dict())
end

RLBase.action_space(env::DojoRLEnv) = env.action_space
RLBase.state_space(env::DojoRLEnv) = env.observation_space
RLBase.is_terminated(env::DojoRLEnv) = env.done

RLBase.reset!(env::DojoRLEnv) = reset(env.dojoenv)

RLBase.reward(env::DojoRLEnv) = error()
RLBase.state(env::DojoRLEnv) = env.state

Random.seed!(env::DojoRLEnv, seed) = Dojo.seed(env.dojoenv, seed)

function (env::DojoRLEnv)(a)
    s, r, d, i = step(env.dojoenv, a)
    env.state = s
    env.reward = r
    env.done = d
    env.info = i
    return nothing
end