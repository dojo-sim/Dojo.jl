using ReinforcementLearning
using Flux
using Flux.Losses

using Random
using Distributions
using Dojo

function RL.Experiment(
    ::Val{:JuliaRL},
    ::Val{:PPO},
    ::Val{:DojoAnt},
    ::Nothing,
    save_dir = nothing,
    seed = 42
)
    rng = MersenneTwister(seed)
    N_ENV = 6
    UPDATE_FREQ = 32
    env_vec = [Dojo.DojoRLEnv("ant") for i in 1:N_ENV]
    for i in 1:N_ENV
        Random.seed!(env_vec[i], hash(seed+i))
    end
    env = MultiThreadEnv(env_vec)

    ns, na = length(state(env[1])), length(action_space(env[1]))
    RLBase.reset!(env; is_force=true)
    
    agent = Agent(
        policy = PPOPolicy(
            approximator = ActorCritic(
                actor = GaussianNetwork(
                    pre = Chain(
                        Dense(ns, 64, relu; init = glorot_uniform(rng)),
                        Dense(64, 64, relu; init = glorot_uniform(rng)),
                    ),
                    μ = Chain(Dense(64, na, tanh; init = glorot_uniform(rng)), vec),
                    logσ = Chain(Dense(64, na; init = glorot_uniform(rng)), vec),
                ),                
                critic = Chain(
                    Dense(ns, 256, relu; init = glorot_uniform(rng)),
                    Dense(256, na; init = glorot_uniform(rng)),
                ),
                optimizer = ADAM(1e-3),
            ),
            γ = 0.99f0,
            λ = 0.95f0,
            clip_range = 0.1f0,
            max_grad_norm = 0.5f0,
            n_epochs = 4,
            n_microbatches = 4,
            actor_loss_weight = 1.0f0,
            critic_loss_weight = 0.5f0,
            entropy_loss_weight = 0.001f0,
            dist = Normal,
            update_freq = UPDATE_FREQ,
        ),    
        trajectory = PPOTrajectory(;
            capacity = UPDATE_FREQ,
            state = Matrix{Float32} => (ns, N_ENV),
            action = Vector{Float32} => (N_ENV,),
            action_log_prob = Vector{Float32} => (N_ENV,),
            reward = Vector{Float32} => (N_ENV,),
            terminal = Vector{Bool} => (N_ENV,),
        ),
    )        
    stop_condition = StopAfterStep(10_000, is_show_progress=!haskey(ENV, "CI"))
    hook = TotalBatchRewardPerEpisode(N_ENV)
    Experiment(agent, env, stop_condition, hook, "# PPO with Dojo Ant")        
end

ex = E`JuliaRL_PPO_DojoAnt`
run(ex)