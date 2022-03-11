using ReinforcementLearning
using Flux
using Flux.Losses

using Random
using Dojo

function RL.Experiment(
    ::Val{:JuliaRL},
    ::Val{:DDPG},
    ::Val{:DojoAnt},
    ::Nothing,
    save_dir = nothing,
    seed = 42
)

    rng = MersenneTwister(seed)
    env = Dojo.DojoRLEnv("ant")
    Random.seed!(env, seed)
    A = action_space(env)
    ns, na = length(state(env)), length(action_space(env))
    @show na
    
    init = glorot_uniform(rng)
    
    create_actor() = Chain(
        Dense(ns, 30, relu; init = init),
        Dense(30, 30, relu; init = init),
        Dense(30, na, tanh; init = init),
    )
    create_critic() = Chain(
        Dense(ns + na, 30, relu; init = init),
        Dense(30, 30, relu; init = init),
        Dense(30, 1; init = init),
    )

    agent = Agent(
        policy = DDPGPolicy(
            behavior_actor = NeuralNetworkApproximator(
                model = create_actor(),
                optimizer = ADAM(),
            ),
            behavior_critic = NeuralNetworkApproximator(
                model = create_critic(),
                optimizer = ADAM(),
            ),
            target_actor = NeuralNetworkApproximator(
                model = create_actor(),
                optimizer = ADAM(),
            ),
            target_critic = NeuralNetworkApproximator(
                model = create_critic(),
                optimizer = ADAM(),
            ),
            γ = 0.99f0,
            ρ = 0.995f0,
            na = na,
            batch_size = 64,
            start_steps = 1000,
            start_policy = RandomPolicy(A; rng = rng),
            update_after = 1000,
            update_freq = 1,
            act_limit = 1.0,
            act_noise = 0.1,
            rng = rng,
        ),
        trajectory = CircularArraySARTTrajectory(
            capacity = 10000,
            state = Vector{Float32} => (ns,),
            action = Float32 => (na, ),
        ),
    )

    stop_condition = StopAfterStep(10_000, is_show_progress=!haskey(ENV, "CI"))
    hook = TotalRewardPerEpisode()
    Experiment(agent, env, stop_condition, hook, "# Dojo Ant with DDPG")    
end

ex = E`JuliaRL_DDPG_DojoAnt`
run(ex)