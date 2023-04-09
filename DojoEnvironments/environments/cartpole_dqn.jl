mutable struct CartpoleDQN 
    mechanism
    storage
end

function cartpole_dqn(;
    horizon=100,
    timestep=0.01, 
    input_scaling=timestep, 
    gravity=-9.81, 
    slider_mass=1,
    pendulum_mass=1,
    pendulum_length=1,
    radius=0.075,
    color=RGBA(0.7, 0.7, 0.7, 1),
    springs=0, 
    dampers=0,
    limits=false,
    joint_limits=Dict(),
    keep_fixed_joints=false, 
    T=Float64)

    mechanism = get_cartpole(;
        timestep, 
        input_scaling, 
        gravity, 
        slider_mass,
        pendulum_mass,
        pendulum_length,
        radius,
        color,
        springs, 
        dampers,
        limits,
        joint_limits,
        keep_fixed_joints, 
        T
    )

    storage = Storage(horizon, length(mechanism.bodies))

    return CartpoleDQN(mechanism, storage)
end

function state_map(::CartpoleDQN, x)
    return x
end

function input_map(::CartpoleDQN, u)
    u = [u;0] # only the cart is actuated

    return u
end

function step!(environment::CartpoleDQN, x, u; k=1, record=false, opts=SolverOptions())
    x = state_map(environment, x)
    u = input_map(environment, u)
    Dojo.step_minimal_coordinates!(environment.mechanism, x, u; opts)
    record && Dojo.save_to_storage!(environment.mechanism, environment.storage, k)

    return
end

function get_state(environment::CartpoleDQN)
    z = Dojo.get_next_state(environment.mechanism)
    x = maximal_to_minimal(environment.mechanism, z)

    return x
end

function visualize(environment::CartpoleDQN; kwargs...)
    visualize(environment.mechanism, environment.storage; kwargs...)
end