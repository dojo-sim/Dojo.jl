mutable struct CartpoleDQN{T,N} <: Environment{T,N}
    mechanism::Mechanism{T}
    storage::Storage{T,N}
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

    return CartpoleDQN{T,horizon}(mechanism, storage)
end

function state_map(::CartpoleDQN, state)
    return state
end

function input_map(::CartpoleDQN, input)
    input = [input;0] # only the cart is actuated
    return input
end

function input_map(::CartpoleDQN, ::Nothing)
    input = zeros(2)
    return input
end

function Dojo.step!(environment::CartpoleDQN, state, input=nothing; k=1, record=false, opts=SolverOptions())
    state = state_map(environment, state)
    input = input_map(environment, input)
    Dojo.step_minimal_coordinates!(environment.mechanism, state, input; opts)
    record && Dojo.save_to_storage!(environment.mechanism, environment.storage, k)

    return
end

function get_state(environment::CartpoleDQN)
    state = get_minimal_state(environment.mechanism)
    return state
end