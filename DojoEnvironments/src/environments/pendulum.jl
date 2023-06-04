mutable struct Pendulum{T,N} <: Environment{T,N}
    mechanism::Mechanism{T}
    storage::Storage{T,N}
end

function pendulum(;
    horizon=100,
    timestep=0.01,
    input_scaling=timestep, 
    gravity=-9.81,
    mass=1,
    link_length=1,
    color=RGBA(1, 0, 0),
    springs=0,
    dampers=0,
    joint_limits=Dict(),
    spring_offset=szeros(1),
    orientation_offset=one(Quaternion),
    T=Float64)

    mechanism = get_pendulum(;
        timestep,
        input_scaling, 
        gravity,
        mass,
        link_length,
        color,
        springs,
        dampers,
        joint_limits,
        spring_offset,
        orientation_offset,
        T
    )

    storage = Storage(horizon, length(mechanism.bodies))

    return Pendulum{T,horizon}(mechanism, storage)
end

function DojoEnvironments.state_map(::Pendulum, state)
    return state
end

function DojoEnvironments.input_map(::Pendulum, input::AbstractVector)
    return input
end

function Dojo.step!(environment::Pendulum, state, input=nothing; k=1, record=false, opts=SolverOptions())
    state = state_map(environment, state)
    input = input_map(environment, input)
    Dojo.step_minimal_coordinates!(environment.mechanism, state, input; opts)
    record && Dojo.save_to_storage!(environment.mechanism, environment.storage, k)

    return
end

function DojoEnvironments.get_state(environment::Pendulum)
    state = get_minimal_state(environment.mechanism)
    return state
end