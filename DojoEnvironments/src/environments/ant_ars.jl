mutable struct AntARS{T,N} <: Environment{T,N}
    mechanism::Mechanism{T}
    storage::Storage{T,N}
end

function ant_ars(;
    horizon=100,
    timestep=0.05, 
    input_scaling=timestep,
    gravity=-9.81, 
    urdf=:ant,
    springs=0, 
    dampers=0, 
    parse_springs=true, 
    parse_dampers=true, 
    joint_limits=Dict([
        (:hip_1, [-30,30] * π / 180), 
        (:ankle_1, [30,70] * π / 180), 
        (:hip_2, [-30,30] * π / 180), 
        (:ankle_2, [-70,-30] * π / 180), 
        (:hip_3, [-30,30] * π / 180), 
        (:ankle_3, [-70,-30] * π / 180), 
        (:hip_4, [-30,30] * π / 180), 
        (:ankle_4, [30,70] * π / 180)]),
    keep_fixed_joints=true, 
    friction_coefficient=0.5,
    contact_feet=true, 
    contact_body=true,
    T=Float64)

    mechanism = get_ant(;
        timestep, 
        input_scaling,
        gravity, 
        urdf,
        springs, 
        dampers, 
        parse_springs, 
        parse_dampers, 
        joint_limits,
        keep_fixed_joints, 
        friction_coefficient,
        contact_feet, 
        contact_body,
        T
    )

    storage = Storage(horizon, length(mechanism.bodies))

    return AntARS{T,horizon}(mechanism, storage)
end

function DojoEnvironments.state_map(::AntARS, state)
    state = state[1:28]
    return state
end

function DojoEnvironments.input_map(::AntARS, input::AbstractVector)
    input = [zeros(6);input] # floating base not actuated
    return input
end

function Dojo.step!(environment::AntARS, state, input=nothing; k=1, record=false, opts=SolverOptions())
    state = state_map(environment, state)
    input = input_map(environment, input)
    Dojo.step_minimal_coordinates!(environment.mechanism, state, input; opts)
    record && Dojo.save_to_storage!(environment.mechanism, environment.storage, k)

    return
end

function DojoEnvironments.get_state(environment::AntARS{T}) where T
    contact_force = T[]
    for contact in environment.mechanism.contacts
        push!(contact_force, max(-1, min(1, contact.impulses[2][1])))
    end
    state = [get_minimal_state(environment.mechanism); contact_force]

    return state
end