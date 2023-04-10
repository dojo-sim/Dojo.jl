mutable struct AntARS 
    mechanism
    storage
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
    limits=true,
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
        limits,
        joint_limits,
        keep_fixed_joints, 
        friction_coefficient,
        contact_feet, 
        contact_body,
        T
    )

    storage = Storage(horizon, length(mechanism.bodies))

    return AntARS(mechanism, storage)
end

function state_map(::AntARS, x)
    return x
end

function input_map(::AntARS, u)
    u = [zeros(6);u] # floating base not actuated

    return u
end

function Dojo.step!(environment::AntARS, x, u; k=1, record=false, opts=SolverOptions())
    x = state_map(environment, x)
    u = input_map(environment, u)
    Dojo.step_minimal_coordinates!(environment.mechanism, x, u; opts)
    record && Dojo.save_to_storage!(environment.mechanism, environment.storage, k)

    return
end

# function Dojo.simulate!(environment::AntARS, controller! = (mechanism, k) -> nothing; kwargs...)
#     simulate!(environment.mechanism, 1:length(environment.storage), environment.storage, controller!; kwargs...)
# end

function get_state(environment::AntARS)
    contact_force = Float64[]
    for contact in environment.mechanism.contacts
        push!(contact_force, max(-1, min(1, contact.impulses[2][1])))
    end
    x = [get_minimal_state(environment.mechanism); contact_force]

    return x
end

function Dojo.visualize(environment::AntARS; kwargs...)
    Dojo.visualize(environment.mechanism, environment.storage; kwargs...)
end