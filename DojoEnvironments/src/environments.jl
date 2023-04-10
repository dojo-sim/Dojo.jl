"""
    Environment{T} 

    simulation object containing a mechanism and a storage

    mechanism: Mechanism
    storage:   Storage
"""
abstract type Environment{T,N} end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, environment::Environment)
    summary(io, environment)
    println(io, "")
    println(io, " mechanism: "*string(environment.mechanism))
    println(io, " storage:   "*string(environment.storage))
end

"""
    get_environment(model; kwargs...)

    construct existing environment 

    model: name of of environment 
    kwargs: environment specific parameters
"""
function get_environment(model; horizon=100, kwargs...)
    return eval(string_to_symbol(model))(; string_to_symbol(kwargs)...)
end

"""
    state_map(environment, state)

    maps the provided state to the environments internal state

    environment: environment 
    state:       provided state
"""
function state_map(::Environment, state)
    return x
end

"""
    input_map(environment, input)

    maps the provided input to the environments internal input

    environment: environment 
    input:       provided input
"""
function input_map(::Environment, u)
    return u
end

"""
    step(environment, x, u; k, record, opts)

    simulates environment one time step 

    environment: Environment 
    x:           state 
    u:           input 
    k:           current timestep
    record:      record step in storage
    opts:        SolverOptions
"""
function Dojo.step!(environment::Environment, x, u=nothing; k=1, record=false, opts=SolverOptions())
    x = state_map(environment, x)
    u = input_map(environment, u)
    Dojo.step_minimal_coordinates!(environment.mechanism, x, u; opts)
    record && Dojo.save_to_storage!(environment.mechanism, environment.storage, k)

    return
end

"""
    simulate!(environment, controller! = (mechanism, k) -> nothing; kwargs...)

    simulates the mechanism of the environment

    environment: Environment
    controller!: Control function
    kwargs:      same as for Dojo.simulate 
"""
function Dojo.simulate!(environment::Environment, controller! = (mechanism, k) -> nothing; kwargs...)
    simulate!(environment.mechanism, 1:length(environment.storage), environment.storage, controller!; kwargs...)
end

"""
    get_state(environment)

    returns the state of the environment

    environment: Environment 
"""
function get_state(environment::Environment)
    return get_minimal_state(environment.mechanism)
end

"""
    visualize(environment; kwargs...)

    visualizes the environment with its internal storage

    environment: Environment
    kwargs:      same as for Dojo.visualize 
"""
function Dojo.visualize(environment::Environment; kwargs...)
    visualize(environment.mechanism, environment.storage; kwargs...)
end
