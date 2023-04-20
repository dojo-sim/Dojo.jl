"""
    Environment{T,N} 

    abstract type for an environment consisting of a mechanism and a storage
"""
abstract type Environment{T,N} end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, environment::Environment)
    summary(io, environment)
end

"""
    get_environment(model; kwargs...)

    construct existing environment 

    model: name of of environment 
    kwargs: environment specific parameters
"""
function get_environment(model; kwargs...)
    return eval(string_to_symbol(model))(; string_to_symbol(kwargs)...)
end

"""
    state_map(environment, state)

    maps the provided state to the environments internal state

    environment: environment 
    state:       provided state
"""
function state_map(::Environment, state)
    return state
end

"""
    input_map(environment, input)

    maps the provided input to the environment's internal input

    environment: environment 
    input:       provided input
"""
function input_map(::Environment, input)
    return input
end

function input_map(environment::Environment, ::Nothing)
    return zeros(input_dimension(environment.mechanism))
end

"""
    set_input!(environment, input)

    sets the provided input to the environment's mechanism

    environment: environment 
    input:       provided input
"""
function Dojo.set_input!(environment::Environment, input)
    set_input!(environment.mechanism, input_map(environment, input))
    return
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
function Dojo.step!(environment::Environment, state, input=nothing; k=1, record=false, opts=SolverOptions())
    state = state_map(environment, state)
    input = input_map(environment, input)
    Dojo.step_minimal_coordinates!(environment.mechanism, state, input; opts)
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
function Dojo.simulate!(environment::Environment{T,N}, controller! = (environment, k) -> nothing; kwargs...) where {T,N}
    controller_wrapper!(mechanism, k) = controller!(environment, k)
    simulate!(environment.mechanism, 1:N, environment.storage, controller_wrapper!; kwargs...)
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
    initialize!(environment; kwargs...)

    initializes the environment's mechanism

    environment: Environment
    kwargs:      same as for DojoEnvironments' mechanisms 
"""
function Dojo.initialize!(environment::Environment, model; kwargs...)
    eval(Symbol(:initialize, :_, string_to_symbol(model), :!))(environment.mechanism; kwargs...)
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
