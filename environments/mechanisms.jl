"""
    get_mechanism(model; kwargs...)
    
    constructs mechanism

    model: name of mechanism 
    kwargs: mechanism specific parameters
"""
function get_mechanism(model; kwargs...)
    # convert potential strings to symbols    
    mech = eval(Symbol(:get, :_, string_to_symbol(model)))(; string_to_symbol(kwargs)...)
    return mech
end

"""
    initialize!(mechanism, model; kwargs)

    state initialization for mechanism 

    mechanism: Mechanism 
    model: name of mechanism 
    kwargs: mechanism specific parameters
"""
function initialize!(mechanism::Mechanism, model; kwargs...)
    eval(Symbol(:initialize, :_, string_to_symbol(model), :!))(mechanism; kwargs...)
end

