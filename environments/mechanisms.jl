"""
    get_mechanism(model; kwargs...)
    
    constructs mechanism

    model: Symbol 
    kwargs: mechanism specific parameters
"""
function get_mechanism(model::Symbol; kwargs...)
    mech = eval(Symbol(:get, :_, model))(; kwargs...)
    return mech
end

"""
    initialize!(mechanism, model; kwargs)

    state initialization for mechanism 

    mechanism: Mechanism 
    model: Symbol 
    kwargs: mechanism specific parameters
"""
function initialize!(mechanism::Mechanism, model::Symbol; kwargs...)
    eval(Symbol(:initialize, :_, model, :!))(mechanism; kwargs...)
end

