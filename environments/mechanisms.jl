"""
     Mechanism constructor. Provides a simple way to construct a example
     mechanisms: atlas, snake, box, etc.
"""
function get_mechanism(model::Symbol; kwargs...)
    mech = eval(Symbol(:get, :_, model))(; kwargs...)
    return mech
end

"""
     Mechanism initialization method. Provides a simple way to set the initial
     conditions (pose and velocity) of the mechanism.
"""
function initialize!(mechanism::Mechanism, model::Symbol; kwargs...)
    eval(Symbol(:initialize, :_, model, :!))(mechanism; kwargs...)
end

