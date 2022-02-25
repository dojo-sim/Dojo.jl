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

include("ant/methods/initialize.jl")

include("atlas/methods/initialize.jl")

include("box/methods/initialize.jl")

include("box2d/methods/initialize.jl")

include("cartpole/methods/initialize.jl")

include("dzhanibekov/methods/initialize.jl")

include("fourbar/methods/initialize.jl")

include("halfcheetah/methods/initialize.jl")

include("hopper/methods/initialize.jl")

include("humanoid/methods/initialize.jl")

include("orbital/methods/initialize.jl")

include("pendulum/methods/initialize.jl")

include("quadruped/methods/initialize.jl")

include("raiberthopper/methods/initialize.jl")

include("rexhopper/methods/initialize.jl")

include("slider/methods/initialize.jl")

include("sphere/methods/initialize.jl")

include("snake/methods/initialize.jl")

include("tennisracket/methods/initialize.jl")

include("tippetop/methods/initialize.jl")

include("twister/methods/initialize.jl")

include("walker2d/methods/initialize.jl")
