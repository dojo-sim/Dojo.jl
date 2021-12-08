"""
     Mechanism constructor. Provides a simple way to construct a example
     mechanisms: atlas, snake, box, etc.
"""
function getmechanism(model::Symbol; kwargs...)
    mech = eval(Symbol(:get, model))(; kwargs...)
    return mech
end

"""
     Mechanism initialization method. Provides a simple way to set the initial
     conditions (pose and velocity) of the mechanism.
"""
function initialize!(mechanism::Mechanism, model::Symbol; kwargs...)
    eval(Symbol(:initialize, model, :!))(mechanism; kwargs...)
end

include("atlas/methods/initialize.jl")

include("box/methods/initialize.jl")

include("cartpole/methods/initialize.jl")

include("dzhanibekov/methods/initialize.jl")

include("halfcheetah/methods/initialize.jl")

include("hopper/methods/initialize.jl")

include("humanoid/methods/initialize.jl")

include("orbital/methods/initialize.jl")

include("pendulum/methods/initialize.jl")

include("quadruped/methods/initialize.jl")

include("slider/methods/initialize.jl")

include("snake/methods/initialize.jl")

include("twister/methods/initialize.jl")

# Utilities
include("../examples/trajectory_optimization/utils.jl")

# # Templates
# include("atlas/methods/template.jl")
# include("halfcheetah/methods/template.jl")
# include("quadruped/methods/template.jl")




