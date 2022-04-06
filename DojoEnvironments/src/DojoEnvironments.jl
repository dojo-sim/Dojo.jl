module DojoEnvironments

using LinearAlgebra
using Random

using MeshCat
using StaticArrays
using Dojo

import Dojo: string_to_symbol, add_limits, RotX, build_robot, initialize!

export
    Environment,
    dynamics,
    dynamics_jacobian_state,
    dynamics_jacobian_input,
    get_environment,
    get_mechanism,
    #step,
    get_observation,
    cost,
    is_done,
    #reset,
    render,
    seed,
    close,
    Space,
    BoxSpace

include("mechanisms.jl")
include("environment.jl")
include("dynamics.jl")
include("utilities.jl")
include("include.jl")

end # module
