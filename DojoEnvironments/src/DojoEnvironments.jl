module DojoEnvironments

# constants
global const X_AXIS = [1;0;0]
global const Y_AXIS = [0;1;0]
global const Z_AXIS = [0;0;1]

using LinearAlgebra
using Random

using MeshCat
using StaticArrays
using Dojo

import Dojo: string_to_symbol, add_limits, RotX, RotY, RotZ, build_robot, initialize!, set_minimal_coordinates_velocities!, velocity_index, vector_rotate, Prototype, rotation_vector, gray_light

export
    Environment,
    dynamics,
    dynamics_jacobian_state,
    dynamics_jacobian_input,
    get_environment,
    get_mechanism,
    get_observation,
    cost,
    is_done,
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
