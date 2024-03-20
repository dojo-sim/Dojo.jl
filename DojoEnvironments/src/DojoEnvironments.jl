module DojoEnvironments

using LinearAlgebra
using Random

using Dojo
import Dojo: string_to_symbol, add_limits, RotX, RotY, RotZ, rotation_vector, SVector

include("mechanisms.jl")
include("environments.jl")
include("utilities.jl")
include("mechanisms/include.jl")
include("environments/include.jl")

# Mechanism
export
    get_mechanism,
    set_springs!,
    set_dampers!,
    set_limits

# Environment
export
    Environment,
    get_environment,
    step!,
    get_state,
    visualize

end # module
