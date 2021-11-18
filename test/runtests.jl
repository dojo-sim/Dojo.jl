using Test
using ForwardDiff
using FiniteDiff
using StaticArrays
using LinearAlgebra
using Random
using SparseArrays
using BenchmarkTools
using Dojo

# const Dojo = Main

include("jacobian.jl")
include("momentum_conservation.jl")
include("energy_conservation.jl")

# @safetestset "Quaternion Tests" begin
#     include("quaternion_test.jl")
# end

# include("diff/diff_test.jl")

# @safetestset "Factorization Test" begin
#     include("factorization_test.jl")
# end

# include("initialization/initialization_test.jl")

# @safetestset "UI Test" begin
#     include("ui_test.jl")
# end

# @safetestset "URDF Tests" begin
#     include("urdf_test.jl")
# end

# @safetestset "Optionals Tests" begin
#     include("optionals_test.jl")
# end

# @safetestset "Dynamics Tests" begin
#     include("dynamics/dynamics_test.jl")
# end

# @safetestset "Example Tests" begin
#     include("example_test.jl")
# end

# using ConstrainedDynamics
# using LinearAlgebra
# @testset "jacobians" begin
#     include("jacobian_test.jl")
# end
