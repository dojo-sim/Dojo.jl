using Test
using ForwardDiff
using FiniteDiff
using StaticArrays
using LinearAlgebra
using Random
using SparseArrays
using BenchmarkTools
using Dojo

include("integrator.jl")
include("minimal_coordinates.jl")
include("jacobian.jl")
include("momentum_conservation.jl")
include("energy_conservation.jl")
include("behavior.jl")
include("joint_limits.jl")
