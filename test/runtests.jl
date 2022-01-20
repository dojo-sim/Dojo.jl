using ForwardDiff
using FiniteDiff
using StaticArrays
using LinearAlgebra
using Random
using SparseArrays
using BenchmarkTools
using Dojo
using Test

include("integrator.jl")
include("minimal_coordinates.jl")
include("jacobian.jl")
include("momentum_conservation.jl")
include("energy_conservation.jl")
include("behaviors.jl")
include("joint_limits.jl")
