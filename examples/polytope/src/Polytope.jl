module Polytope

using BenchmarkTools
using FiniteDiff
using GeometryBasics
using JLD2
using LinearAlgebra
using MeshCat
using Meshing
using Parameters
using Pkg
using Plots
using Quaternions
using Random
using StaticArrays

global const OSF_PATH = joinpath("/home/simon/research/repos/osf-pytorch")
# global OSF_PATH = joinpath("/home/simon/research/repos/osf-pytorch")
# include("setup.jl")
using PyCall

include("collider.jl")
include("utils.jl")
include("visualize.jl")
include("integrator.jl")
end # module
