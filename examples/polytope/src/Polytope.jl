module Polytope

using BenchmarkTools
using GeometryBasics
using GLPK
using JLD2
using LinearAlgebra
using MeshCat
using Meshing
using Pkg
using Plots
using Polyhedra
using Random
using StaticArrays

global OSF_PATH = joinpath("/home/simon/research/repos/osf-pytorch")
include("setup.jl")
using PyCall

include("utils.jl")
end # module
