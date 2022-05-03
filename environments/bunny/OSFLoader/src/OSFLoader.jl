module OSFLoader

using FiniteDiff
using JLD2
using LinearAlgebra
using Pkg
using Plots
using PyCall
using StaticArrays
using Statistics
using MeshIO
using GeometryBasics

global OSF_PATH = joinpath("/home/simon/research/repos/osf-pytorch")
global OSF_CONFIG_FOLDER = joinpath("configs/osf")

function osf_loader_dir()
    joinpath("..", @__DIR__)
end

include("load.jl")
include("utils.jl")

export
    load_density_script

export
    get_nerf_object,
    density_query,
    grid_particles,
    slice_particles,
    grid_density,
    slice_density,
    inertia_properties,
    sample_soft,
    finite_difference_gradient

end # module
