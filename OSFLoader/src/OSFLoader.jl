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
using Meshing
using MeshCat
using GeometryBasics
using Quaternions

import MeshIO.save

global OSF_PATH = joinpath(@__DIR__, "../osf-pytorch")
global OSF_CONFIG_FOLDER = joinpath("configs/osf")

function osf_loader_dir()
    joinpath(@__DIR__, "..")
end

include("normalizer.jl")
include("density.jl")
include("load.jl")
include("mesh.jl")
include("collider.jl")
include("scan.jl")

export
    inertia_properties,
    sample_soft,
    finite_difference_gradient

export
    get_nerf_object,
    density_query,
    grid_particles,
    slice_particles,
    grid_density,
    slice_density,
    nerf_density

export
    load_python_script

export
    nerf_mesh

export
    particle_normalization

export
    nerf_scan

end # module
