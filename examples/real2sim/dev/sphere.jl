# Utils
function module_dir()
    return joinpath(@__DIR__, "..", "..", "..")
end

# Activate package
using Pkg
Pkg.activate(module_dir())

# Load packages
using Plots
using Random
using MeshCat

# Open visualizer
vis = Visualizer()
open(vis)

# Include new files
include(joinpath(module_dir(), "examples", "loader.jl"))
include(joinpath(module_dir(), "env", "sphere", "deps", "texture.jl"))

mech = getmechanism(:sphere, Δt = 0.01, g = -9.81, radius = 0.5, cf = 0.1);
initialize!(mech, :sphere, x=[0,0,1.], v=[0,0.5,0.], ω=[10,0,0.])
storage = simulate!(mech, 4.0, record=true, verbose=true,
    opts=InteriorPointOptions(btol=1e-6, rtol=1e-6))
visualize(mech, storage, vis=vis)
sphere_texture!(vis, mech)


mech = getmechanism(:box, Δt = 0.01, g = -9.81, cf = 0.1);
initialize!(mech, :box, x=[0,0,0.5], v=[1,1.5,1.], ω=[5,4,2.])
storage = simulate!(mech, 0.64, record=true, verbose=true,
    opts=InteriorPointOptions(btol=1e-6, rtol=1e-6, verbose=true, max_ls=3))
visualize(mech, storage, vis=vis)
sphere_texture!(vis, mech)

const Dojo = Main

################################################################################
# Generate & Save Dataset
################################################################################
function filename(; N::Int = 10, cf = 0.1, radius = 0.5)
    "N:$(N)_cf:$(cf)_radius:$(radius),jld2"
end

function build_triplets(trajs::AbstractVector)
    triplets = []
    for traj in trajs
        push!(triplets, build_triplets(traj)...)
    end
    return triplets
end

function build_triplets(traj::Storage{T,N}) where {T,N}
    triplets = []
    for t = 1:N-2
        t1 = [traj.x[1][t]; vector(traj.q[1][t])]
        t2 = [traj.x[1][t+1]; vector(traj.q[1][t+1])]
        t3 = [traj.x[1][t+2]; vector(traj.q[1][t+2])]
        triplet = [t1, t2, t3]
        push!(triplets, triplet)
    end
    return triplets
end

function generate_dataset(; N::Int = 10, H = 2.0, Δt = 0.01, g = -9.81, cf = 0.1, radius = 0.5,
        xlims=[zeros(3), ones(3)],
        vlims=[-1*ones(3), ones(3)],
        ωlims=[-1*ones(3), ones(3)])

    mechanism = getmechanism(:sphere, Δt = Δt, g = g, cf = cf, radius = radius)
    trajs = []
    for i = 1:N
        x = xlims[1] + rand(3) .* (xlims[2] - xlims[1])
        v = vlims[1] + rand(3) .* (vlims[2] - vlims[1])
        ω = ωlims[1] + rand(3) .* (ωlims[2] - ωlims[1])
        initialize!(mech, :sphere, x=x, v=v, ω=ω)
        storage = simulate!(mech, H, record=true, opts=InteriorPointOptions(btol=1e-6, rtol=1e-6))
        push!(trajs, storage)
        visualize(mech, storage, vis=vis)
        sphere_texture!(vis, mech)
    end
    pars = Dict(:N => N, :H => H, :Δt => Δt, :g => g, :cf => cf, :radius => radius)
    triplets = build_triplets(trajs)
    jldsave(joinpath(@__DIR__, "dataset", filename(N=N, cf=cf, radius=radius));
        pars=pars, trajs=trajs, triplets=triplets)
    return nothing
end

################################################################################
# Load Dataset
################################################################################
function open_dataset(; N::Int=10, cf = 0.1, radius = 0.5)
    dataset = jldopen(joinpath(@__DIR__, "dataset", filename(N=N, cf=cf, radius=radius)))
    pars = dataset["pars"]
    trajs = dataset["trajs"]
    triplets = dataset["triplets"]
    JLD2.close(dataset)
    return pars, trajs, triplets
end


################################################################################
# Optimization Loss: Evaluation & Gradient
################################################################################
function loss(triplets, cf, radius; N::Int = 10, H = 2.0, Δt = 0.01, g = -9.81)
    mechanism = getmechanism(:sphere, Δt = Δt, g = g, cf = cf, radius = radius)
    l = 0.0
    ∇ = zeros(2)
    for triplet in triplets
        li, ∇i += loss(mechanism, triplet)
        l += li
        ∇ += ∇i
    end
    return l, ∇
end

function loss(mechanism::Mechanism, triplet)
    u = zeros(controldim(mechanism))
    z = triplet[1]
    ztrue = triplet[2]
    z = step(mechanism, z, u)
    l = norm(z - ztrue)
    ∇ = FiniteDiff.finite_difference_gradient( -  )
    return l, ∇
end





################################################################################
# Optimization Algorithm: L-BFGS
################################################################################


################################################################################
# Visualization
################################################################################





################################################################################
# Generate & Save Dataset
################################################################################
generate_dataset()

################################################################################
# Load Dataset
################################################################################
pars, trajs, triplets = open_dataset()

################################################################################
# Optimization Objective: Evaluation & Gradient
################################################################################


################################################################################
# Optimization Algorithm: L-BFGS
################################################################################


################################################################################
# Visualization
################################################################################
