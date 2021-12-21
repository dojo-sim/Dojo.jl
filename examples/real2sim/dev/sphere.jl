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

mech = getmechanism(:sphere, Δt=0.01, g=-9.81, radius=0.5, cf=0.1);
initialize!(mech, :sphere, x=[0,0,1.], v=[0,0.5,0.], ω=[10,0,0.])
storage = simulate!(mech, 4.0, record=true, verbose=true,
    opts=InteriorPointOptions(btol=1e-6, rtol=1e-6))
visualize(mech, storage, vis=vis)
sphere_texture!(vis, mech)


# mech = getmechanism(:box, Δt=0.01, g=-9.81, cf=0.1);
# initialize!(mech, :box, x=[0,0,0.5], v=[1,1.5,1.], ω=[5,4,2.])
# storage = simulate!(mech, 5.0, record=true,
#     opts=InteriorPointOptions(btol=1e-6, rtol=1e-6, verbose=false))
# visualize(mech, storage, vis=vis)
#
# collect(mech.ineqconstraints)[1].constraints[1].cf
# collect(mech.ineqconstraints)[1].constraints[1].offset
#  # = sones(3)
# collect(mech.ineqconstraints)[1].constraints[1].p
#  # = sones(3)


simulator_data_dim(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}) where {T,Nn,Ne,Nb,Ni} = 5Ni

function unpack_simulator_data(data::AbstractVector)
    cf = data[SVector{1,Int}(1:1)]
    off = data[SVector{1,Int}(2:2)]
    p = data[SVector{3,Int}(3:5)]
    return cf, off, p
end

function set_simulator_data!(mechanism::Mechanism, data::AbstractVector)
    c = 0
    for ineqc in mechanism.ineqconstraints
        cf, off, p = unpack_simulator_data(data[c .+ (1:5)]); c += 5
        ineqc.constraints[1].cf = cf[1]
        ineqc.constraints[1].offset = [szeros(2); off]
        ineqc.constraints[1].p = p
    end
    return nothing
end

function get_simulator_data(mechanism::Mechanism{T}) where T
    data = Vector{T}()
    for ineqc in mechanism.ineqconstraints
        cf = ineqc.constraints[1].cf
        off = ineqc.constraints[1].offset[3]
        p = ineqc.constraints[1].p
        push!(data, cf, off, p...)
    end
    return data
end

function simulator_data_jacobian(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}) where {T,Nn,Ne,Nb,Ni}
    mechanism = deepcopy(mechanism)
    nsol = soldim(mechanism)
    ndata = simulator_data_dim(mechanism)
    data = get_simulator_data(mechanism)

    function res(data)
        set_simulator_data!(mechanism, data)
        setentries!(mechanism)
        return Vector(full_vector(mechanism.system))
    end
    A = FiniteDiff.finite_difference_jacobian(data -> res(data), data)
    return A
end



################################################################################
# Generate & Save Dataset
################################################################################
function filename(; N::Int=10, cf=0.1, radius=0.5, p=szeros(3))
    "N:$(N)_cf:$(cf)_radius:$(radius)_p:$(p),jld2"
end

function build_pairs(mechanism::Mechanism, trajs::AbstractVector)
    pairs = []
    for traj in trajs
        push!(pairs, build_pairs(mechanism, traj)...)
    end
    return pairs
end

function build_pairs(mechanism::Mechanism, traj::Storage{T,N}) where {T,N}
    pairs = []
    z = getMaxState(storage)
    for t = 1:N-1
        z1 = z[t]
        z2 = z[t+1]
        pair = [z1, z2]
        push!(pairs, pair)
    end
    return pairs
end

function generate_dataset(; N::Int=10, H=2.0, Δt=0.01, g=-9.81,
		cf=0.1, radius=0.5, p=szeros(3),
        xlims=[zeros(3), ones(3)],
        vlims=[-1*ones(3), ones(3)],
        ωlims=[-1*ones(3), ones(3)],
        opts=InteriorPointOptions(btol=1e-6, rtol=1e-6))

    mechanism = getmechanism(:sphere, Δt=Δt, g=g, cf=cf, radius=radius)
    trajs = []
    for i = 1:N
        x = xlims[1] + rand(3) .* (xlims[2] - xlims[1])
        v = vlims[1] + rand(3) .* (vlims[2] - vlims[1])
        ω = ωlims[1] + rand(3) .* (ωlims[2] - ωlims[1])
        initialize!(mech, :sphere, x=x, v=v, ω=ω)
        storage = simulate!(mech, H, record=true, opts=opts)
        push!(trajs, storage)
        visualize(mech, storage, vis=vis)
        sphere_texture!(vis, mech)
    end
	data = get_simulator_data(mechanism)
    params = Dict(:N => N, :H => H, :Δt => Δt, :g => g, :data => data)
    pairs = build_pairs(mechanism, trajs)
    jldsave(joinpath(@__DIR__, "dataset", filename(N=N, cf=cf, radius=radius));
        params=params, trajs=trajs, pairs=pairs)
    return nothing
end

################################################################################
# Load Dataset
################################################################################
function open_dataset(; N::Int=10, cf=0.1, radius=0.5, p=szeros(3))
    dataset = jldopen(joinpath(@__DIR__, "dataset", filename(N=N, cf=cf, radius=radius)))
    params = dataset["params"]
    trajs = dataset["trajs"]
    pairs = dataset["pairs"]
    JLD2.close(dataset)
    return params, trajs, pairs
end


################################################################################
# Optimization Loss: Evaluation & Gradient
################################################################################
function getSimulatorMaxGradients(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}) where {T,Nn,Ne,Nb,Ni}
	Δt = mechanism.Δt
	nu = controldim(mechanism)
	nsd = simulator_data_dim(mechanism)
	neqcs = eqcdim(mechanism)
	datamat = simulator_data_jacobian(mechanism)
	solmat = full_matrix(mechanism.system)

	∇data_z = - solmat \ datamat
	∇data_vϕ = ∇data_z[neqcs .+ (1:6Nb),:]
	∇data_z̄ = zeros(13Nb,nsd)
	for (i, body) in enumerate(mechanism.bodies)
		# Fill in gradients of v25, ϕ25
		∇data_z̄[13*(i-1) .+ [4:6; 11:13],:] += ∇data_vϕ[6*(i-1) .+ (1:6),:]

		# Fill in gradients of x3, q3
		q2 = body.state.q2[1]
		ϕ25 = body.state.ϕsol[2]
		∇data_z̄[13*(i-1) .+ (1:3),:] += ∂integrator∂v(Δt) * ∇data_vϕ[6*(i-1) .+ (1:3),:]
		∇data_z̄[13*(i-1) .+ (7:10),:] += ∂integrator∂ϕ(q2, ϕ25, Δt) * ∇data_vϕ[6*(i-1) .+ (4:6),:]
	end
	return ∇data_z̄
end

function getSimulatorMaxGradients!(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, z::AbstractVector{T}, u::AbstractVector{T};
		opts=InteriorPointOptions()) where {T,Nn,Ne,Nb,Ni}
	step!(mechanism, z, u, opts=opts)
	∇data_z̄ = getSimulatorMaxGradients(mechanism)
	return ∇data_z̄
end

function loss(mechanism, pairs, data; Δt=0.01, g=-9.81,
		opts=InteriorPointOptions(btol=1e-6, rtol=1e-6))
	nsd = simulator_data_dim(mechanism)
	mechanism = getmechanism(:sphere, Δt=Δt, g=g)
	set_simulator_data!(mechanism, data)

	l = 0.0
    ∇ = zeros(nsd)


    for pair in pairs
        li, ∇i = loss(mechanism, pair, opts=opts)
        l += li
        ∇ += ∇i
    end
	n = length(pairs)
    return l / n, ∇ ./ n
end

function loss(mechanism::Mechanism, pair; opts=InteriorPointOptions(btol=1e-6, rtol=1e-6))
    u = zeros(controldim(mechanism))
    z1 = pair[1]
    z2true = pair[2]
    z2 = step!(mechanism, z1, u, opts=opts)
    ∇z2 = getSimulatorMaxGradients!(mechanism, z1, u, opts=opts)
    l = 0.5 * (z2 - z2true)'*(z2 - z2true)
    ∇ = ∇z2' * (z2 - z2true)
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
generate_dataset(H = 0.5, N = 50,
	xlims = [[0,0,0], [1,1,0.2]],
	ωlims = [-5ones(3), 5ones(3)],
	opts=InteriorPointOptions(btol=3e-4, rtol=3e-4))

################################################################################
# Load Dataset
################################################################################
params0, trajs0, pairs0 = open_dataset()


data0 = params0[:data]
data1 = [0.1, 0.49, 0,0,0]
data2 = [0.1, 0.51, 0,0,0]
loss(mech, pairs0, data_0, opts=InteriorPointOptions(btol=3e-4, rtol=3e-4))
loss(mech, pairs0, data0, opts=InteriorPointOptions(btol=3e-4, rtol=3e-4))
loss(mech, pairs0, data1, opts=InteriorPointOptions(btol=3e-4, rtol=3e-4))
loss(mech, pairs0, data2, opts=InteriorPointOptions(btol=3e-4, rtol=3e-4))

# data_0 = [0.1, 1.0, 0,0,0]
for k = 1:10
	l, ∇ = loss(mech, pairs0, data_0, opts=InteriorPointOptions(btol=3e-4, rtol=3e-4))
	data_0[1:2] += 0.1 * ∇[1:2] ./ max(norm(∇[1:2], Inf), 1e1)
	@show scn.(∇[1:2])
	@show scn.(data_0[1:2])
end
data_0


################################################################################
# Optimization Objective: Evaluation & Gradient
################################################################################


################################################################################
# Optimization Algorithm: L-BFGS
################################################################################


################################################################################
# Visualization
################################################################################
