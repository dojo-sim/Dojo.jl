################################################################################
# Simulator Data
################################################################################

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
function filename(model::Symbol; kwargs...)
    eval(Symbol(:filename, model))(; kwargs...)
end

function filenamesphere(; N::Int=10, cf=0.1, radius=0.5, p=szeros(3))
    "sphere_dim_N:$(N)_cf:$(cf)_radius:$(radius)_p:$(p).jld2"
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
    z = getMaxState(traj)
    for t = 1:N-1
        z1 = z[t]
        z2 = z[t+1]
        pair = [z1, z2]
        push!(pairs, pair)
    end
    return pairs
end

function generate_dataset(model::Symbol; N::Int=10, H=2.0, Δt=0.05, g=-9.81,
		cf=0.1, radius=0.5, p=szeros(3),
        xlims=[zeros(3), ones(3)],
        vlims=[-1*ones(3), ones(3)],
        ωlims=[-1*ones(3), ones(3)],
        opts=InteriorPointOptions(btol=1e-6, rtol=1e-6))

    mechanism = getmechanism(model, Δt=Δt, g=g, cf=cf, radius=radius)
    trajs = []
    for i = 1:N
        x = xlims[1] + rand(3) .* (xlims[2] - xlims[1])
        v = vlims[1] + rand(3) .* (vlims[2] - vlims[1])
        ω = ωlims[1] + rand(3) .* (ωlims[2] - ωlims[1])
        initialize!(mechanism, model, x=x, v=v, ω=ω)
        storage = simulate!(mechanism, H, record=true, opts=opts)
        push!(trajs, storage)
        visualize(mechanism, storage, vis=vis)
        sphere_texture!(vis, mechanism)
    end
	data = get_simulator_data(mechanism)
    params = Dict(:N => N, :H => H, :Δt => Δt, :g => g, :data => data)
    pairs = build_pairs(mechanism, trajs)
    jldsave(joinpath(@__DIR__, "dev", "dataset", filename(model, N=N, cf=cf, radius=radius));
        params=params, trajs=trajs, pairs=pairs)
    return nothing
end


################################################################################
# Load Dataset
################################################################################
function open_dataset(model::Symbol ; N::Int=10, cf=0.1, radius=0.5, p=szeros(3))
    dataset = jldopen(joinpath(@__DIR__, "dev", "dataset", filename(model, N=N, cf=cf, radius=radius)))
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

function loss(mechanism, pairs, data; Δt=0.05, g=-9.81,
		opts=InteriorPointOptions(btol=1e-6, rtol=1e-6), n_sample = 20)
	nsd = simulator_data_dim(mechanism)
	mechanism = getmechanism(:sphere, Δt=Δt, g=g)
	set_simulator_data!(mechanism, data)

	l = 0.0
    ∇ = zeros(nsd)
	n = length(pairs)

	mask = rand(1:n, n_sample)
    for i in mask
        li, ∇i = loss(mechanism, pairs[i], opts=opts)
        l += li
        ∇ += ∇i
    end
    return l / n_sample, ∇ ./ n_sample
end

function loss(mechanism::Mechanism, pair; opts=InteriorPointOptions(btol=1e-6, rtol=1e-6))
	nu = controldim(mechanism)
	nz = maxCoordDim(mechanism)
	u = zeros(nu)
    z1 = pair[1]
    z2true = pair[2]
    z2 = step!(mechanism, z1, u, opts=opts)
    ∇z2 = getSimulatorMaxGradients!(mechanism, z1, u, opts=opts)
    l = 0.5 * (z2 - z2true)'*(z2 - z2true)
    ∇ = - ∇z2' * (z2 - z2true)
    return l, ∇
end
