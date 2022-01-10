################################################################################
# Simulator Data
################################################################################

simdata_dim(mechanism::Mechanism) = sum(simdata_dim.(mech.ineqconstraints))
simdata_dim(ineqc::InequalityConstraint) = sum(simdata_dim.(ineqc.constraints))
simdata_dim(bound::ContactBound) = 7 # [cf, p, offset]
simdata_dim(bound::LinearContactBound) = 7 # [cf, p, offset]
simdata_dim(bound::ImpactBound) = 6 # [p, offset]

function set_simulator_data!(bound::ContactBound, data::AbstractVector)
	bound.cf = data[1]
    bound.offset = data[SVector{3,Int}(2:4)]
    bound.p = data[SVector{3,Int}(5:7)]
    return nothing
end

function set_simulator_data!(bound::LinearContactBound, data::AbstractVector)
	bound.cf = data[1]
    bound.offset = data[SVector{3,Int}(2:4)]
    bound.p = data[SVector{3,Int}(5:7)]
    return nothing
end

function set_simulator_data!(bound::ImpactBound, data::AbstractVector)
    bound.offset = data[SVector{3,Int}(2:4)]
    bound.p = data[SVector{3,Int}(5:7)]
    return nothing
end

function set_simulator_data!(mechanism::Mechanism, data::AbstractVector)
    c = 0
    for ineqc in mechanism.ineqconstraints
		for bound in ineqc.constraints
			N = simdata_dim(bound)
	        set_simulator_data!(bound, data[c .+ (1:N)]); c += N
		end
    end
    return nothing
end

function get_simulator_data(bound::ContactBound)
	return [bound.cf; bound.offset; bound.p]
end

function get_simulator_data(boundset::LinearContactBound)
	return [bound.cf; bound.off; bound.p]
end

function get_simulator_data(bound::ImpactBound)
	return [bound.offset; bound.p]
end

function get_simulator_data(mechanism::Mechanism{T}) where T
    data = zeros(T,simdata_dim(mechanism))
	e = 0
    for ineqc in mechanism.ineqconstraints
		for bound in ineqc.constraints
			N = simdata_dim(bound)
        	data[e .+ (1:N)] = get_simulator_data(bound); e += N
		end
    end
    return data
end

################################################################################
# Simulator Data Jacobian
################################################################################

function simdata_jacobian(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}) where {T,Nn,Ne,Nb,Ni}
	nsol = soldim(mechanism)
    ndata = simdata_dim(mechanism)
    data = get_simulator_data(mechanism)

	A = zeros(nsol, ndata)
	neqcs = eqcdim(mechanism)
	n = sum(length.(mech.eqconstraints.values))
	rb = neqcs # row ineqc
	ri = neqcs + 6Nb # row body
	c = 0 # column
	for ineqc in mechanism.ineqconstraints
		for bound in ineqc.constraints
			N = length(ineqc)
			N½ = Int(N/2)
			Ns = simdata_dim(bound)
			∇ineqc, ∇body = ∂g∂simdata(mechanism, ineqc) # assumes only one bound per ineqc
			A[rb + 3 .+ (1:3), c .+ (1:Ns)] -= ∇body # only works for one body for now
			A[ri + N½ .+ (1:N½), c .+ (1:Ns)] -= ∇ineqc
			ri += N
			c += Ns
		end
	end
	return A
end

function ∂g∂simdata(mechanism, ineqc::InequalityConstraint{T,N,Nc,Cs}) where {T,N,Nc,Cs<:Tuple{ContactBound{T,N}}}
    bound = ineqc.constraints[1]
	p = bound.p
	offset = bound.offset
    body = getbody(mechanism, ineqc.parentid)
    x2, v25, q2, ϕ25 = fullargssol(body.state)
    x3, q3 = posargs3(body.state, mechanism.Δt)
	s = ineqc.ssol[2]
	γ = ineqc.γsol[2]

	# Contribution to Ineqconstraint
	∇cf = SA[0,γ[1],0,0]
	∇off = [-bound.ainv3; szeros(T,1,3); -bound.Bx * skew(vrotate(ϕ25, q3))]
	∇p = [bound.ainv3 * ∂vrotate∂p(bound.p, q3); szeros(T,1,3); bound.Bx * skew(vrotate(ϕ25, q3)) * ∂vrotate∂p(bound.p, q3)]
	∇ineqc = [∇cf ∇off ∇p]

	# Contribution to Body dynamics
	∇cf = szeros(T,3)
	X = forcemapping(bound)
	# this what we differentiate: Qᵀγ = - skew(p - vrotate(offset, inv(q3))) * VRmat(q3) * LᵀVᵀmat(q3) * X' * γ
	∇off = - ∂pskew(VRmat(q3) * LᵀVᵀmat(q3) * X' * γ) * -∂vrotate∂p(offset, inv(q3))
	∇p = - ∂pskew(VRmat(q3) * LᵀVᵀmat(q3) * X' * γ)
	∇body = [∇cf ∇off ∇p]
	return ∇ineqc, ∇body
end


################################################################################
# Generate & Save Dataset
################################################################################
function filename(model::Symbol; kwargs...)
    eval(Symbol(:filename, model))(; kwargs...)
end

function filenamesphere(; N::Int=10, cf=0.1, radius=0.5)
    "sphere_dim_N:$(N)_cf:$(cf)_radius:$(radius).jld2"
end

function filenamebox2d(; N::Int=10, cf=0.1, radius=0.05, side=0.50)
    "box2d_dim_N:$(N)_cf:$(cf)_radius:$(radius)_side:$(side).jld2"
end

function filenamebox(; N::Int=10, cf=0.1, radius=0., side=0.50)
    "box_dim_N:$(N)_cf:$(cf)_radius:$(radius)_side:$(side).jld2"
end

function filenamehardwarebox(; N::Int=10, S::Int=1)
    "hardwarebox_dim_N:$(N)_sample_S:$(S).jld2"
end

function initial_state(model::Symbol; kwargs...)
    eval(Symbol(:initial_state_, model))(; kwargs...)
end

function initial_state_sphere(;
	xlims=[[0,0,0], [1,1,0.2]],
	vlims=[-1ones(3), 1ones(3)],
	ωlims=[-5ones(3), 5ones(3)],)
	x = xlims[1] + rand(3) .* (xlims[2] - xlims[1])
	v = vlims[1] + rand(3) .* (vlims[2] - vlims[1])
	ω = ωlims[1] + rand(3) .* (ωlims[2] - ωlims[1])
	return Dict(:x => x, :v => v , :ω => ω)
end

function initial_state_box2d(;
		xlims=[[0,0.2], [1,0.4]],
        vlims=[-ones(2), ones(2)],
		θlims=[-π, π],
        ωlims=[-10, 10],)
	x = xlims[1] + rand(2) .* (xlims[2] - xlims[1])
	v = vlims[1] + rand(2) .* (vlims[2] - vlims[1])
	θ = θlims[1] + rand() * (θlims[2] - θlims[1])
	ω = ωlims[1] + rand() * (ωlims[2] - ωlims[1])
	return Dict(:x => x, :v => v , :θ => θ, :ω => ω)
end

function initial_state_box(;
		xlims=[[0,0,0.2], [1,1,0.4]],
        vlims=[-2ones(3), [2,2,-1.]],
        ωlims=[-6ones(3), 6ones(3)],)
	x = xlims[1] + rand(3) .* (xlims[2] - xlims[1])
	v = vlims[1] + rand(3) .* (vlims[2] - vlims[1])
	ω = ωlims[1] + rand(3) .* (ωlims[2] - ωlims[1])
	return Dict(:x => x, :v => v , :q => rand(UnitQuaternion), :ω => ω)
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

function generate_dataset(model::Symbol;
		N::Int=10, H=2.0, Δt=0.05, g=-9.81,
		opts=InteriorPointOptions(btol=1e-6, rtol=1e-6),
		init_kwargs=Dict(), # xlims, vlims, θlims, ωlims...
		mech_kwargs=Dict(), # cf, radius, side...
		sleep_ratio = 0.0,
		show_contact = true,
		)
    mechanism = getmechanism(model, Δt=Δt, g=g; mech_kwargs...);
    trajs = []
    for i = 1:N
		state = initial_state(model; init_kwargs...)
        initialize!(mechanism, model; state...)
        storage = simulate!(mechanism, H, record=true, opts=opts)
        push!(trajs, storage)
        visualize(mechanism, storage, vis=vis, show_contact=show_contact)
		sleep(H*sleep_ratio)
    end

	data = get_simulator_data(mechanism)
    params = Dict(:N => N, :H => H, :Δt => Δt, :g => g, :data => data)
    pairs = build_pairs(mechanism, trajs)
    jldsave(joinpath(@__DIR__, "data", "dataset", filename(model; N = N, mech_kwargs...));
        params=params, trajs=trajs, pairs=pairs)
    return nothing
end


################################################################################
# Load Dataset
################################################################################
function open_dataset(model::Symbol; kwargs...)
    dataset = jldopen(joinpath(@__DIR__, "data", "dataset", filename(model; kwargs...)))
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
	nsd = simdata_dim(mechanism)
	neqcs = eqcdim(mechanism)
	datamat = simdata_jacobian(mechanism)
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

function clean_loss(model::Symbol, pairs, data; Δt=0.05, g=-9.81,
		opts=InteriorPointOptions(btol=1e-6, rtol=1e-6), n_sample = 20, rot = 0.0)
	mechanism = getmechanism(model, Δt=Δt, g=g)
	nsd = simdata_dim(mechanism)
	set_simulator_data!(mechanism, data)

	l = 0.0
	∇ = zeros(nsd)
    H = zeros(nsd, nsd)
	n = length(pairs)

	offset = Int(floor(rot))
	Random.seed!(0)
	mask = shuffle!(Vector(1:n))
	mask = [mask; mask]
	mask = mask[offset%n .+ (1:n_sample)]
    for i in mask
        li, ∇i, Hi = loss(mechanism, pairs[i], opts=opts)
        l += li
        ∇ += ∇i
		H += Hi
    end
    return l / n_sample, ∇ ./ n_sample, H ./ n_sample
end

function loss(mechanism::Mechanism, pair; opts=InteriorPointOptions(btol=1e-6, rtol=1e-6))
	nu = controldim(mechanism)
	nz = maxCoordDim(mechanism)
	u = zeros(nu)
    z1 = pair[1]
    z2true = pair[2]
    z2 = step!(mechanism, z1, u, opts=opts)
    ∇z2 = getSimulatorMaxGradients(mechanism)
	# @show scn.(z2 - z2true)
	Q = Diagonal([ones(3); 1e-1ones(3); ones(4); 1e-1ones(3)])
    l = 0.5 * (z2 - z2true)'* Q *(z2 - z2true)
    ∇ = - ∇z2' * Q * (z2 - z2true)
	H = ∇z2'* Q * ∇z2 # quasi newton approx oft he Hessian
    return l, ∇, H
end
