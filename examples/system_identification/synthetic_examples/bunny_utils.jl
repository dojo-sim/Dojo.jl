################################################################################
# Generate & Save Dataset
################################################################################
function initial_state_bunny(;
	   xlims=[[0,0,0.2], [1,1,0.4]],
       vlims=[[-2,-2,-0.5], [2,2,0.]],
       ωlims=[-6ones(3), 6ones(3)],)
	x = xlims[1] + rand(3) .* (xlims[2] - xlims[1])
	v = vlims[1] + rand(3) .* (vlims[2] - vlims[1])
	ω = ωlims[1] + rand(3) .* (ωlims[2] - ωlims[1])
	return Dict(:position => x, :velocity => v , :orientation => normalize(rand(Quaternion{Float64})), :angular_velocity => ω)
end

function datafilenamebunny(; N::Int=10, friction_coefficient=0.3)
    "bunny_N_$(N)_friction_coefficient_$(friction_coefficient).jld2"
end

function generate_dataset(model::Symbol;
		N::Int=10, H=2.0, timestep=0.01, gravity=-9.81,
		opts=SolverOptions(btol=1e-6, rtol=1e-6),
		init_kwargs=Dict(), # xlims, vlims, θlims, ωlims...
		mech_kwargs=Dict(), # friction_coefficient, radius, side...
		vis::Visualizer=Visualizer(),
		sleep_ratio = 0.0,
		show_contact = true,
		)
    mechanism = get_mechanism(model, timestep=timestep, gravity=gravity; mech_kwargs...);
    trajs = []
    for i = 1:N
		state = initial_state(model; init_kwargs...)
        initialize!(mechanism, model; state...)
        storage = simulate!(mechanism, H, record=true, opts=opts)
        push!(trajs, storage)
        visualize(mechanism, storage, vis=vis, show_contact=show_contact)
		sleep(H*sleep_ratio)
    end

	data = get_data(mechanism)
    params = Dict(:N => N, :H => H, :timestep => timestep, :gravity => gravity, :data => data)
    jldsave(joinpath(@__DIR__, "..", "data", "dataset", datafilename(model; N = N, mech_kwargs...));
        params=params, trajs=trajs)
    return nothing
end


################################################################################
# Load Dataset
################################################################################
function open_dataset(model::Symbol; kwargs...)
    dataset = jldopen(joinpath(@__DIR__, "..", "data", "dataset", datafilename(model; kwargs...)))
    params = dataset["params"]
    trajs = dataset["trajs"]
    JLD2.close(dataset)
    return params, trajs
end


################################################################################
# Optimization Loss: Evaluation & Gradient
################################################################################
"""
    get_contact_gradients!(mechanism, z, θ; opts)

    return gradients for mechanism with respect to contact data
    note: this requires simulating the mechanism for one time step

    mechanism: Mechanism
    z: state
    u: input
	θ: data
    opts: SolverOptions
"""
function get_contact_gradients!(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, z::AbstractVector{T}, θ::AbstractVector{T};
		opts=SolverOptions()) where {T,Nn,Ne,Nb,Ni}
	z_next = contact_step!(mechanism, z, θ)
    jacobian_state, jacobian_contact = get_contact_gradients(mechanism)
    return z_next, jacobian_state, jacobian_contact
end

function contact_step!(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, z::AbstractVector{T}, θ::AbstractVector{T};
		opts=SolverOptions()) where {T,Nn,Ne,Nb,Ni}
	set_data!(mechanism.contacts, θ)
	nu = input_dimension(mechanism)
	step!(mechanism, z, zeros(nu), opts=opts)
end

function get_contact_gradients(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}) where {T,Nn,Ne,Nb,Ni}
	timestep= mechanism.timestep
	nu = input_dimension(mechanism)
	nc = sum(data_dim.(mechanism.contacts))

	for entry in mechanism.data_matrix.nzval # reset matrix
		entry.value .= 0.0
	end
	jacobian_data!(mechanism.data_matrix, mechanism)
	nodes = [mechanism.joints; mechanism.bodies; mechanism.contacts]
	dimrow = length.(nodes)
	dimcol = data_dim.(nodes)
	index_row = [1+sum(dimrow[1:i-1]):sum(dimrow[1:i]) for i in 1:length(dimrow)]
	index_col = [1+sum(dimcol[1:i-1]):sum(dimcol[1:i]) for i in 1:length(dimcol)]

	index_state = [index_col[body.id][[14:16; 8:10; 17:19; 11:13]] for body in mechanism.bodies] # ∂ x2 v15 q2 ϕ15
	# index_control = [index_col[joint.id][1:input_dimension(joint)] for joint in mechanism.joints] # ∂ u
	index_contact = [index_col[contact.id][1:data_dim(contact)] for contact in mechanism.contacts] # ∂ θ

	datamat = full_matrix(mechanism.data_matrix, dimrow, dimcol)
	solmat = full_matrix(mechanism.system)

	# data Jacobian
	data_jacobian = solmat \ datamat #TODO: use pre-factorization

	# Jacobian
	jacobian_state = zeros(12Nb,12Nb)
	# jacobian_control = zeros(12Nb,nu)
	jacobian_contact = zeros(12Nb,nc)
	for (i, body) in enumerate(mechanism.bodies)
		id = body.id
		# Fill in gradients of v25, ϕ25
		jacobian_state[12*(i-1) .+ [4:6; 10:12],:] += data_jacobian[index_row[id], vcat(index_state...)]
		# jacobian_control[12*(i-1) .+ [4:6; 10:12],:] += data_jacobian[index_row[id], vcat(index_control...)]
		jacobian_contact[12*(i-1) .+ [4:6; 10:12],:] += data_jacobian[index_row[id], vcat(index_contact...)]

		# Fill in gradients of x3, q3
		x2 = body.state.x2
		q2 = body.state.q2
		v25 = body.state.vsol[2]
		ϕ25 = body.state.ϕsol[2]
		q3 = next_orientation(q2, ϕ25, timestep)
		jacobian_state[12*(i-1) .+ (1:3), :] += linear_integrator_jacobian_velocity(x2, v25, timestep) * data_jacobian[index_row[id][1:3], vcat(index_state...)]
		jacobian_state[12*(i-1) .+ (1:3), 12*(i-1) .+ (1:3)] += linear_integrator_jacobian_position(x2, v25, timestep)
		jacobian_state[12*(i-1) .+ (7:9), :] += LVᵀmat(q3)' * rotational_integrator_jacobian_velocity(q2, ϕ25, timestep) * data_jacobian[index_row[id][4:6], vcat(index_state...)]
		jacobian_state[12*(i-1) .+ (7:9), 12*(i-1) .+ (7:9)] += LVᵀmat(q3)' * rotational_integrator_jacobian_orientation(q2, ϕ25, timestep, attjac=true)

		# jacobian_control[12*(i-1) .+ (1:3),:] += linear_integrator_jacobian_velocity(x2, v25, timestep) * data_jacobian[index_row[id][1:3], vcat(index_control...)]
		# jacobian_control[12*(i-1) .+ (7:9),:] += LVᵀmat(q3)' * rotational_integrator_jacobian_velocity(q2, ϕ25, timestep) * data_jacobian[index_row[id][4:6], vcat(index_control...)]

		jacobian_contact[12*(i-1) .+ (1:3),:] += linear_integrator_jacobian_velocity(x2, v25, timestep) * data_jacobian[index_row[id][1:3], vcat(index_contact...)]
		jacobian_contact[12*(i-1) .+ (7:9),:] += LVᵀmat(q3)' * rotational_integrator_jacobian_velocity(q2, ϕ25, timestep) * data_jacobian[index_row[id][4:6], vcat(index_contact...)]
	end

	return jacobian_state, jacobian_contact
end


function loss(mechanism::Mechanism, θ::AbstractVector{T}, traj::Storage{T,N},
		indices::UnitRange{Int}; opts=SolverOptions(btol=1e-6, rtol=1e-6), derivatives::Bool=false) where {T,N}

	ni = length(indices)
	nz = maximal_dimension(mechanism, attjac=true)
	nd = length(θ)
	cost = 0.0
	grad = zeros(nd)
	hess = zeros(nd,nd)

	d_contact = zeros(nz,nd)
	for i in indices
		z = get_maximal_state(traj, i)
		z_true = get_maximal_state(traj, i+1)
		Q = Diagonal([ones(3); 1e-1ones(3); ones(4); 1e-1ones(3)])

		if derivatives
			z_pred, ∂_state, ∂_contact = get_contact_gradients!(mechanism, z, θ, opts=opts)
			d_contact = ∂_contact + ∂_state * d_contact
			attjac = attitude_jacobian(z_pred, 1)
			grad += - (attjac * d_contact)' * Q * (z_pred - z_true)
			hess += (attjac * d_contact)' * Q * (attjac * d_contact)
		else
			z_pred = contact_step!(mechanism, z, θ, opts=opts)
		end
		cost += 0.5 * (z_pred - z_true)'* Q *(z_pred - z_true)
	end
	if derivatives
		return cost/ni, grad/ni, hess/ni
	else
		return cost/ni
	end
end

function clean_loss(model::Symbol, trajs, data; timestep=0.01, g=-9.81,
		opts=SolverOptions(btol=1e-6, rtol=1e-6), n_sample=20, rot=0.0)

	mechanism = get_mechanism(model, timestep=timestep, gravity=gravity)
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

# length(1:10)
#
# L1 = [loss(mech, [x; θ0[2:4]], trajs0[1], 100:199) for x in 0:0.05:1]
# L2 = [loss(mech, [θ0[1]; x; θ0[3:4]], trajs0[1], 100:199) for x in -0.1:0.005:0.1]
# L3 = [loss(mech, [θ0[1:2]; x; θ0[4]], trajs0[1], 100:199) for x in -0.1:0.005:0.1]
# L4 = [loss(mech, [θ0[1:3]; x], trajs0[1], 100:199) for x in -0.1:0.005:0.1]
# plot(0:0.05:1, [l[1] for l in L1])
# plot(-0.1:0.005:0.1, [l[1] for l in L2])
# plot(-0.1:0.005:0.1, [l[1] for l in L3])
# plot(-0.1:0.005:0.1, [l[1] for l in L4])
#
# loss(mech, θ0, trajs0[1], 1:11, derivatives=true)
