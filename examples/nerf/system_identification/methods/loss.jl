################################################################################
# Optimization Loss: Evaluation & Gradient
################################################################################
"""
    get_body_contact_gradients!(mechanism, z, θ; opts)

    return gradients for mechanism with respect to body & contact data
    note: this requires simulating the mechanism for one time step

    mechanism: Mechanism
    z: state
    u: input
	θ: data for the whole mechanism
    opts: SolverOptions
"""
function get_body_contact_gradients!(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, z::AbstractVector{T},
		θ::AbstractVector{T}; opts=SolverOptions()) where {T,Nn,Ne,Nb,Ni}
	z_next = body_contact_step!(mechanism, z, θ)
    jacobian_state, jacobian_body, jacobian_contact = get_body_contact_gradients(mechanism)
    return z_next, jacobian_state, jacobian_body, jacobian_contact
end

function body_contact_step!(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, z::AbstractVector{T},
		θ::AbstractVector{T}; opts=SolverOptions()) where {T,Nn,Ne,Nb,Ni}
	set_data!(mechanism, θ)
	nu = input_dimension(mechanism)
	step!(mechanism, z, zeros(nu), opts=opts)
end

function get_body_contact_gradients(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}) where {T,Nn,Ne,Nb,Ni}
	timestep= mechanism.timestep
	nu = input_dimension(mechanism)
	nb = 7 * length(mechanism.bodies)
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

	index_state = [index_col[body.id][[14:16; 8:10; 17:19; 11:13]] for body in mechanism.bodies] # ∂ x2 v15 q2 ω15
	# index_control = [index_col[joint.id][1:input_dimension(joint)] for joint in mechanism.joints] # ∂ u
	index_body = [index_col[body.id][1:7] for body in mechanism.bodies] # ∂ θbody
	index_contact = [index_col[contact.id][1:data_dim(contact)] for contact in mechanism.contacts] # ∂ θcontact

	datamat = full_matrix(mechanism.data_matrix, dimrow, dimcol)
	solmat = full_matrix(mechanism.system)

	# data Jacobian
	data_jacobian = solmat \ datamat #TODO: use pre-factorization

	# Jacobian
	jacobian_state = zeros(12Nb,12Nb)
	# jacobian_control = zeros(12Nb,nu)
	jacobian_body = zeros(12Nb,nb)
	jacobian_contact = zeros(12Nb,nc)
	for (i, body) in enumerate(mechanism.bodies)
		id = body.id
		# Fill in gradients of v25, ω25
		jacobian_state[12*(i-1) .+ [4:6; 10:12],:] += data_jacobian[index_row[id], vcat(index_state...)]
		# jacobian_control[12*(i-1) .+ [4:6; 10:12],:] += data_jacobian[index_row[id], vcat(index_control...)]
		jacobian_body[12*(i-1) .+ [4:6; 10:12],:] += data_jacobian[index_row[id], vcat(index_body...)]
		jacobian_contact[12*(i-1) .+ [4:6; 10:12],:] += data_jacobian[index_row[id], vcat(index_contact...)]

		# Fill in gradients of x3, q3
		x2 = body.state.x2
		q2 = body.state.q2
		v25 = body.state.vsol[2]
		ω25 = body.state.ωsol[2]
		q3 = next_orientation(q2, ω25, timestep)
		jacobian_state[12*(i-1) .+ (1:3), :] += linear_integrator_jacobian_velocity(x2, v25, timestep) * data_jacobian[index_row[id][1:3], vcat(index_state...)]
		jacobian_state[12*(i-1) .+ (1:3), 12*(i-1) .+ (1:3)] += linear_integrator_jacobian_position(x2, v25, timestep)
		jacobian_state[12*(i-1) .+ (7:9), :] += LVᵀmat(q3)' * rotational_integrator_jacobian_velocity(q2, ω25, timestep) * data_jacobian[index_row[id][4:6], vcat(index_state...)]
		jacobian_state[12*(i-1) .+ (7:9), 12*(i-1) .+ (7:9)] += LVᵀmat(q3)' * rotational_integrator_jacobian_orientation(q2, ω25, timestep, attjac=true)

		# jacobian_control[12*(i-1) .+ (1:3),:] += linear_integrator_jacobian_velocity(x2, v25, timestep) * data_jacobian[index_row[id][1:3], vcat(index_control...)]
		# jacobian_control[12*(i-1) .+ (7:9),:] += LVᵀmat(q3)' * rotational_integrator_jacobian_velocity(q2, ω25, timestep) * data_jacobian[index_row[id][4:6], vcat(index_control...)]

		jacobian_body[12*(i-1) .+ (1:3),:] += linear_integrator_jacobian_velocity(x2, v25, timestep) * data_jacobian[index_row[id][1:3], vcat(index_body...)]
		jacobian_body[12*(i-1) .+ (7:9),:] += LVᵀmat(q3)' * rotational_integrator_jacobian_velocity(q2, ω25, timestep) * data_jacobian[index_row[id][4:6], vcat(index_body...)]

		jacobian_contact[12*(i-1) .+ (1:3),:] += linear_integrator_jacobian_velocity(x2, v25, timestep) * data_jacobian[index_row[id][1:3], vcat(index_contact...)]
		jacobian_contact[12*(i-1) .+ (7:9),:] += LVᵀmat(q3)' * rotational_integrator_jacobian_velocity(q2, ω25, timestep) * data_jacobian[index_row[id][4:6], vcat(index_contact...)]
	end

	return jacobian_state, jacobian_body, -jacobian_contact
end

function loss(mechanism::Mechanism, θ::AbstractVector{T}, traj::Storage{T,N},
		indices::UnitRange{Int}; opts=SolverOptions(btol=1e-6, rtol=1e-6), derivatives::Bool=false) where {T,N}

	ni = length(indices)
	nz = maximal_dimension(mechanism, attjac=true)
	Nb = length(mechanism.bodies)
	nb = 7 * Nb
	nc = sum(data_dim.(mechanism.contacts))
	cost = 0.0
	grad = zeros(nb+nc)
	hess = zeros(nb+nc,nb+nc)

	d_body_contact = zeros(nz,nb+nc)
	z_prev = get_maximal_state(traj, indices[1])
	Z = [deepcopy(z_prev)]
	for i in indices
		z = z_prev
		z_true = get_maximal_state(traj, i+1)
		Q = Diagonal(vcat(fill([ones(3); 1e-1ones(3); ones(4); 1e-1ones(3)], Nb)...))

		if derivatives
			z_pred, ∂_state, ∂_body, ∂_contact = get_body_contact_gradients!(mechanism, z, θ, opts=opts)
			d_body_contact = [∂_body ∂_contact] + ∂_state * d_body_contact
			attjac = attitude_jacobian(z_pred, Nb)
			grad += (attjac * d_body_contact)' * Q * (z_pred - z_true)
			hess += (attjac * d_body_contact)' * Q * (attjac * d_body_contact)
		else
			z_pred = body_contact_step!(mechanism, z, θ, opts=opts)
		end
		cost += 0.5 * (z_pred - z_true)'* Q *(z_pred - z_true)
		z_prev = z_pred
		push!(Z, deepcopy(z_prev))
	end
	if derivatives
		return cost/ni, grad/ni, hess/ni
	else
		return cost/ni, Z
	end
end
