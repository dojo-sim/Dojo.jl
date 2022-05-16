################################################################################
# Optimization Loss: Evaluation & Gradient
################################################################################
function getSimulatorMaxGradients(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}) where {T,Nn,Ne,Nb,Ni}
	timestep= mechanism.timestep
	nu = input_dimension(mechanism)
	nsd = simdata_dim(mechanism)
	njoints = joint_dimension(mechanism)
	datamat = simdata_jacobian(mechanism)
	solmat = full_matrix(mechanism.system)

	data_jacobian = - solmat \ datamat
	∇data_vϕ = data_jacobian[njoints .+ (1:6Nb),:]
	data_jacobian̄ = zeros(13Nb,nsd)
	for (i, body) in enumerate(mechanism.bodies)
		# Fill in gradients of v25, ϕ25
		data_jacobian̄[13*(i-1) .+ [4:6; 11:13],:] += ∇data_vϕ[6*(i-1) .+ (1:6),:]

		# Fill in gradients of x3, q3
		x2 = body.state.x2
		q2 = body.state.q2
		v25 = body.state.vsol[2]
		ϕ25 = body.state.ϕsol[2]
		data_jacobian̄[13*(i-1) .+ (1:3),:] += linear_integrator_jacobian_velocity(x2, v25, timestep) * ∇data_vϕ[6*(i-1) .+ (1:3),:]
		data_jacobian̄[13*(i-1) .+ (7:10),:] += rotational_integrator_jacobian_velocity(q2, ϕ25, timestep) * ∇data_vϕ[6*(i-1) .+ (4:6),:]
	end
	return data_jacobian̄
end

function getSimulatorMaxGradients!(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, z::AbstractVector{T}, u::AbstractVector{T};
		opts=SolverOptions()) where {T,Nn,Ne,Nb,Ni}
	step!(mechanism, z, u, opts=opts)
	data_jacobian̄ = getSimulatorMaxGradients(mechanism)
	return data_jacobian̄
end

function clean_loss(model::Symbol, pairs, data; timestep=0.05, g=-9.81,
		opts=SolverOptions(btol=1e-6, rtol=1e-6), n_sample = 20, rot = 0.0)
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

function loss(mechanism::Mechanism, pair; opts=SolverOptions(btol=1e-6, rtol=1e-6))
	nu = input_dimension(mechanism)
	nz = maximal_dimension(mechanism)
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
