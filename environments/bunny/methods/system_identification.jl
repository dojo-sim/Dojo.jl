using Dojo
using Plots
using BenchmarkTools


vis = Visualizer()
open(vis)


mech = get_bunny(timestep=0.01)
mech.contacts[1].model.collision.collider.options = ColliderOptions()

initialize!(mech, :bunny,
    position=[0,0,0.6],
    # orientation=Quaternion(0,1,0,0.0,true),
    orientation=Quaternion(1,0,0,0.0,true),
    velocity=[0,0.5,5.0],
    angular_velocity=[0.5,10.0,3.0])

constraint(mech, mech.contacts[1])

@elapsed storage = simulate!(mech, 5.0,
    opts=SolverOptions(verbose=false, rtol=1e-4))
visualize(mech, storage, vis=vis)
mech
jacobian_data!(mech.data_matrix, mech)


function body_constraint_jacobian_contact_data(mechanism::Mechanism, body::Body{T},
        contact::SoftContactConstraint{T,N,Nc,Cs}) where {T,N,Nc,Cs}
    Nd = data_dim(contact)
	∇cf = szeros(T,6,1) # sliding friction
	∇p = szeros(T,6,3) # collider origin
	return [∇cf ∇p]
end

function contact_constraint_jacobian_contact_data(mechanism::Mechanism,
		contact::SoftContactConstraint{T,N,Nc,Cs}, body::Body{T}) where {T,N,Nc,Cs<:SoftContact{T,N}}
	timestep = mechanism.timestep
	model = contact.model
    collision = model.collision
	cf = [collision.collider.options.sliding_friction]
	p = collision.collider_origin
	xp, vp, qp, ϕp = next_configuration_velocity(body.state, timestep)
	xc, vc, qc, ϕc = next_configuration_velocity(mechanism.origin.state, timestep)
	function set_cf!(model, cf)
		m = deepcopy(model)
		m.collision.collider.options.sliding_friction = cf
		return m
	end
	function set_p!(model, p)
		m = deepcopy(model)
		m.collision.collider_origin = p
		return m
	end

	∇cf = FiniteDiff.finite_difference_jacobian(
		cf -> constraint(set_cf!(model, cf[1]), xp, vp, qp, ϕp, xc, vc, qc, ϕc, timestep; recompute=true),
		cf)

	∇p = FiniteDiff.finite_difference_jacobian(
		p -> constraint(set_p!(model, p), xp, vp, qp, ϕp, xc, vc, qc, ϕc, timestep; recompute=true),
		p)
	return [∇cf ∇p]
end

function contact_constraint_jacobian_body_data(mechanism::Mechanism,
		contact::SoftContactConstraint{T,N,Nc,Cs}, body::Body{T}) where {T,N,Nc,Cs}
    Nd = data_dim(body)
	# [m; j; v15; ϕ15; x2; vector(q2)]
	∇ = szeros(T,6,13)

	timestep = mechanism.timestep
	model = contact.model
    collision = model.collision
	xp, vp, qp, ϕp = next_configuration_velocity(body.state, timestep)
	xc, vc, qc, ϕc = next_configuration_velocity(mechanism.origin.state, timestep)
	x2p = current_position(body.state)
	q2p = current_orientation(body.state)

	∇x2p = FiniteDiff.finite_difference_jacobian(
		x2p -> constraint(model, next_position(x2p, vp, timestep), vp, qp, ϕp, xc, vc, qc, ϕc, timestep; recompute=true),
		x2p)
	∇q2p = FiniteDiff.finite_difference_jacobian(
		q2p -> constraint(model, xp, vp, next_orientation(Quaternion(q2p...), ϕp, timestep), ϕp, xc, vc, qc, ϕc, timestep; recompute=true),
		vector(q2p)) * LVᵀmat(q2p)
	return [∇ ∇x2 ∇q2]
end

################################################################################
# Analytical Jacobian
################################################################################
# Controller
function ctrl!(mechanism, k)
	nu = Dojo.input_dimension(mechanism)
	if Dojo.input_dimension(mechanism.joints[1]) == 6
		u = 0.2 * [szeros(6); mechanism.timestep * sones(nu-6)]
	else
		u = 0.2 * mechanism.timestep * sones(nu)
	end
	Dojo.set_input!(mechanism, u)
	return
end

function test_data_system(model::Symbol;
		ϵ=1.0e-6,
		tsim=0.1,
		ctrl=(m, k)->nothing,
        timestep=0.01,
		gravity=[0.0; 0.0; -9.81],
		verbose=false,
		T=Float64,
		kwargs...)

    # mechanism
    mechanism = Dojo.get_mechanism(model,
		timestep=timestep,
		gravity=gravity;
		kwargs...)
    Dojo.initialize!(mechanism, model)
    # simulate
    Dojo.simulate!(mechanism, tsim, ctrl!,
        record=false,
		verbose=false,
		opts=Dojo.SolverOptions(rtol=ϵ, btol=ϵ))

	# Finite Difference
	Nd = Dojo.data_dim(mechanism,
		attjac=false)
	data0 = Dojo.get_data(mechanism)
	sol0 = Dojo.get_solution(mechanism)
	datajac0 = Dojo.finite_difference_data_jacobian(mechanism, data0, sol0)
	attjac0 = Dojo.data_attitude_jacobian(mechanism)
	datajac0 *= attjac0

	# Analytical
	D = Dojo.create_data_matrix(mechanism.joints, mechanism.bodies, mechanism.contacts)
	Dojo.jacobian_data!(D, mechanism)
	nodes = [mechanism.joints; mechanism.bodies; mechanism.contacts]
	dimrow = length.(nodes)
	dimcol = Dojo.data_dim.(nodes)
	datajac1 = Dojo.full_matrix(D, dimrow, dimcol)

	# Test
	norm(datajac0 - datajac1, Inf)
	return datajac0, datajac1
end

test_data_system(:bunny)
get_data(mech.contacts[1].model)

J0, J1 = test_solmat(:bunny, tsim=0.3)

plot(Gray.(abs.(J0)))
plot(Gray.(abs.(J1)))
plot(Gray.(100*abs.(J1 + J0)))
