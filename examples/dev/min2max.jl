
function minimal_to_maximal_jacobian(mechanism::Mechanism{T,Nn,Ne,Nb,Ni},
		x::AbstractVector{Tx}; attjac::Bool=false) where {T,Nn,Ne,Nb,Ni,Tx}
	J = zeros(maximal_dimension(mechanism, attjac=false), minimal_dimension(mechanism))
	z = minimal_to_maximal(mechanism, x)
	off = 0
	for id in reverse(mechanism.system.dfs_list)
		(id > Ne) && continue # only treat joints
		joint = mechanism.joints[id]
		nu = control_dimension(joint)
		idx = collect(off .+ (1:(2nu)))

		J[:, idx] = position_velocity_jacobian(mechanism, joint, z, x[idx])

		off += 2nu
	end
	attjac && (J = attitude_jacobian(z, Nb)' * J)
	return J
end

function chain_jacobian(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, joint::JointConstraint, pnode::Node,
		cnode::Node) where {T,Nn,Ne,Nb,Ni,Tx}
	xmin = minimal_coordinates_velocites(joint, pnode, cnode)
	pz =
	cz =

	return nothing
end


mech = get_pendulum()
initialize!(mech, :pendulum)
pnode = mech.origin
cnode = mech.bodies[1]
chain_jacobian(mech, pnode, cnode)




joint = mech.joints[1]
tra = mech.joints[1].translational
rot = mech.joints[1].rotational
xa = srand(3)
va = srand(3)
qa = UnitQuaternion(rand(4)...)
ϕa = srand(3)
xb = srand(3)
vb = srand(3)
qb = UnitQuaternion(rand(4)...)
ϕb = srand(3)
timestep = mech.timestep
Δvϕ = minimal_velocities_new(joint, xa, va, qa, ϕa, xb, vb, qb, ϕb, timestep)
Δv = minimal_velocities_new(tra, xa, va, qa, ϕa, xb, vb, qb, ϕb, timestep)
Δϕ = minimal_velocities_new(rot, xa, va, qa, ϕa, xb, vb, qb, ϕb, timestep)
Δvϕ1 = minimal_velocities_new(joint, xa, va, qa, ϕa, xb, vb, qb, ϕb, 10*timestep)
Δv1 = minimal_velocities_new(tra, xa, va, qa, ϕa, xb, vb, qb, ϕb, 10*timestep)
Δϕ1 = minimal_velocities_new(rot, xa, va, qa, ϕa, xb, vb, qb, ϕb, 0.01*timestep)
Δvϕ - [Δv; Δϕ]

set_minimal_velocities_new!


minimal_velocities_jacobian_configuration_new(:parent,
	tra, xa, va, qa, ϕa, xb, vb, qb, ϕb, timestep)
minimal_velocities_jacobian_configuration_new(:child,
	tra, xa, va, qa, ϕa, xb, vb, qb, ϕb, timestep)
minimal_velocities_jacobian_velocity_new(:parent,
	tra, xa, va, qa, ϕa, xb, vb, qb, ϕb, timestep)
minimal_velocities_jacobian_velocity_new(:child,
	tra, xa, va, qa, ϕa, xb, vb, qb, ϕb, timestep)




function ctrl!(mech, k)
	set_control!(mech, 0.0*sones(control_dimension(mech))*mech.timestep)
end
# vis = Visualizer()
# open(vis)
mechanism = get_mechanism(:pendulum, timestep = 0.02, gravity = -0.0 * 9.81)
initialize!(mechanism, :pendulum, ϕ1 = 0.5 * π, ω1 = 2.0)
storage = simulate!(mechanism, 0.07, ctrl!, record = true, verbose = false)
visualize(mechanism, storage, vis=vis)


joint = mechanism.joints[1]
tra = joint.translational
timestep = mechanism.timestep
bodya = mechanism.origin
bodyb = mechanism.bodies[1]

xa = bodya.state.x2[1]
qa = bodya.state.q2[1]
va = bodya.state.v15
ϕa = bodya.state.ϕ15
# va = bodya.state.vsol[2]
# ϕa = bodya.state.ϕsol[2]

xb = bodyb.state.x2[1]
qb = bodyb.state.q2[1]
vb = bodyb.state.v15
ϕb = bodyb.state.ϕ15
# vb = bodyb.state.vsol[2]
# ϕb = bodyb.state.ϕsol[2]

Δx = minimal_coordinates(joint.translational, xa, qa, xb, qb)
Δθ = minimal_coordinates(joint.rotational, xa, qa, xb, qb)
Δv = minimal_velocities_new(joint.translational, xa, va, qa, ϕa, xb, vb, qb, ϕb, 1e-5*timestep)
Δϕ = minimal_velocities_new(joint.rotational, xa, va, qa, ϕa, xb, vb, qb, ϕb, 1e-5*timestep)
vb2, ϕb2 = set_minimal_velocities_new(joint, xa, va, qa, ϕa, xb, qb, 1e-5*timestep,
	Δx=Δx, Δθ=Δθ, Δv=Δv, Δϕ=Δϕ)
norm(xb - xb2, Inf)
norm(vb - vb2, Inf)
norm(ϕb - ϕb2, Inf)
norm(vector(qb) - vector(qb2), Inf)


Arotᵀ = zerodimstaticadjoint(nullspace_mask(joint.rotational))
Δq = axis_angle_to_quaternion(Arotᵀ*Δθ)
