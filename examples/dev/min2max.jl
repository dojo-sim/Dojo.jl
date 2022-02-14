
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
Δvϕ - [Δv; Δϕ]

minimal_velocities_jacobian_configuration_new(:parent,
	tra, xa, va, qa, ϕa, xb, vb, qb, ϕb, timestep)
minimal_velocities_jacobian_configuration_new(:child,
	tra, xa, va, qa, ϕa, xb, vb, qb, ϕb, timestep)
minimal_velocities_jacobian_velocity_new(:parent,
	tra, xa, va, qa, ϕa, xb, vb, qb, ϕb, timestep)
minimal_velocities_jacobian_velocity_new(:child,
	tra, xa, va, qa, ϕa, xb, vb, qb, ϕb, timestep)
