
# function chain_jacobian(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, joint::JointConstraint, pnode::Node,
# 		cnode::Node) where {T,Nn,Ne,Nb,Ni,Tx}
# 	xmin = minimal_coordinates_velocites(joint, pnode, cnode)
# 	pz =
# 	cz =
#
# 	return nothing
# end

s = Dict{Vector{Int}, Matrix}([1,1] => ones(3,3))
s[2] = 3
parent()





mech = get_quadruped()
initialize!(mech, :quadruped)
x = get_minimal_state(mech)
J = minimal_to_maximal_jacobian2(mech, x)
Main.@profiler minimal_to_maximal_jacobian2(mech, x)
@benchmark J = minimal_to_maximal_jacobian2(mech, x)
spy(J)
mech.system.dfs_list

parents(mech.system, 6)
parents(mech.system, 9)
parents(mech.system, 9)


pnode = mech.origin
cnode = mech.bodies[1]
# chain_jacobian(mech, pnode, cnode)




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
timestep= mech.timestep
Δvϕ = minimal_velocities(joint, xa, va, qa, ϕa, xb, vb, qb, ϕb, timestep)
Δv = minimal_velocities(tra, xa, va, qa, ϕa, xb, vb, qb, ϕb, timestep)
Δϕ = minimal_velocities(rot, xa, va, qa, ϕa, xb, vb, qb, ϕb, timestep)
Δvϕ1 = minimal_velocities(joint, xa, va, qa, ϕa, xb, vb, qb, ϕb, 10*timestep)
Δv1 = minimal_velocities(tra, xa, va, qa, ϕa, xb, vb, qb, ϕb, 10*timestep)
Δϕ1 = minimal_velocities(rot, xa, va, qa, ϕa, xb, vb, qb, ϕb, 0.01*timestep)
Δvϕ - [Δv; Δϕ]



minimal_velocities_jacobian_configuration(:parent,
	tra, xa, va, qa, ϕa, xb, vb, qb, ϕb, timestep)
minimal_velocities_jacobian_configuration(:child,
	tra, xa, va, qa, ϕa, xb, vb, qb, ϕb, timestep)
minimal_velocities_jacobian_velocity(:parent,
	tra, xa, va, qa, ϕa, xb, vb, qb, ϕb, timestep)
minimal_velocities_jacobian_velocity(:child,
	tra, xa, va, qa, ϕa, xb, vb, qb, ϕb, timestep)




function ctrl!(mech, k)
	set_input!(mech, 0.0*sones(input_dimension(mech))*mech.timestep)
end
# vis=visualizer()
# open(vis)
mechanism = get_mechanism(:pendulum, timestep=0.02, gravity=-0.0 * 9.81)
initialize!(mechanism, :pendulum, ϕ1 = 0.5 * π, ω1 = 2.0)
storage = simulate!(mechanism, 0.07, ctrl!, record=true, verbose=false)
visualize(mechanism, storage, vis=vis)


joint = mechanism.joints[1]
tra = joint.translational
timestep= mechanism.timestep
pbody = mechanism.origin
cbody = mechanism.bodies[1]

xa = pbody.state.x2
qa = pbody.state.q2
va = pbody.state.v15
ϕa = pbody.state.ϕ15
# va = pbody.state.vsol[2]
# ϕa = pbody.state.ϕsol[2]

xb = cbody.state.x2
qb = cbody.state.q2
vb = cbody.state.v15
ϕb = cbody.state.ϕ15
# vb = cbody.state.vsol[2]
# ϕb = cbody.state.ϕsol[2]

Δx = minimal_coordinates(joint.translational, xa, qa, xb, qb)
Δθ = minimal_coordinates(joint.rotational, xa, qa, xb, qb)
Δv = minimal_velocities(joint.translational, xa, va, qa, ϕa, xb, vb, qb, ϕb, 1e-5*timestep)
Δϕ = minimal_velocities(joint.rotational, xa, va, qa, ϕa, xb, vb, qb, ϕb, 1e-5*timestep)
vb2, ϕb2 = set_minimal_velocities(joint, xa, va, qa, ϕa, xb, qb, 1e-5*timestep, Δv=Δv, Δϕ=Δϕ)
norm(xb - xb2, Inf)
norm(vb - vb2, Inf)
norm(ϕb - ϕb2, Inf)
norm(vector(qb) - vector(qb2), Inf)


Arotᵀ = zerodimstaticadjoint(nullspace_mask(joint.rotational))
Δq = axis_angle_to_quaternion(Arotᵀ*Δθ)
