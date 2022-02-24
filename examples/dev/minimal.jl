













function minimal_velocities(joint::Translational, xa::AbstractVector,
        va::AbstractVector,  qa::UnitQuaternion, ωa::AbstractVector,
        xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ωb::AbstractVector)
	vertices = joint.vertices
    pbcb_w = vrotate(-vertices[2], qb)
    pbca_w = xa - (xb + vrotate(vertices[2], qb))
    # Δvw = V(pb,B/A)w - V(pa,A/A)w
	@show nullspace_mask(joint)
    Δvw = vb + skew(pbcb_w) * vrotate(ωb, qb) - (va + skew(pbca_w) * vrotate(ωa, qa)) # in world frame
	Δv = vrotate(Δvw, inv(qa)) # in the a frame
	@show Δv
	@show Δvw
	# @show pbcb_w
	# @show pbca_w
	# @show + skew(pbcb_w) * vrotate(ωb, qb) - (va + skew(pbca_w) * vrotate(ωa, qa))
    return nullspace_mask(joint) * Δv
end

function set_minimal_velocities(joint::Translational, xa::AbstractVector,
        va::AbstractVector, qa::UnitQuaternion, ωa::AbstractVector,
        xb::AbstractVector, qb::UnitQuaternion, ωb::AbstractVector;
        Δv::AbstractVector=szeros(control_dimension(joint)))
    # Δv is expressed in along the joint's nullspace axes, in pnode's frame

    vertices = joint.vertices
    pbcb_w = vrotate(-vertices[2], qb)
    pbca_w = xa - (xb + vrotate(vertices[2], qb))
    Aᵀ = zerodimstaticadjoint(nullspace_mask(joint))
    Δv = Aᵀ * Δv # in pnode's frame
    Δvw = vrotate(Δv, qa) # in the world frame
    # Δvw = V(pb,B/A)w - V(pa,A/A)w
    vb = Δvw - skew(pbcb_w) * vrotate(ϕb, qb) + (va + skew(pbca_w) * vrotate(ϕa, qa)) # in world frame
	@show Δv
	@show Δvw
	# @show pbcb_w
	# @show pbca_w
	# @show - skew(pbcb_w) * vrotate(ϕb, qb) + (va + skew(pbca_w) * vrotate(ϕa, qa))
    return vb
end


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

xb = storage.x[1][4]
qb = storage.q[1][4]
vb = storage.vl[1][4]
ϕb = storage.ωl[1][4]
# vb = storage.x[1].vsol[2]
# ϕb = storage.x[1].ϕsol[2]
#
# Δv = minimal_velocities(tra, xa, va, qa, ϕa, xb, vb, qb, ϕb)
# vb2 = set_minimal_velocities(tra, xa, va, qa, ϕa, xb, qb, ϕb; Δv=Δv)
# vb2
# vb
#
# mechanism.bodies[1].state


z = get_maximal_state(mechanism)
x = get_minimal_state(mechanism)
u = zeros(control_dimension(mechanism))
norm(minimal_to_maximal(mechanism, x) - z, Inf)
norm(maximal_to_minimal(mechanism, z) - x, Inf)

(minimal_to_maximal(mechanism, x) - z)[[1:3; 7:13]]

z[4:6]
minimal_to_maximal(mechanism, x)[4:6]

set_minimal_coordinates_velocities!(mechanism, mechanism.joints[1], xmin=x)
get_maximal_state(mechanism)
Δx0 = minimal_coordinates(mechanism.joints[1].translational, mechanism.origin, mechanism.bodies[1])
Δθ0 = minimal_coordinates(mechanism.joints[1].rotational, mechanism.origin, mechanism.bodies[1])
Δv0 = minimal_velocities(mechanism.joints[1].translational, mechanism.origin, mechanism.bodies[1])
Δϕ0 = minimal_velocities(mechanism.joints[1].rotational, mechanism.origin, mechanism.bodies[1])
xmin0 = [Δx0; Δθ0; Δv0; Δϕ0]
set_minimal_coordinates_velocities!(mechanism, mechanism.joints[1], xmin=xmin0)
Δx1 = minimal_coordinates(mechanism.joints[1].translational, mechanism.origin, mechanism.bodies[1])
Δθ1 = minimal_coordinates(mechanism.joints[1].rotational, mechanism.origin, mechanism.bodies[1])
Δv1 = minimal_velocities(mechanism.joints[1].translational, mechanism.origin, mechanism.bodies[1])
Δϕ1 = minimal_velocities(mechanism.joints[1].rotational, mechanism.origin, mechanism.bodies[1])
xmin1 = [Δx1; Δθ1; Δv1; Δϕ1]



storage.v[1]
storage.ω[1]




















function minimal_velocities(joint::JointConstraint, xa::AbstractVector,
        va::AbstractVector, qa::UnitQuaternion, ωa::AbstractVector,
        xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ωb::AbstractVector, timestep)
	rot = joint.rotational
	tra = joint.translational
	pa = tra.vertices[1]
    pb = tra.vertices[2]
	qoffset = rot.qoffset
	Arot = nullspace_mask(rot)
	Atra = nullspace_mask(tra)

	# Coordinates
	q = inv(qoffset) * inv(qa) * qb
    Δθ = Arot * rotation_vector(q)
    Δx = Atra * displacement(tra, xa, qa, xb, qb)

	xa10 = next_position(xa, -va, timestep)
	qa10 = next_orientation(qa, -ϕa, timestep)
	xb10 = next_position(xb, -vb, timestep)
	qb10 = next_orientation(qb, -ϕb, timestep)
	# Previous step coordinates
	q10 = inv(qoffset) * inv(qa10) * qb10
    Δθ10 = Arot * rotation_vector(q10)
    Δx10 = Atra * displacement(tra, xa10, qa10, xb10, qb10)

	# Finite difference
	Δϕ = (Δθ - Δθ10) / timestep
	Δv = (Δx - Δx10) / timestep
	return Δx, Δθ, Δv, Δϕ
end

function set_minimal_velocities(joint::JointConstraint, xa::AbstractVector,
        va::AbstractVector, qa::UnitQuaternion, ωa::AbstractVector, timestep;
		Δx=szeros(control_dimension(joint)),
		Δθ=szeros(control_dimension(joint)),
		Δv=szeros(control_dimension(joint)),
		Δϕ=szeros(control_dimension(joint)),
		)
	rot = joint.rotational
	tra = joint.translational
	pa = tra.vertices[1]
    pb = tra.vertices[2]
	qoffset = rot.qoffset
	Arotᵀ = zerodimstaticadjoint(nullspace_mask(rot))
	Atraᵀ = zerodimstaticadjoint(nullspace_mask(tra))

	# Δθ is expressed in along the joint's nullspace axes, in pnode's offset frame
    Δq = axis_angle_to_quaternion(Arotᵀ*Δθ)
    qb = qa * qoffset * Δq
    # Δx is expressed in along the joint's nullspace axes, in pnode's frame
    xb = xa + vrotate(pa + Atraᵀ*Δx, qa) - vrotate(pb, qb)

	xa10 = next_position(xa, -va, timestep)
	qa10 = next_orientation(qa, -ϕa, timestep)

	# Finite difference
	Δx10 = Δx - Δv * timestep
	Δθ10 = Δθ - Δϕ * timestep

	# Δθ is expressed in along the joint's nullspace axes, in pnode's offset frame
    Δq10 = axis_angle_to_quaternion(Arotᵀ*Δθ10)
    qb10 = qa10 * qoffset * Δq10
    # Δx is expressed in along the joint's nullspace axes, in pnode's frame
    xb10 = xa10 + vrotate(pa + Atraᵀ*Δx10, qa10) - vrotate(pb, qb10)

	# Finite difference
	vb = (xb - xb10) / timestep
	ϕb = angular_velocity(qb10, qb, timestep)

	return xb, vb, qb, ϕb
	# set_maximal_configuration!(cnode; x=xb, q=qb)
    # set_maximal_velocity!(cnode; v=vb, ω=ϕb)
	# return nothing
end


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

# xb = storage.x[1][4]
# qb = storage.q[1][4]
# vb = storage.vl[1][4]
# ϕb = storage.ωl[1][4]
# # vb = storage.x[1].vsol[2]
# # ϕb = storage.x[1].ϕsol[2]
#

Δx, Δθ, Δv, Δϕ = minimal_velocities(joint, xa, va, qa, ϕa, xb, vb, qb, ϕb, timestep)
Δx
Δv
Δθ = Δθ
Δϕ = Δϕ
xb2, vb2, qb2, ϕb2 = set_minimal_velocities(joint, xa, va, qa, ϕa, timestep;
	Δx=Δx, Δθ=Δθ, Δv=Δv, Δϕ=Δϕ)
norm(xb - xb2, Inf)
norm(vb - vb2, Inf)
norm(ϕb - ϕb2, Inf)
norm(vector(qb) - vector(qb2), Inf)


Arotᵀ = zerodimstaticadjoint(nullspace_mask(joint.rotational))
Δq = axis_angle_to_quaternion(Arotᵀ*Δθ)
