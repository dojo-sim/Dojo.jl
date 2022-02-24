function minimal_to_maximal(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, x::AbstractVector{Tx}) where {T,Nn,Ne,Nb,Ni,Tx}
	off = 0
	for id in mechanism.root_to_leaves
		(id > Ne) && continue # only treat joints
		joint = mechanism.joints[id]
		nu = control_dimension(joint)
		set_minimal_coordinates_velocities!(mechanism, joint, xmin=x[off .+ SUnitRange(1, 2nu)])
		off += 2nu
	end
	z = get_maximal_state(mechanism)
	return z
end

function maximal_to_minimal(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, z::AbstractVector{Tz}) where {T,Nn,Ne,Nb,Ni,Tz}
	x = []
	for id in mechanism.root_to_leaves
		(id > Ne) && continue # only treat joints
		joint = mechanism.joints[id]
		c = zeros(Tz,0)
		v = zeros(Tz,0)
		ichild = joint.child_id - Ne
		for element in [joint.translational, joint.rotational]
			xb, vb, qb, ϕb = unpack_maximal_state(z, ichild)
			if joint.parent_id != 0
				iparent = joint.parent_id - Ne
				xa, va, qa, ϕa = unpack_maximal_state(z, iparent)
			else
				xa, va, qa, ϕa = current_configuration_velocity(mechanism.origin.state)
			end
			push!(c, minimal_coordinates(element, xa, qa, xb, qb)...)
			push!(v, minimal_velocities(element, xa, va, qa, ϕa, xb, vb, qb, ϕb, mechanism.timestep)...)
		end
		push!(x, [c; v]...)
	end
	x = [x...]
	return x
end

function unpack_maximal_state(z::AbstractVector, i::Int)
	zi = z[(i-1) * 13 .+ (1:13)]
	x2 = zi[SUnitRange(1,3)]
	v15 = zi[SUnitRange(4,6)]
	q2 = UnitQuaternion(zi[7:10]..., false)
	ϕ15 = zi[SUnitRange(11,13)]
	return x2, v15, q2, ϕ15
end

function pack_maximal_state!(z::AbstractVector, 
	x2::AbstractVector, v15::AbstractVector, q2::UnitQuaternion, ϕ15::AbstractVector, 
	i::Int)

	z[(i-1) * 13 .+ (1:3)] = x2
	z[(i-1) * 13 .+ (4:6)] = v15
	z[(i-1) * 13 .+ (7:10)] = vector(q2)
	z[(i-1) * 13 .+ (11:13)] = ϕ15

	return nothing
end
