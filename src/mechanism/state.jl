"""
	minimal_to_maximal(mechanism, y)

	convert minimal to maximal representation

	mechanism: Mechanism
	y: minimal state
"""
function minimal_to_maximal(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, y::AbstractVector) where {T,Nn,Ne,Nb,Ni}
	temp_mechanism = deepcopy(mechanism)
	nus = input_dimensions(mechanism)

	for id in temp_mechanism.root_to_leaves
		(id > Ne) && continue # only treat joints
		joint = temp_mechanism.joints[id]
		nu = input_dimension(joint)
		off = sum(nus[SUnitRange(1,id-1)])*2
		set_minimal_coordinates_velocities!(temp_mechanism, joint, xmin=y[SUnitRange(off+1,off+2nu)]) # TODO does this actually set a state and not just convert min to max?
	end

	return get_maximal_state(temp_mechanism)
end

function joint_indices(mechanism::Mechanism{T,Nn,Ne}) where {T,Nn,Ne}
    ind = []
    off = 0
    for joint in mechanism.joints
        nu = input_dimension(joint)
        push!(ind, Vector(off + nu .+ (1:nu)))
        off += 2nu
    end
    
    return vcat(ind...)
end

"""
	maximal_to_minimal(mechanism, z)

	convert maximal to minimal representation

	mechanism: Mechanism
	z: maximal state
"""
function maximal_to_minimal(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, z::AbstractVector{Tz}) where {T,Nn,Ne,Nb,Ni,Tz}
	x = []
	for joint in mechanism.joints
		c = zeros(Tz,0)
		v = zeros(Tz,0)
		i_child = joint.child_id - Ne
		for element in (joint.translational, joint.rotational)
			xb, vb, qb, ωb = unpack_maximal_state(z, i_child)
			if joint.parent_id != 0
				i_parent = joint.parent_id - Ne
				xa, va, qa, ωa = unpack_maximal_state(z, i_parent)
			else
				xa, va, qa, ωa = current_configuration_velocity(mechanism.origin.state)
			end
			push!(c, minimal_coordinates(element, xa, qa, xb, qb)...)
			push!(v, minimal_velocities(element, xa, va, qa, ωa, xb, vb, qb, ωb, mechanism.timestep)...)
		end
		push!(x, [c; v]...)
	end
	x = [x...]
	
	return x
end

function unpack_maximal_state(z::AbstractVector, i::Int)
	zi = z[SUnitRange((i-1)*13+1,i*13)]
	x2 = zi[SUnitRange(1,3)]
	v15 = zi[SUnitRange(4,6)]
	q2 = Quaternion(zi[SUnitRange(7,10)]...)
	ω15 = zi[SUnitRange(11,13)]
	return x2, v15, q2, ω15
end

function pack_maximal_state!(z::AbstractVector, 
	x2::AbstractVector, v15::AbstractVector, q2::Quaternion, ω15::AbstractVector, 
	i::Int)

	z[SUnitRange((i-1)*13+1,(i-1)*13+3)] = x2
	z[SUnitRange((i-1)*13+4,(i-1)*13+6)] = v15
	z[SUnitRange((i-1)*13+7,(i-1)*13+10)] = vector(q2)
	z[SUnitRange((i-1)*13+11,(i-1)*13+13)] = ω15

	return nothing
end
