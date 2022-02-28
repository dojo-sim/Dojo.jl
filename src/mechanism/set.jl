# maximal
""" 
    set_maximal_state(mechanism, z) 

    set the maximal state of a mechanism 

    mechanism: Mechanism 
    z: state 
"""
function set_maximal_state!(mechanism::Mechanism, z::AbstractVector)
    off = 0
    for body in mechanism.bodies
        x2, v15, q2, ϕ15 = unpack_data(z[off+1:end]); off += 13
        q2 = UnitQuaternion(q2..., false)
        body.state.v15 = v15
        body.state.ϕ15 = ϕ15
        body.state.x2 = x2
        body.state.q2 = q2
		initialize_state!(mechanism) # set x1, q1 and zeroes out F2 τ2
    end
	# warm-start solver
	for body in mechanism.bodies
		set_velocity_solution!(body)
	end
end

function initialize_state!(mechanism::Mechanism)
    for body in mechanism.bodies initialize_state!(body, mechanism.timestep) end
end

# inputs
""" 
    set_input(mechanism, u) 

    set input for each joint in mechanism 

    mechanism: Mechanism 
    u: input 
"""
function set_input!(mechanism::Mechanism{T}, u::AbstractVector) where T
	joints = mechanism.joints
	# set the controls in the equality constraints
	off = 0
	for joint in joints
		nu = input_dimension(joint)
		set_input!(joint, SVector{nu,T}(u[off .+ (1:nu)]))
		off += nu
	end
	# apply the controls to each body's state
	for joint in joints
		input_impulse!(joint, mechanism)
	end
end

function set_input!(mechanism::Mechanism, dict::Dict)
    for (id, joint) in pairs(mechanism.joints)
        set_input!(joint, dict[id])
    end
end

# minimal
""" 
    set_minimal_state(mechanism, y) 

    set the maximal state of a mechanism 

    mechanism: Mechanism 
    y: state 
"""
function set_minimal_state!(mechanism::Mechanism, y::AbstractVector)
    z = minimal_to_maximal(mechanism, y) 
    set_maximal_state!(mechanism, z)
end
    
function set_minimal_coordinates!(mechanism::Mechanism, dict::Dict)
    for (id, joint) in pairs(mechanism.joints)
        set_minimal_coordinates!(mechanism, joint, dict[id])
    end
end

function set_minimal_velocities!(mechanism::Mechanism, dict::Dict)
    for (id, joint) in pairs(mechanism.joints)
        set_minimal_velocities!(mechanism, joint, dict[id])
    end
end

# velocity
""" 
    zero_velocity!(mechanism) 

    set all mechanism body velocities to zero 

    mechanism: Mechanism 
"""
function zero_velocity!(mechanism::Mechanism)
    for (i, body) in enumerate(mechanism.bodies)
        try
            set_maximal_velocities!(body, v=zeros(3), ω=zeros(3))
            set_previous_configuration!(body, mechanism.timestep)
        catch
            nothing
        end
    end
end

# springs + dampers
function set_spring_offset!(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, x::AbstractVector) where {T,Nn,Ne,Nb,Ni}
	off = 0
	for id in mechanism.root_to_leaves
		(id > Ne) && continue # only treat joints
		joint = mechanism.joints[id]
        N̄ = 3 - length(joint)
        joint.spring_offset = x[off .+ (1:N̄)]
        off += 2N̄
	end
	return nothing
end