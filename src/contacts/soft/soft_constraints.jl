function impulse_map(relative::Symbol, model::SoftContact{T}, pbody::Node, cbody::Node, timestep) where T
    # configurations
    x2p, q2p = current_configuration(pbody.state)
    x2c, q2c = current_configuration(cbody.state)

    # find barycenter
	collision = model.collision
	ψ, barycenter, normal = model.ψ, model.barycenter, model.normal

    # mapping
    XF = force_mapping(relative, model, x2p, q2p, x2c, q2c)
    if relative == :parent
		QF = skew(parent_origin(model.collision) + barycenter) * rotation_matrix(inv(q2p)) * XF
		Qτ = rotation_matrix(inv(q2p))
    elseif relative == :child
		p = x2p + Dojo.vector_rotate(parent_origin(model.collision) + barycenter, q2p)
		pc = x2c + Dojo.vector_rotate(child_origin(model.collision), q2c)
		QF = rotation_matrix(inv(q2c)) * skew(p - pc) * XF
		Qτ = -rotation_matrix(inv(q2c))
    end
	Z = szeros(T,3,3)
	return [
		[XF Z];
		[QF Qτ]]
end

function impulse_map_jacobian(relative::Symbol, jacobian::Symbol, model::SoftContact{T},
        pbody::Node, cbody::Node, λ, timestep) where T
    return szeros(T,6,7)
end

# off-diagonal terms for linear system
function off_diagonal_jacobians(mechanism, body::Body, contact::SoftContactConstraint{T,N,Nc,Cs}) where {T,N,Nc,Cs}
    return -impulse_map(mechanism, contact, body), constraint_jacobian_velocity(mechanism, contact, body)
end

# off-diagonal terms for linear system
function off_diagonal_jacobians(mechanism, contact::SoftContactConstraint{T,N,Nc,Cs}, body::Body) where {T,N,Nc,Cs}
    return constraint_jacobian_velocity(mechanism, contact, body), -impulse_map(mechanism, contact, body)
end

# linear system entries
function set_matrix_vector_entries!(mechanism::Mechanism, matrix_entry::Entry, vector_entry::Entry, contact::SoftContactConstraint)
    matrix_entry.value = constraint_jacobian(contact)
    vector_entry.value = -constraint(mechanism, contact)
    return
end

# reset variables using cone-specific neutral vector
function reset!(contact::SoftContactConstraint{T,N}; scale=0.0) where {T,N}
    contact.impulses[1] = szeros(T,N)
    contact.impulses[2] = szeros(T,N)
    return
end

# initialization
function initialize!(mechanism::Mechanism, contact::SoftContactConstraint)
	pbody = get_body(mechanism, contact.parent_id)
	cbody = get_body(mechanism, contact.child_id)
	x2p, q2p = current_configuration(pbody.state)
	x2c, q2c = current_configuration(cbody.state)
	model = contact.model
	collision = model.collision
	model.ψ, model.barycenter, model.normal = overlap(collision, x2p, q2p, x2c, q2c)
    return nothing
end

# complementarity
complementarity(mechanism, contact::SoftContactConstraint{T,N}; scaling=false) where {T,N} = szeros(T,N)
complementarityμ(mechanism, contact::SoftContactConstraint{T,N}; scaling=false) where {T,N} = szeros(T,N)

# cone line search
function cone_line_search!(α, mechanism, contact::SoftContactConstraint,
        vector_entry::Entry, τort, τsoc; scaling::Bool=false)
    return α
end

# centering
function centering!(ν, νaff, n, mechanism, contact::SoftContactConstraint, vector_entry::Entry, αaff)
    return ν, νaff, n
end

# candidate step
function candidate_step!(α::T, contact::SoftContactConstraint{T}, vector_entry::Entry, scale) where T
    contact.impulses[2] = contact.impulses[1] + 1 / (2^scale) * α * vector_entry.value
    return
end

# update
function update!(contact::SoftContactConstraint)
    contact.impulses[1] = contact.impulses[2]
    return
end

# get data
data_dim(model::SoftContact) = 5 # [bounciness, friction_coefficient, p]
function get_data(model::SoftContact{T}) where T
	options = model.collision.collider.options
	return [
		log(options.impact_damper/options.impact_spring);
		options.sliding_friction;
		model.collision.collider_origin]
end
function set_data!(model::SoftContact, data::AbstractVector)
	options = model.collision.collider.options
	options.impact_damper = exp(data[1]) * options.impact_spring
	options.sliding_friction = data[2]
    model.collision.collider_origin = data[SUnitRange(3,5)]
    return nothing
end

coulomb_direction(v, smoothing=1e3, regularizer=1e-3) = - atan(smoothing * norm(v)) * v/(regularizer + norm(v))
function ∂coulomb_direction∂v(v, smoothing=1e3, regularizer=1e-3)
    ∇ = - 1 / (1 + smoothing^2 * v'*v) * smoothing * v/(regularizer + norm(v)) * v'/(norm(v)+1e-20) +
        - atan(smoothing * norm(v)) * (1/(norm(v) + regularizer) * Diagonal(sones(3)) - v*v' ./ ((norm(v)+1e-20) * (norm(v) + regularizer)^2))
    return ∇
end

# constructor
function soft_contact_constraint(body::Body{T},
        normal::AbstractVector{T},
        collider::Collider;
        friction_coefficient::T=1.0,
        collider_origin::AbstractVector{T}=szeros(T, 3),
        name::Symbol=Symbol("contact_" * randstring(4))) where T

	collider.options.sliding_friction = friction_coefficient
    model = SoftContact(body, normal, collider,
        parent_origin=collider_origin)
    contact = SoftContactConstraint((model, body.id, 0); name=name)
    return contact
end

# constructor
function SoftBody(collider::Collider, inner_mesh_path::String, outer_mesh_path::String;
        position_offset::AbstractVector=szeros(3),
        axis_offset::Quaternion=one(Quaternion),
        scale::AbstractVector=sones(3),
        name::Symbol=Symbol("body_" * randstring(4)),
        color=RGBA(0.2,0.2,0.2,1.0))
    T = promote_type(quateltype.((position_offset, axis_offset))...)

    shape1 = Mesh(inner_mesh_path,
        position_offset = position_offset - collider.center_of_mass,
        axis_offset = axis_offset,
        color=color)
    shape2 = Mesh(outer_mesh_path,
        position_offset = position_offset - collider.center_of_mass,
        axis_offset = axis_offset,
        color=RGBA(0.9,0.9,0.9,0.5))
    shape_vec = Vector{Shape{T}}([shape1, shape2])
    shapes = Shapes(shape_vec)

    return Body(collider.mass, collider.inertia; name=name, shape=shapes)
end
