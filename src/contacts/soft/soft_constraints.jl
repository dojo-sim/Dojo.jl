# impulses
function impulses!(mechanism::Mechanism{T}, body::Body, contact::SoftContactConstraint) where T
    body.state.d -= impulse_map(mechanism, contact, body) * contact.impulses[2]
    return
end

function impulse_map(mechanism, contact::SoftContactConstraint, body::Body)
    relative = (body.id == contact.parent_id ? :parent : :child)
    pbody = get_body(mechanism, contact.parent_id)
    cbody = get_body(mechanism, contact.child_id)
    return impulse_map(relative, contact.model, pbody, cbody, mechanism.timestep)
end

function impulse_map_jacobian_configuration(mechanism, body::Body{T}, contact::SoftContactConstraint{T}) where T
    relative = (body.id == contact.parent_id ? :parent : :child)

    return impulse_map_jacobian(relative, relative, contact.model,
        get_body(mechanism, contact.parent_id),
        get_body(mechanism, contact.child_id),
        contact.impulses[2],
        mechanism.timestep)
end

function impulses_jacobian_velocity!(mechanism, body::Body, contact::SoftContactConstraint)
    timestep = mechanism.timestep
    body.state.D -= impulse_map_jacobian_configuration(mechanism, body, contact) * integrator_jacobian_velocity(body, timestep)
    return
end

function impulse_map(relative::Symbol, model::SoftContact{T}, pbody::Node, cbody::Node, timestep) where T
    # configurations
    x2p, q2p = current_configuration(pbody.state)
    x2c, q2c = current_configuration(cbody.state)

    # find barycenter
	collision = model.collision
    ψ, barycenter, contact_normal = overlap(collision, x2p, q2p, x2c, q2c)

    # mapping
    # X = force_mapping(relative, model, xp, qp, xc, qc)
    X = force_mapping(relative, model, x2p, q2p, x2c, q2c)
    # # c = contact_point(relative, model.collision, xp, qp, xc, qc)
    # c = x2p + vector_rotate(model.contact_origin + barycenter, q2p)
    # # @show scn.(c)
    # if relative == :parent
    #     r = c - xp
    #     Q = rotation_matrix(inv(q2p)) * skew(r) * X
    # elseif relative == :child
    #     r = c - xc
    #     Q = rotation_matrix(inv(q2c)) * skew(r) * X
    # end

    Q = 1.0*skew(model.contact_origin + barycenter) * rotation_matrix(inv(q2p)) * X
    return [
                X szeros(T,3,3);
                Q X;
           ]
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
function initialize!(contact::SoftContactConstraint)
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
    # @show scn(norm(vector_entry.value))
    # @show scn.(vector_entry.value[1:3])
    # @show scn.(contact.impulses[2][1:3])
    contact.impulses[2] = contact.impulses[1] + 1 / (2^scale) * α * vector_entry.value
    return
end


# update
function update!(contact::SoftContactConstraint)
    contact.impulses[1] = contact.impulses[2]
    return
end

# get data
get_data(model::SoftContact{T}) where T = [model.collision.contact_radius; model.collision.contact_origin]
function set_data!(model::SoftContact, data::AbstractVector)
	# model.collider.options.sliding_friction = data[1]
    model.collision.contact_radius = data[1]
    model.collision.contact_origin = data[SVector{3,Int}(2:4)]
    return nothing
end

# constructor
function soft_contact_constraint(body::Body{T},
        normal::AbstractVector{T},
        collider::Collider;
        friction_coefficient::T=1.0,
        contact_origin::AbstractVector{T}=szeros(T, 3),
        contact_radius::T=0.0,
        name::Symbol=Symbol("contact_" * randstring(4))) where T

    model = SoftContact(body, normal, collider,
        contact_origin=contact_origin,
        contact_radius=contact_radius)
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
        color=RGBA(0.9,0.9,0.9,0.3))
    shape_vec = Vector{Shape{T}}([shape1, shape2])
    shapes = Shapes(shape_vec)

    return Body(collider.mass, collider.inertia; name=name, shape=shapes)
end


coulomb_direction(v, smoothing=1e3, regularizer=1e-3) = - atan(smoothing * norm(v)) * v/(regularizer + norm(v))
function ∂coulomb_direction∂v(v, smoothing=1e3, regularizer=1e-3)
    ∇ = - 1 / (1 + smoothing^2 * v'*v) * smoothing * v/(regularizer + norm(v)) * v'/(norm(v)+1e-20) +
        - atan(smoothing * norm(v)) * (1/(norm(v) + regularizer) * Diagonal(sones(3)) - v*v' ./ ((norm(v)+1e-20) * (norm(v) + regularizer)^2))
    return ∇
end

v = srand(3)
coulomb_direction(v)
J0 = ∂coulomb_direction∂v(v)
J1 = FiniteDiff.finite_difference_jacobian(v-> coulomb_direction(v), v)
norm(J0 - J1, Inf) < 1e-5
