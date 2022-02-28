""" 
    get_sdf(contact, x, q) 

    returns the signed distance for a contact

    contact: ContactConstraint
    x: body position 
    q: body orientation
"""
function get_sdf(contact::ContactConstraint{T,N,Nc,Cs}, x::AbstractVector{T},
    q::UnitQuaternion{T}) where {T,N,Nc,Cs<:Contact{T,N}}
    model = contact.model
    return model.surface_normal_projector * (x + vector_rotate(model.contact_point, q) - model.offset)
end

function get_sdf(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, storage::Storage{T,N}) where {T,Nn,Ne,Nb,Ni,N}
    d = []
    for contact in mechanism.contacts
        ibody = get_body(mechanism, contact.parent_id).id - Ne
        push!(d, [get_sdf(contact, storage.x[ibody][i], storage.q[ibody][i]) for i = 1:N])
    end
    return d
end

""" 
    contact_location(contact, x, q) 

    location of contact point in world coordinates 

    contact: ContactConstraint 
    x: body position 
    q: body orientation 
"""
function contact_location(contact::ContactConstraint{T,N,Nc,Cs}, x::AbstractVector{T},
    q::UnitQuaternion{T}) where {T,N,Nc,Cs<:Contact{T,N}}
    model = contact.model
    return x + vector_rotate(model.contact_point,q) - model.offset
end

function contact_location(contact::ContactConstraint{T,N,Nc,Cs},
    body::Body) where {T,N,Nc,Cs<:Contact{T,N}}
    x = body.state.x2
    q = body.state.q2
    return contact_location(contact, x, q)
end

function contact_location(mechanism::Mechanism, contact::ContactConstraint)
    body = mechanism.bodies[findfirst(x -> x.id == contact.parent_id, mechanism.bodies)]
    return contact_location(contact, body)
end

function contact_location(mechanism::Mechanism)
    return [contact_location(mech, contact) for contact in mechanism.contacts]
end



