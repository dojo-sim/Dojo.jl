function sdf(contact::ContactConstraint{T,N,Nc,Cs}, x::AbstractVector{T},
    q::UnitQuaternion{T}) where {T,N,Nc,Cs<:Tuple{<:Contact{T,N}}}
    cont = contact.model
    return cont.surface_normal_projector * (x + vrotate(cont.p, q) - cont.offset)
end

function get_sdf(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, storage::Storage{T,N}) where {T,Nn,Ne,Nb,Ni,N}
    d = []
    for contact in mechanism.contacts
        ibody = get_body(mechanism, contact.parent_id).id - Ne
        push!(d, [sdf(contact, storage.x[ibody][i], storage.q[ibody][i]) for i = 1:N])
    end
    return d
end

function contact_location(mechanism::Mechanism)
    return [contact_location(mech, contact) for contact in mechanism.contacts]
end

function contact_location(mechanism::Mechanism, contact::ContactConstraint)
    body = mechanism.bodies[findfirst(x -> x.id == contact.parent_id, mechanism.bodies)]
    return contact_location(contact, body)
end

function contact_location(contact::ContactConstraint{T,N,Nc,Cs},
    body::Body) where {T,N,Nc,Cs<:Tuple{<:Contact{T,N}}}
    x = body.state.x2[1]
    q = body.state.q2[1]
    return contact_location(contact, x, q)
end

function contact_location(contact::ContactConstraint{T,N,Nc,Cs}, x::AbstractVector{T},
    q::UnitQuaternion{T}) where {T,N,Nc,Cs<:Tuple{<:Contact{T,N}}}
    cont = contact.model
    return x + vrotate(cont.p,q) - cont.offset
end
