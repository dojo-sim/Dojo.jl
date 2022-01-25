function sdf(ineqc::ContactConstraint{T,N,Nc,Cs}, x::AbstractVector{T},
    q::UnitQuaternion{T}) where {T,N,Nc,Cs<:Tuple{<:Contact{T,N}}}
    cont = ineqc.constraints[1]
    return cont.ainv3 * (x + vrotate(cont.p, q) - cont.offset)
end

function get_sdf(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, storage::Storage{T,N}) where {T,Nn,Ne,Nb,Ni,N}
    d = []
    for ineqc in mechanism.ineqconstraints
        ibody = get_body(mechanism, ineqc.parentid).id - Ne
        push!(d, [sdf(ineqc, storage.x[ibody][i], storage.q[ibody][i]) for i = 1:N])
    end
    return d
end

function contact_location(mechanism::Mechanism)
    return [contact_location(mech, ineqc) for ineqc in mechanism.ineqconstraints]
end

function contact_location(mechanism::Mechanism, ineqc::ContactConstraint)
    body = mechanism.bodies[findfirst(x -> x.id == ineqc.parentid, mechanism.bodies)]
    return contact_location(ineqc, body)
end

function contact_location(ineqc::ContactConstraint{T,N,Nc,Cs},
    body::Body) where {T,N,Nc,Cs<:Tuple{<:Contact{T,N}}}
    x = body.state.x2[1]
    q = body.state.q2[1]
    return contact_location(ineqc, x, q)
end

function contact_location(ineqc::ContactConstraint{T,N,Nc,Cs}, x::AbstractVector{T},
    q::UnitQuaternion{T}) where {T,N,Nc,Cs<:Tuple{<:Contact{T,N}}}
    cont = ineqc.constraints[1]
    return x + vrotate(cont.p,q) - cont.offset
end
