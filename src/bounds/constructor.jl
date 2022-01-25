"""
    Generate contact inequality constraints attached to a list of bodies. You need to provide:
    - the normal for each contact point
    - the coefficient of friction for each contact point (optional for `ImpactBound`)
    - the offset vector p with respect to the center of the body for each contact point (optional)
    - the altitude offset for each contact point (optional)
    - the contact type: `:contact`, `:linear_contact`, `:impact`
"""
function contactconstraint(bodies::AbstractVector{<:Body{T}},
        normal::AbstractVector{<:AbstractVector};
        cf::AbstractVector{T} = ones(length(normal)),
        p::AbstractVector = [szeros(T, 3) for i=1:length(normal)],
        offset::AbstractVector = [szeros(T, 3) for i=1:length(normal)],
        names::Vector{Symbol} = [Symbol("contact_" * randstring(4)) for i = 1:length(normal)],
        contact_type::Symbol = :contact) where {T}

    n = length(normal)
    @assert n == length(bodies) == length(normal) == length(cf) == length(p) == length(offset)
    ineqcs = Vector{InequalityConstraint}()
    for i = 1:n
        ineqc = contactconstraint(bodies[i], normal[i], cf=cf[i], p=p[i],
            offset=offset[i], name=names[i], contact_type=contact_type)
        push!(ineqcs, ineqc)
    end
    ineqcs = [ineqcs...] # vector typing
    return ineqcs
end

function contactconstraint(body::Body{T},
        normal::AbstractVector{<:AbstractVector};
        cf::AbstractVector{T} = ones(length(normal)),
        p::AbstractVector = [szeros(T, 3) for i=1:length(normal)],
        offset::AbstractVector = [szeros(T, 3) for i=1:length(normal)],
        names::Vector{Symbol} = [Symbol("contact_" * randstring(4)) for i = 1:length(normal)],
        contact_type::Symbol = :contact) where {T}
    n = length(normal)
    @assert n == length(normal) == length(cf) == length(p) == length(offset)
    return contactconstraint(fill(body, n), normal, cf=cf, p=p, offset=offset,
        names=names, contact_type=contact_type)
end

"""
    Generate contact inequality constraint attached to one body. You need to provide:
    - the normal for the contact point
    - the coefficient of friction for the contact point
    - the offset vector p with respect to the center of the body for the contact point (optional)
    - the altitude offset for the contact point (optional)
"""
function contactconstraint(body::Body{T},
        normal::AbstractVector{T};
        cf::T = 1.0,
        p::AbstractVector{T} = szeros(T, 3),
        offset::AbstractVector{T} = szeros(T, 3),
        name::Symbol = Symbol("contact_" * randstring(4)),
        contact_type::Symbol = :contact) where {T}

    if contact_type == :contact
        bound = ContactBound(body, normal, cf, p=p, offset=offset)
    elseif contact_type == :linear_contact
        bound = LinearContactBound(body, normal, cf, p=p, offset=offset)
    elseif contact_type == :impact
        bound = ImpactBound(body, normal, p=p, offset=offset)
    else
        @warn "unknown contact_type"
    end
    ineqcs = InequalityConstraint((bound, body.id, nothing); name=name)
    return ineqcs
end
