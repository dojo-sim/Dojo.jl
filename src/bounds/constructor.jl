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
        names::AbstractVector{String} = fill("", length(normal)),
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
        names::AbstractVector{String} = fill("", length(normal)),
        contact_type::Symbol = :contact) where {T}
    n = length(normal)
    @assert n == length(normal) == length(cf) == length(p) == length(offset)
    return contactconstraint(fill(body, n), normal, cf, p=p, offset=offset,
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
        name::String = "",
        contact_type::Symbol = :contact) where {T}

    if contact_type == :contact
        bound = ContactBound(body, normal, cf=cf, p=p, offset=offset)
    elseif contact_type == :linear_contact
        bound = LinearContactBound(body, normal, cf=cf, p=p, offset=offset)
    elseif contact_type == :impact
        bound = ImpactBound(body, normal, p=p, offset=offset)
    else
        @warn "unknown contact_type"
    end
    ineqcs = InequalityConstraint((bound, body.id, nothing); name=name)
    return ineqcs
end



#
# #######################################################
# #######################################################
# #######################################################
# #######################################################
#
#
#
# """
#     Generate linear contact inequality constraints attached to a list of bodies. You need to provide:
#     - the normal for each contact point
#     - the coefficient of friction for each contact point
#     - the offset vector p with respect to the center of the body for each contact point (optional)
#     - the altitude offset for each contact point (optional)
# """
# function linearcontactconstraint(bodies::AbstractVector{<:Body{T}}, normal::AbstractVector{<:AbstractVector}, cf::AbstractVector{T};
#         p::AbstractVector = [szeros(T, 3) for i=1:length(normal)],
#         offset::AbstractVector = [szeros(T, 3) for i=1:length(normal)]) where {T}
#
#     n = length(normal)
#     @assert n == length(bodies) == length(normal) == length(cf) == length(p) == length(offset)
#     linineqcs = Vector{InequalityConstraint}()
#     for i = 1:n
#         linineqc = linearcontactconstraint(bodies[i], normal[i], cf[i], p = p[i], offset = offset[i])
#         push!(linineqcs, linineqc)
#     end
#     linineqcs = [linineqcs...] # vector typing
#     return linineqcs
# end
#
# function linearcontactconstraint(body::Body{T}, normal::AbstractVector{<:AbstractVector}, cf::AbstractVector{T};
#         p::AbstractVector = [szeros(T, 3) for i=1:length(normal)],
#         offset::AbstractVector = [szeros(T, 3) for i=1:length(normal)]) where {T}
#     n = length(normal)
#     @assert n == length(normal) == length(cf) == length(p) == length(offset)
#     return linearcontactconstraint(fill(body, n), normal, cf, p = p, offset = offset)
# end
#
# """
#     Generate linear contact inequality constraint attached to one body. You need to provide:
#     - the normal for the contact point
#     - the coefficient of friction for the contact point
#     - the offset vector p with respect to the center of the body for the contact point (optional)
#     - the altitude offset for the contact point (optional)
# """
# function linearcontactconstraint(body::Body{T}, normal::AbstractVector{T}, cf::T;
#         p::AbstractVector{T} = szeros(T, 3),
#         offset::AbstractVector{T} = szeros(T, 3)) where {T}
#
#     linbound = LinearContactBound(body, normal, cf, p = p, offset = offset)
#     linineqcs = InequalityConstraint((linbound, body.id, nothing))
#     return linineqcs
# end
#
#
#
#
# #######################################################
# #######################################################
# #######################################################
# #######################################################
#
#
# """
#     Generate impact inequality constraints attached to a list of bodies. You need to provide:
#     - the normal for each contact point
#     - the coefficient of friction for each contact point
#     - the offset vector p with respect to the center of the body for each contact point (optional)
#     - the altitude offset for each contact point (optional)
# """
# function impactconstraint(bodies::AbstractVector{<:Body{T}}, normal::AbstractVector{<:AbstractVector};
#         p::AbstractVector = [szeros(T, 3) for i=1:length(normal)],
#         offset::AbstractVector = [szeros(T, 3) for i=1:length(normal)]) where {T}
#
#     n = length(normal)
#     @assert n == length(bodies) == length(normal) == length(p) == length(offset)
#     impineqcs = Vector{InequalityConstraint}()
#     for i = 1:n
#         impineqc = impactconstraint(bodies[i], normal[i], p = p[i], offset = offset[i])
#         push!(impineqcs, impineqc)
#     end
#     impineqcs = [impineqcs...] # vector typing
#     return impineqcs
# end
#
# function impactconstraint(body::Body{T}, normal::AbstractVector{<:AbstractVector};
#         p::AbstractVector = [szeros(T, 3) for i=1:length(normal)],
#         offset::AbstractVector = [szeros(T, 3) for i=1:length(normal)]) where {T}
#     n = length(normal)
#     @assert n == length(normal) == length(p) == length(offset)
#     return impactconstraint(fill(body, n), normal, p = p, offset = offset)
# end
#
# """
#     Generate impact inequality constraint attached to one body. You need to provide:
#     - the normal for the contact point
#     - the offset vector p with respect to the center of the body for the contact point (optional)
#     - the altitude offset for the contact point (optional)
# """
# function impactconstraint(body::Body{T}, normal::AbstractVector{T};
#         p::AbstractVector{T} = szeros(T, 3),
#         offset::AbstractVector{T} = szeros(T, 3)) where {T}
#
#     impbound = ImpactBound(body, normal, p = p, offset = offset)
#     impineqcs = InequalityConstraint((impbound, body.id, nothing))
#     return impineqcs
# end
