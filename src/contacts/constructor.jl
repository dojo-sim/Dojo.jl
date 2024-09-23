"""
    ContactConstraint{T} <: Constraint{T}

    constraint containing information for contact node.

    id: unique identifying number 
    name: unique identifying name 
    model: type of contact model: ImpactContact, LinearContact, NonlinearContact 
    parent_id: identifying number of Body experiencing contact 
    child_id: always 0
    impulses: contact impulses applied to Body 
    impulses_dual: dual contact impulses, used by solver to enforce correct contact behaviors 
"""
mutable struct ContactConstraint{T,N,Nc,Cs,N½} <: Constraint{T,N}
    # ID
    id::Int64
    name::Symbol

    # contact model
    model::Cs

    # neighbor IDs
    parent_id::Int
    child_id::Int

    # variables
    impulses::Vector{SVector{N½,T}}
    impulses_dual::Vector{SVector{N½,T}}

    function ContactConstraint(data; 
        name::Symbol=Symbol("contact_" * randstring(4)))

        model, parent_id, child_id = data
        T = typeof(model).parameters[1]

        N = length(model)
        N½ = Int64(N/2)

        impulses = [neutral_vector(model) for i = 1:2]
        impulses_dual = [neutral_vector(model) for i = 1:2]
        new{T,N,1,typeof(model),N½}(getGlobalID(), name, model, parent_id, child_id, impulses, impulses_dual)
    end
end

# function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, constraint::ContactConstraint)
#     summary(io, constraint)
#     println(io, "")
#     println(io, "id:            "*string(constraint.id))
#     println(io, "name:          "*string(constraint.name))
#     println(io, "model:         "*string(constraint.model))
#     println(io, "parent_id:     "*string(constraint.parent_id))
#     println(io, "child_id:      "*string(constraint.child_id))
#     println(io, "impulses:      "*string(constraint.impulses))
#     println(io, "impulses_dual: "*string(constraint.impulses_dual))
# end

"""
    contact_constraint(bodies::Vector{Body}) 

    generate ContactConstraints for each Body in list

    normals: surface normal for each contact point
    friction_coefficients: value of coefficient of friction for each contact point (optional for ImpactContact)
    contact_origins: the offset with respect to the center of Body for each contact point (optional)
    contact_radius: radius for each contact (optional)
    contact_type: :nonlinear, :linear, :impact
"""
function contact_constraint(bodies::Vector{Body{T}},
        normals::AbstractVector{<:AbstractVector};
        friction_coefficients::AbstractVector=ones(T,length(normals)),
        contact_origins::AbstractVector=[szeros(T, 3) for i=1:length(normals)],
        contact_radii::AbstractVector=[0.0 for i=1:length(normals)],
        contact_offsets::AbstractVector=[szeros(T, 3)  for i=1:length(normals)],
        names::Vector{Symbol}=[Symbol("contact_" * randstring(4)) for i = 1:length(normals)],
        contact_type::Symbol=:nonlinear) where T

    n = length(normals)
    @assert n == length(bodies) == length(normals) == length(friction_coefficients) == length(contact_origins) == length(contact_radii) == length(contact_offsets)
    contacts = Vector{ContactConstraint}()
    for i = 1:n
        contact = contact_constraint(bodies[i], normals[i], 
            friction_coefficient=friction_coefficients[i], 
            contact_origin=contact_origins[i],
            contact_radius=contact_radii[i], 
            contact_offset=contact_offsets[i],
            name=names[i], 
            contact_type=contact_type)
        push!(contacts, contact)
    end
    contacts = [contacts...] # vector typing
    return contacts
end

function contact_constraint(body::Body{T},
        normals::AbstractVector{<:AbstractVector};
        friction_coefficients::AbstractVector=ones(T,length(normals)),
        contact_origins::AbstractVector=[szeros(T, 3) for i=1:length(normals)],
        contact_radii::AbstractVector=[0.0 for i=1:length(normals)],
        contact_offsets::AbstractVector=[szeros(T, 3) for i=1:length(normals)],
        names::Vector{Symbol}=[Symbol("contact_" * randstring(4)) for i = 1:length(normals)],
        contact_type::Symbol=:nonlinear) where T
    n = length(normals)
    @assert n == length(normals) == length(friction_coefficients) == length(contact_origins) == length(contact_radii)
    return contact_constraint(fill(body, n), normals;
        friction_coefficients, contact_origins, contact_radii, contact_offsets, names, contact_type)
end

function contact_constraint(body::Body{T},
        normal::AbstractVector;
        friction_coefficient=T(1),
        contact_origin::AbstractVector=szeros(T, 3),
        contact_radius=T(0),
        contact_offset::AbstractVector=szeros(T, 3),
        name::Symbol=Symbol("contact_" * randstring(4)),
        contact_type::Symbol=:nonlinear) where T

    if contact_type == :nonlinear
        model = NonlinearContact(body, normal, friction_coefficient; 
            contact_origin, contact_radius, contact_offset)
    elseif contact_type == :linear
        model = LinearContact(body, normal, friction_coefficient; 
            contact_origin, contact_radius, contact_offset)
    elseif contact_type == :impact
        model = ImpactContact(body, normal; 
            contact_origin, contact_radius, contact_offset)
    else
        @warn "unknown contact_type"
    end
    contacts = ContactConstraint((model, body.id, 0); name)
    return contacts
end
