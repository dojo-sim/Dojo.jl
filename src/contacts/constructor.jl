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
end

function ContactConstraint(data::Tuple; 
    name::Symbol=Symbol("contact_" * randstring(4)))

    model, parent_id, child_id = data
    T = typeof(model).parameters[1]

    N = length(model)
    N½ = Int64(N/2)

    impulses = [neutral_vector(model) for i = 1:2]
    impulses_dual = [neutral_vector(model) for i = 1:2]
    ContactConstraint{T,N,1,typeof(model),N½}(getGlobalID(), name, model, parent_id, child_id, impulses, impulses_dual)
end

function ContactConstraint(data::AbstractVector{<:Tuple}; 
    names::AbstractVector{Symbol}=[Symbol("contact_" * randstring(4)) for i = 1:length(data)])

    [ContactConstraint(data[i]; name=names[i]) for i in eachindex(data)]
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

function ContactConstraint(contact_type::Symbol, body::Body{T}, normal::AbstractVector;
    friction_coefficient=T(1),
    contact_origin::AbstractVector=szeros(T, 3),
    contact_radius=T(0),
    name::Symbol=Symbol("contact_" * randstring(4))) where T

    if contact_type == :nonlinear
        data = NonlinearContact(body, normal, friction_coefficient; 
            contact_origin, contact_radius)
    elseif contact_type == :linear
        data = LinearContact(body, normal, friction_coefficient; 
            contact_origin, contact_radius)
    elseif contact_type == :impact
        data = ImpactContact(body, normal; 
            contact_origin, contact_radius)
    else
        @warn "unknown contact_type"
    end

    ContactConstraint(data; name)
end


function ContactConstraint(contact_type::Symbol, body::Body{T}, normals::AbstractVector{<:AbstractVector};
    friction_coefficients::AbstractVector=ones(T,length(normals)),
    contact_origins::AbstractVector=[szeros(T, 3) for i=1:length(normals)],
    contact_radii::AbstractVector=[0.0 for i=1:length(normals)],
    names::Vector{Symbol}=[Symbol("contact_" * randstring(4)) for i = 1:length(normals)]) where T

    @assert length(normals) == length(friction_coefficients) == length(contact_origins) == length(contact_radii) == length(names)
    [ContactConstraint(contact_type, body, normals[i]; friction_coefficient=friction_coefficients[i], contact_origin=contact_origins[i], contact_radius=contact_radii[i], name=names[i]) for i in eachindex(normals)]
end

"""
    contact_constraint(bodies::Vector{Body}) 

    generate ContactConstraints for each Body in list

    normals: surface normal for each contact point
    friction_coefficients: value of coefficient of friction for each contact point (optional for ImpactContact)
    contact_origins: the offset with respect to the center of Body for each contact point (optional)
    contact_radius: radius for each contact (optional)
    contact_type: :nonlinear, :linear, :impact
"""
function ContactConstraint(contact_type::Symbol, bodies::Vector{Body{T}}, normals::AbstractVector{<:AbstractVector};
    friction_coefficients::AbstractVector=ones(T,length(normals)),
    contact_origins::AbstractVector=[szeros(T, 3) for i=1:length(normals)],
    contact_radii::AbstractVector=[0.0 for i=1:length(normals)],
    names::Vector{Symbol}=[Symbol("contact_" * randstring(4)) for i = 1:length(normals)]) where T

    @assert length(bodies) == length(normals) == length(friction_coefficients) == length(contact_origins) == length(contact_radii) == length(names)
    [ContactConstraint(contact_type, bodies[i], normals[i]; friction_coefficient=friction_coefficients[i], contact_origin=contact_origins[i], contact_radius=contact_radii[i], name=names[i]) for i in eachindex(bodies)]
end