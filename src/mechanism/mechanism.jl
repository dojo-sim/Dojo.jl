mutable struct Mechanism{T,Nn,Ne,Nb,Ni}
    origin::Origin{T}
    joints::Vector{<:JointConstraint{T}}
    bodies::Vector{<:Body{T}}
    contacts::Vector{<:ContactConstraint{T}}

    system::System{Nn}
    residual_entries::Vector{Entry}
    matrix_entries::SparseMatrixCSC{Entry,Int64}
    diagonal_inverses::Vector{Entry}

    timestep::T
    g::T
    μ::T
end

function Mechanism(origin::Origin{T}, bodies::Vector{<:Body{T}}, joints::Vector{<:JointConstraint{T}}, contacts::Vector{<:ContactConstraint{T}};
    spring=0.0, damper=0.0, timestep::T=0.01, g::T=-9.81) where T

    # reset ids
    resetGlobalID()

    # check body inertia parameters
    check_body.(bodies)

    # dimensions
    Ne = length(joints)
    Nb = length(bodies)
    Ni = length(contacts)
    Nn = Ne + Nb + Ni

    # nodes
    nodes = [joints; bodies; contacts]

    # ids
    oldnewid = Dict([node.id=>i for (i,node) in enumerate(nodes)]...)

    for node in nodes
        node.id = oldnewid[node.id]
        if typeof(node) <: Constraint
            node.parentid = get(oldnewid, node.parentid, nothing)
            node.childids = [get(oldnewid, childid, nothing) for childid in node.childids]
        end
    end

    # graph system
    system = create_system(origin, joints, bodies, contacts)
    residual_entries = deepcopy(system.vector_entries)
    matrix_entries = deepcopy(system.matrix_entries)
    diagonal_inverses = deepcopy(system.diagonal_inverses)

    # springs and dampers
    joints = set_spring_damper_values!(joints, spring, damper)

    Mechanism{T,Nn,Ne,Nb,Ni}(origin, joints, bodies, contacts, system, residual_entries, matrix_entries, diagonal_inverses, timestep, g, 0.0)
end

Mechanism(origin::Origin{T}, bodies::Vector{<:Body{T}}, joints::Vector{<:JointConstraint{T}}; kwargs...) where T = Mechanism(origin, bodies, joints, ContactConstraint{T}[]; kwargs...)

function Mechanism(filename::String, floating::Bool=false, T=Float64; kwargs...)
    # parse urdf
    origin, links, joints, loopjoints = parse_urdf(filename, floating, T)

    # create mechanism
    mechanism = Mechanism(origin, links, [joints; loopjoints]; kwargs...)

    # initialize mechanism
    set_parsed_values!(mechanism, loopjoints)

    return mechanism
end

Base.length(mechanism::Mechanism) =
    sum(Vector{Int}(length.(mechanism.joints))) +
    sum(Vector{Int}(length.(mechanism.bodies))) +
    sum(Vector{Int}(length.(mechanism.contacts)))
