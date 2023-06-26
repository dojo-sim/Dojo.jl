"""
    JointConstraint{T} <: Constraint{T}

    constraint restricting translational and rotational degrees of freedom between two Body objects.

    id: a unique identifying number
    name: a unique identifying name
    translational: Translational
    rotational: Rotational
    spring: flag for joint springs on
    damper: flag for joint dampers on
    parent_id: identifying number for parent Body{T}
    child_id: identifying number for child Body{T}
    minimal_index: indices for minimal coordinates
    impulses: joint impulses that maintain constraint between two Body{T} objects
"""
mutable struct JointConstraint{T,N,Nc,TJ,RJ} <: Constraint{T,N}
    # ID
    id::Int64
    name::Symbol

    # joint constraints
    translational::TJ
    rotational::RJ

    # springs and dampers
    spring::Bool
    damper::Bool

    # neighbor IDs
    parent_id::Int
    child_id::Int

    # indices
    minimal_index::SVector{Nc,SVector{2,Int64}} # indices for minimal coordinates, assumes joints # Nc = 2 THIS IS SPECIAL CASED

    # impulses
    impulses::Vector{SVector{N,T}}

    function JointConstraint(data;
        name::Symbol=Symbol("joint_" * randstring(4)))

        @assert data[1][2] == data[2][2] # check parent ids
        @assert data[1][3] == data[2][3] # check child ids

        # joints
        translational = data[1][1]
        rotational = data[2][1]

        # IDs
        parent_id = data[1][2]
        child_id = data[1][3]

        # data dype
        T = typeof(data[1][1]).parameters[1]

        # set springs & dampers off
        spring = false
        damper = false

        minimal_index = Vector{Int64}[]
        N = 0
        for joint_data in data
            joint = joint_data[1]

            # set spring & damper on
            joint.spring != 0 && (spring = true)
            joint.damper != 0 && (damper = true)

            # minimal-coordaintes indices
            Nλ = joint_length(joint)
            Nset = impulses_length(joint)
            if isempty(minimal_index)
                push!(minimal_index, [1;3-Nλ])
            else
                push!(minimal_index, [last(minimal_index)[2]+1; last(minimal_index)[2]+3-Nλ])
            end
            N += Nset
        end

        Nc = 2
        impulses = [zeros(T, N) for i=1:2]

        return new{T,N,Nc,typeof(translational),typeof(rotational)}(getGlobalID(), name, translational, rotational, spring, damper, parent_id, child_id, minimal_index, impulses)
    end
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, constraint::JointConstraint)
    summary(io, constraint)
    println(io, "")
    println(io, "id:            "*string(constraint.id))
    println(io, "name:          "*string(constraint.name))
    println(io, "spring:        "*string(constraint.spring))
    println(io, "damper:        "*string(constraint.damper))
    println(io, "parent_id:     "*string(constraint.parent_id))
    println(io, "child_id:      "*string(constraint.child_id))
    println(io, "minimal_index: "*string(constraint.minimal_index))
    println(io, "impulses:      "*string(constraint.impulses))
end
