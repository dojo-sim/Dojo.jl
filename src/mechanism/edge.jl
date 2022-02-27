"""
    Constraint{T}

    Abstract type for graph edge object
"""
abstract type Constraint{T,N} <: Node{T} end
Base.length(::Constraint{T,N}) where {T,N} = N