abstract type Node{T} end

abstract type Constraint{T,N} <: Node{T} end
Base.length(::Constraint{T,N}) where {T,N} = N