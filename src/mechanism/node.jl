abstract type Node{T} end
abstract type Constraint{T,N} <: Node{T} end

# TODO do id differently?
CURRENTID = -1
getGlobalID() = (global CURRENTID -= 1; return CURRENTID + 1)
resetGlobalID() = (global CURRENTID = -1; return)

Base.length(::Constraint{T,N}) where {T,N} = N
