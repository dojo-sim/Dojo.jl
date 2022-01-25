abstract type Component{T} end
abstract type Constraint{T,N} <: Component{T} end

# TODO do id differently?
CURRENTID = -1
getGlobalID() = (global CURRENTID -= 1; return CURRENTID + 1)
resetGlobalID() = (global CURRENTID = -1; return)

Base.eltype(::Type{<:Component{E}}) where {E} = @isdefined(E) ? E : Any
Base.length(::Constraint{T,N}) where {T,N} = N
getid(component::Component) = component.id