abstract type Component{T} end
abstract type AbstractBody{T} <: Component{T} end
abstract type AbstractConstraint{T,N} <: Component{T} end


# TODO do id differently?
CURRENTID = -1
getGlobalID() = (global CURRENTID -= 1; return CURRENTID + 1)
resetGlobalID() = (global CURRENTID = -1; return)

Base.eltype(::Type{<:Component{E}}) where {E} = @isdefined(E) ? E : Any
Base.length(::AbstractConstraint{T,N}) where {T,N} = N
getid(component::Component) = component.id

function Base.getindex(dict::UnitDict{Base.OneTo{K},<:Component}, key::String) where K
    for component in dict.values
        component.name == key && return component
    end
    
    return
end
function Base.getindex(dict::UnitDict{UnitRange{K},<:Component}, key::String) where K
    for component in dict.values
        component.name == key && return component
    end
    
    return
end