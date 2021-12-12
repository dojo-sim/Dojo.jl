mutable struct Body{T} <: Component{T}
    id::Int64
    name::String
    m::T
    J::SMatrix{3,3,T,9}
    state::State{T}
    shape::Shape{T}
    parentid::Int64
    childid::Int64

    function Body(m::T, J::SMatrix{3,3,T,9}; name::String="", shape::Shape=EmptyShape()) where T
        new{T}(getGlobalID(), name, m, J, State{T}(), shape, 0, 0)
    end

    function Body{T}(contents...) where T
        new{T}(getGlobalID(), contents...)
    end
end

function Base.deepcopy(b::Body{T}) where T
    contents = []
    for i = 2:getfieldnumber(b)
        push!(contents, deepcopy(getfield(b, i)))
    end

    return Body{T}(contents...)
end

Base.length(::Body) = 6
Base.zero(::Body{T}) where T = szeros(T, 6, 6)

