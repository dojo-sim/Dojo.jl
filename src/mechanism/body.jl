mutable struct Body{T} <: Node{T}
    id::Int64
    name::Symbol
    m::T
    J::SMatrix{3,3,T,9}
    state::State{T}
    shape::Shape{T}

    function Body(m::T, J::AbstractMatrix; name::Symbol=:origin, shape::Shape=EmptyShape()) where T
        new{T}(getGlobalID(), name, m, J, State{T}(), shape)
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



