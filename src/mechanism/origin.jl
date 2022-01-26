mutable struct Origin{T} <: Node{T}
    id::Int64
    name::Symbol
    state::State{T}
    shape::Shape{T}

    function Origin{T}(; name::Symbol=:origin, state=State{T}(), shape::Shape=EmptyShape()) where T
        new{T}(getGlobalID(), name, state, shape)
    end

    function Origin(body::Body{T}) where T
        new{T}(body.id, body.name, body.state, body.shape)
    end
end
