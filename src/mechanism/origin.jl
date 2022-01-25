mutable struct Origin{T} <: Component{T}
    id::Int64
    name::Symbol
    state::State{T}
    shape::Shape{T}

    function Origin{T}(; name::Symbol=:origin, state=State{T}(), shape::Shape=EmptyShape()) where T
        new{T}(getGlobalID(), name, state, shape)
    end

    function Origin(; name::Symbol=:origin, state=State{Float64}(), shape::Shape=EmptyShape())
        Origin{Float64}(; name=name, state=state, shape=shape)
    end

    function Origin(body::Body{T}) where T
        new{T}(body.id, body.name, body.state, body.shape)
    end
end
