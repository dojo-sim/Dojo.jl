mutable struct Origin{T} <: Component{T}
    id::Int64
    name::String
    state::State{T}
    shape::Shape{T}

    function Origin{T}(; name::String="", state=State{T}(), shape::Shape=EmptyShape()) where T
        new{T}(getGlobalID(), name, state, shape)
    end

    function Origin(; name::String="", shape::Shape=EmptyShape())
        Origin{Float64}(; name=name, state=state, shape=shape)
    end

    function Origin(body::Body{T}) where T
        new{T}(body.id, body.name, body.state, body.shape)
    end
end