mutable struct Origin{T} <: Component{T}
    id::Int64
    name::String
    shape::Shape{T}

    function Origin{T}(; name::String="", shape::Shape=EmptyShape()) where T
        new{T}(getGlobalID(), name, shape)
    end

    function Origin(; name::String="", shape::Shape=EmptyShape())
        Origin{Float64}(; name = name, shape = shape)
    end

    function Origin(body::Body{T}) where T
        new{T}(body.id, body.name, body.shape)
    end
end