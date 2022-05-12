"""
    Origin{T} <: Node{T}

    Global reference frame

    id: always 0 
    name: always :origin 
    state: State; defaults to zero values
    shape: Shape; defaults to EmptyShape
"""
mutable struct Origin{T} <: Node{T}
    id::Int64
    name::Symbol
    state::State{T}
    shape::Shape{T}

    function Origin{T}(; 
        name::Symbol=:origin, 
        state=State{T}(), 
        shape::Shape=EmptyShape()) where T
        new{T}(getGlobalID(), name, state, shape)
    end

    function Origin(body::Body{T}) where T
        new{T}(body.id, body.name, body.state, body.shape)
    end
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, origin::Origin)
    summary(io, origin)
    println(io, "")
    println(io, " id:      "*string(origin.id))
    println(io, " name:    "*string(origin.name))
end