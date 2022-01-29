mutable struct Body{T} <: Node{T}
    id::Int64
    name::Symbol
    mass::T
    inertia::SMatrix{3,3,T,9}
    state::State{T}
    shape::Shape{T}

    function Body(mass::T, inertia::AbstractMatrix; name::Symbol=:origin, shape::Shape=EmptyShape()) where T
        new{T}(getGlobalID(), name, mass, inertia, State{T}(), shape)
    end
end

Base.length(::Body) = 6

function check_body(body::Body)
    initialize!(body.state)
    if norm(body.mass) == 0 || norm(body.inertia) == 0
        @info "Bad inertial properties detected"
    end
end



