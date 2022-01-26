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
end

Base.length(::Body) = 6

function check_body(body::Body)
    initialize!(body.state)
    if norm(body.m) == 0 || norm(body.J) == 0
        @info "Bad inertial properties detected"
    end
end



