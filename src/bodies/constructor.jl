"""
    Body{T} <: Node{T}

    A rigid body object

    id: unique identifying number 
    name: unique identifying name 
    mass: inertial property (kilograms)
    inertia: inertia matrix (kilograms meter^2)
    state: State; representation of the system's: position, linear velocity, orientation, and angular velocity 
    shape: Shape; geometry information about the Body
"""
mutable struct Body{T} <: Node{T}
    id::Int64
    name::Symbol
    mass::T
    inertia::SMatrix{3,3,T,9}
    state::State{T}
    shape::Shape{T}

    function Body(mass::T, inertia::AbstractMatrix; 
        name::Symbol=:origin, 
        shape::Shape=EmptyShape()) where T
        new{T}(getGlobalID(), name, mass, inertia, State{T}(), shape)
    end
end

Base.length(::Body) = 6

# warn that body has poor inertial properties
function check_body(body::Body)
    if norm(body.mass) == 0 || norm(body.inertia) == 0
        @info "Bad inertial properties detected"
    end
    #TODO check condition number and potentially adapt
end


