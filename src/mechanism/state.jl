"""
$(TYPEDEF)

The `State` contains the position and velocity information of a [`Body`](@ref).
# Important attributes
* `x1`: Continuous-time position (vector)
* `q1`: Continuous-time oritentation (quaternion)
* `v15`: Continuous-time velocity (vector)
* `ϕ15`: Continuous-time angular velocity (vector)
* `x2`: Knot-point positions (vector of vectors)
* `q2`: Knot-point orientations (vector of vectors)
"""
mutable struct State{T}
    order::Integer

    # Continuous states
    x1::SVector{3,T}
    q1::UnitQuaternion{T}
    v15::SVector{3,T}
    ϕ15::SVector{3,T}

    # Knot points
    x2::Vector{SVector{3,T}}
    q2::Vector{UnitQuaternion{T}}
    F2::Vector{SVector{3,T}}
    τ2::Vector{SVector{3,T}}

    # Current solution estimate [before step;after step]
    vsol::Vector{SVector{3,T}}
    ϕsol::Vector{SVector{3,T}}

    # Current equations of motion and derivative evaluation
    d::SVector{6,T}
    D::SMatrix{6,6,T,36}

    function State{T}() where T
        x1 = szeros(T, 3)
        q1 = one(UnitQuaternion{T})
        v15 = szeros(T, 3)
        ϕ15 = szeros(T, 3)

        x2 = [szeros(T, 3)]
        q2 = [one(UnitQuaternion{T})]
        F2 = [szeros(T, 3)]
        τ2 = [szeros(T, 3)]

        vsol = [szeros(T, 3) for i=1:2]
        ϕsol = [szeros(T, 3) for i=1:2]

        d = szeros(T, 6)
        D = szeros(T, 6, 6)

        new{T}(0, x1, q1, v15, ϕ15, x2, q2, F2, τ2, vsol, ϕsol, d, D)
    end
end

# function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, state::State{T}) where {T}
#     summary(io, state)
#     println(io,"")
#     println(io,"x1:   "*string(state.x1))
#     println(io,"q1:   "*string(state.q1))
#     println(io,"ϕ15:   "*string(state.ϕ15))
#     println(io,"x2:   "*string(state.x2))
#     println(io,"q2:   "*string(state.q2))
#     println(io,"F2:   "*string(state.F2))
#     println(io,"τ2:   "*string(state.τ2))
#     println(io,"vsol: "*string(state.vsol))
#     println(io,"ϕsol: "*string(state.ϕsol))
#     println(io,"d:    "*string(state.d))
# end


function initknotpoints!(state::State, order)
    state.order = order

    state.x2 = [state.x2[1] for i = 1:order]
    state.q2 = [state.q2[1] for i = 1:order]
    state.F2 = [state.F2[1] for i = 1:order]
    state.τ2 = [state.τ2[1] for i = 1:order]

    return
end
