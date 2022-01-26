mutable struct State{T}
    # previous state
    x1::SVector{3,T}
    q1::UnitQuaternion{T}
    v15::SVector{3,T}
    ϕ15::SVector{3,T}

    # current state
    x2::Vector{SVector{3,T}}
    q2::Vector{UnitQuaternion{T}}
    F2::Vector{SVector{3,T}}
    τ2::Vector{SVector{3,T}}

    # solution estimate [before step; after step]
    vsol::Vector{SVector{3,T}}
    ϕsol::Vector{SVector{3,T}}

    # dynamics and Jacobian evaluation
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

        new{T}(x1, q1, v15, ϕ15, x2, q2, F2, τ2, vsol, ϕsol, d, D)
    end
end

function initialize!(state::State, order=1)
    state.x2 = [state.x2[1] for i = 1:order]
    state.q2 = [state.q2[1] for i = 1:order]
    state.F2 = [state.F2[1] for i = 1:order]
    state.τ2 = [state.τ2[1] for i = 1:order]
end
