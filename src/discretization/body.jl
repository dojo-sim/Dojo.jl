@inline function g(mechanism::Mechanism{T,Nn,Ne,Nb}, body::Body{T}) where {T,Nn,Ne,Nb}
    state = body.state
    Δt = mechanism.Δt

    x1, q1 = posargs1(state)
    x2, q2 = posargs2(state)
    x3, q3 = posargs3(state, Δt)

    # dynamics
    ezg = SA{T}[0; 0; -mechanism.g]
    D1x = - 1 / Δt * body.m * (x2 - x1) + Δt/2 * body.m * ezg
    D2x =   1 / Δt * body.m * (x3 - x2) + Δt/2 * body.m * ezg
    D1q = -4 / Δt * LVᵀmat(q2)' * Lmat(q1) * Vᵀmat() * body.J * Vmat() * Lmat(q1)' * vector(q2)
    D2q = -4 / Δt * LVᵀmat(q2)' * Tmat() * Rmat(q3)' * Vᵀmat() * body.J * Vmat() * Lmat(q2)' * vector(q3)

    dynT = D2x + D1x
    dynR = D2q + D1q

    state.d = [dynT; dynR]

    # inputs
    state.d -= [state.F2[1]; 2.0 * state.τ2[1]]

    # impulses
    for id in connections(mechanism.system, body.id)
        Ne < id <= Ne+Nb && continue # body
        impulses!(mechanism, body, getcomponent(mechanism, id))
    end

    # Regularize the angular velocity when necessary.
    # for body in mechanism.bodies
    #     angular_damping!(mechanism, body)
    # end

    return state.d
end

@inline function ∂g∂z(mechanism::Mechanism{T,Nn,Ne,Nb}, body::Body{T}) where {T,Nn,Ne,Nb}
    state = body.state
    Δt = mechanism.Δt
    x1, q1 = posargs1(state)
    x2, q2 = posargs2(state)
    x3, q3 = posargs3(state, Δt)

    # dynamics
    dynT = I(3) * body.m / Δt

    rot_q3(q) = -4 / Δt * LVᵀmat(q2)' * Tmat() * Rmat(UnitQuaternion(q..., false))' * Vᵀmat() * body.J * Vmat() * Lmat(q2)' * q
    dynR = FiniteDiff.finite_difference_jacobian(rot_q3, vector(q3)) #* ∂integrator∂ϕ(q2, state.ϕsol[2], Δt)

    Z33 = szeros(T, 3, 3)
    Z34 = szeros(T, 3, 4)

    state.D = [[dynT; Z33] [Z34; dynR]] * ∂i∂v(body, mechanism.Δt)

    # inputs
    nothing

    # impulses
    for id in connections(mechanism.system, body.id)
        Ne < id <= Ne+Nb && continue # body
        ∂impulses∂v!(mechanism, body, getcomponent(mechanism, id))
    end

    # regularize the angular velocity when necessary.
    # for body in mechanism.bodies
    #     ∂angular_damping!(mechanism, body)
    # end

    return state.D
end

@inline function ∂i∂v(body::Body{T}, Δt) where {T}
    state = body.state
    x2, v25, q2, ϕ25 = fullargssol(state)
    ∂i∂v(q2, ϕ25, Δt)
end

@inline function ∂i∂z(body::Body{T}, Δt; attjac::Bool=true) where {T}
    state = body.state
    x2, v25, q2, ϕ25 = fullargssol(state)
    ∂i∂z(q2, ϕ25, Δt, attjac=attjac)
end
