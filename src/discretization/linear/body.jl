# @inline function g(mechanism::Mechanism{T,Nn,Ne,Nb}, body::Body{T}) where {T,Nn,Ne,Nb}
#     state = body.state
#     Δt = mechanism.Δt

#     ezg = SA{T}[0; 0; -mechanism.g * Δt]
#     dynT = body.m * ((state.vsol[2] - state.v15) + ezg) - state.F2[1]

#     J = body.J
#     ω1 = state.ϕ15
#     ω2 = state.ϕsol[2]
#     sq1 = sqrt(4 / Δt^2 - ω1' * ω1)
#     sq2 = sqrt(4 / Δt^2 - ω2' * ω2)
#     dynR = Δt * skewplusdiag(ω2, sq2) * (J * ω2) - Δt * skewplusdiag(-1.0 * ω1, sq1) * (J * ω1) - 2 * state.τ2[1]

#     state.d = [dynT;dynR]

#     for connectionid in connections(mechanism.system, body.id)
#         Ne < connectionid <= Ne+Nb && continue # body
#         constraintForceMapping!(mechanism, body, getcomponent(mechanism, connectionid))
#     end

#     return state.d
# end

# @inline function ∂g∂ʳself(mechanism::Mechanism{T,Nn,Ne,Nb}, body::Body{T}) where {T,Nn,Ne,Nb}
#     state = body.state
#     Δt = mechanism.Δt
#     J = body.J
#     ω2 = state.ϕsol[2]
#     sq = sqrt(4 / Δt^2 - ω2' * ω2)

#     dynT = I * body.m
#     dynR = Δt * (skewplusdiag(ω2, sq) * J - J * ω2 * (ω2' / sq) - skew(J * ω2))

#     Z = szeros(T, 3, 3)

#     state.D = [[dynT; Z] [Z; dynR]]

#     for connectionid in connections(mechanism.system, body.id)
#         Ne < connectionid <= Ne+Nb && continue # body
#         ∂constraintForceMapping!(mechanism, body, getcomponent(mechanism, connectionid))
#     end

#     return state.D
# end

@inline function g(mechanism::Mechanism{T,Nn,Ne,Nb}, body::Body{T}) where {T,Nn,Ne,Nb}
    state = body.state
    Δt = mechanism.Δt

    x1, q1 = posargs1(state)
    x2, q2 = posargs2(state)
    x3, q3 = posargs3(state, Δt)
    
    ezg = SA{T}[0; 0; -mechanism.g]
    D1x = - 1/Δt * body.m * (x2 - x1) + Δt/2 * body.m * ezg
    D2x =   1/Δt * body.m * (x3 - x2) + Δt/2 * body.m * ezg
    D1q =   2/Δt * LVᵀmat(q1)' * Tmat() * Rmat(q2)' * Vᵀmat() * body.J * Vmat() * Lmat(q1)' * vector(q2)
    D2q =   2/Δt * LVᵀmat(q3)' * Lmat(q2) * Vᵀmat() * body.J * Vmat() * Lmat(q2)' * vector(q3)

    dynT = D2x + D1x - state.F2[1]
    dynR = D2q + D1q - state.τ2[1]

    state.d = [dynT;dynR]

    for connectionid in connections(mechanism.system, body.id)
        Ne < connectionid <= Ne+Nb && continue # body
        constraintForceMapping!(mechanism, body, getcomponent(mechanism, connectionid))
    end

    return state.d
end

@inline function ∂g∂ʳself(mechanism::Mechanism{T,Nn,Ne,Nb}, body::Body{T}) where {T,Nn,Ne,Nb}
    state = body.state
    Δt = mechanism.Δt
    J = body.J
    x1, q1 = posargs1(state)
    x2, q2 = posargs2(state)
    x3, q3 = posargs3(state, Δt)

    dynT = I * body.m

    rot_q3(q) = 2/Δt * LVᵀmat(UnitQuaternion(q..., false))' * Lmat(q2) * Vᵀmat() * body.J * Vmat() * Lmat(q2)' * vector(UnitQuaternion(q..., false))
    dynR = FiniteDiff.finite_difference_jacobian(rot_q3, vector(q3)) * ∂integrator∂ϕ(q2, state.ϕsol[2], Δt)

    Z = szeros(T, 3, 3)

    state.D = [[dynT; Z] [Z; dynR]]

    for connectionid in connections(mechanism.system, body.id)
        Ne < connectionid <= Ne+Nb && continue # body
        ∂constraintForceMapping!(mechanism, body, getcomponent(mechanism, connectionid))
    end

    return state.D
end