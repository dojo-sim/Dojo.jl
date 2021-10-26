@inline function g(mechanism::Mechanism{T,Nn,Ne,Nb}, body::Body{T}) where {T,Nn,Ne,Nb}
    state = body.state
    Δt = mechanism.Δt

    # ezg = SA{T}[0; 0; -mechanism.g]
    ezg = SA{T}[0; 0; -mechanism.g * Δt]
    # dynT = body.m * ((state.vsol[2] - state.vc) / Δt + ezg) - state.Fk[1]
    dynT = body.m * ((state.vsol[2] - state.vc) + ezg) - state.Fk[1]

    J = body.J
    ω1 = state.ωc
    ω2 = state.ωsol[2]
    sq1 = sqrt(4 / Δt^2 - ω1' * ω1)
    sq2 = sqrt(4 / Δt^2 - ω2' * ω2)
    # dynR = skewplusdiag(ω2, sq2) * (J * ω2) - skewplusdiag(ω1, sq1) * (J * ω1) - 2 * state.τk[1]
    dynR = Δt * skewplusdiag(ω2, sq2) * (J * ω2) - Δt * skewplusdiag(ω1, sq1) * (J * ω1) - 2 * state.τk[1]

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
    ω2 = state.ωsol[2]
    sq = sqrt(4 / Δt^2 - ω2' * ω2)

    # dynT = I * body.m / Δt
    dynT = I * body.m
    dynR = Δt * (skewplusdiag(ω2, sq) * J - J * ω2 * (ω2' / sq) - skew(J * ω2))

    Z = szeros(T, 3, 3)

    state.D = [[dynT; Z] [Z; dynR]]

    for connectionid in connections(mechanism.system, body.id)
        Ne < connectionid <= Ne+Nb && continue # body
        ∂constraintForceMapping!(mechanism, body, getcomponent(mechanism, connectionid))
    end
    

    return state.D
end
