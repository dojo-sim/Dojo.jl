@inline function g(mechanism::Mechanism{T,Nn,Ne,Nb}, body::Body{T}) where {T,Nn,Ne,Nb}
    state = body.state
    Δt = mechanism.Δt

    ezg = SA{T}[0; 0; -mechanism.g * Δt]
    dynT = body.m * ((state.vsol[2] - state.vc) + ezg) - state.Fk[1]

    J = body.J
    ω1 = state.ωc
    ω2 = state.ωsol[2]
    dynR = sqrt(1.0 - ω2' * ω2) * J * ω2 + cross(ω2, J * ω2) - (sqrt(1.0 - ω1' * ω1) * J * ω1 + cross(ω1, J * ω1)) + 0.5 * Δt^2.0 * state.τk[1]

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

    dynT = I * body.m
    dynR = 2.0 / Δt * (sqrt(1.0 - ω2' * ω2) * J - (J * ω2 * ω2') / sqrt(1.0 - ω2' * ω2) + skew(ω2) * J - skew(J * ω2))

    Z = szeros(T, 3, 3)

    state.D = [[dynT; Z] [Z; dynR]]

    for connectionid in connections(mechanism.system, body.id)
        Ne < connectionid <= Ne+Nb && continue # body
        ∂constraintForceMapping!(mechanism, body, getcomponent(mechanism, connectionid))
    end

    return state.D
end

# Random.seed!(100)
# Δt = 0.01
# ω2 = rand(3)
# ω1 = rand(3)
# q2 = UnitQuaternion(rand(4)...)
# q1 = Δt/2 * Lmat(q2) * [sqrt(4/Δt^2 - ω1'*ω1); -ω1]
# q3 = Δt/2 * Lmat(q2) * [sqrt(4/Δt^2 - ω2'*ω2); ω2]
#
# q21_ = Δt/2 * Lmat(q1) * [sqrt(4/Δt^2 - ω1'*ω1); ω1]
# q23_ = Δt/2 * Lmat(q3) * [sqrt(4/Δt^2 - ω2'*ω2); -ω2]
#
# q2_ = [q2.w, q2.x, q2.y, q2.z]
# norm(q2_ - q21_)
# norm(q2_ - q23_)
#
#
#
# sq2 = sqrt(4 / Δt^2 - ω2' * ω2)
# D2 = Δt * skewplusdiag(ω2, sq2) * (body2.J * ω2)
# D2_ = body2.J * ω2 * sq2 + skew(ω2) * body2.J * ω2
# norm(D2_ * Δt - D2)
#
#
# Vmat(UnitQuaternion(q2_, false) \ UnitQuaternion(q3..., false))
# ω2 * Δt_
#
# mech.eqconstraints[2].constraints[2].qoffset
