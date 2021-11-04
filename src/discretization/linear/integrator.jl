METHODORDER = 1 # This refers to the interpolating spline
getGlobalOrder() = (global METHODORDER; return METHODORDER)

# Convenience functions
@inline getx3(state::State, Δt) = state.xk[1] + state.vsol[2]*Δt
function getq3(state::State, Δt) 
    ϕ = state.ωsol[2]
    f = [ϕ; sqrt(1.0 - ϕ' * ϕ)]
    state.qk[1] * UnitQuaternion(f...) 
end

@inline posargsc(state::State) = (state.xc, state.qc)
@inline fullargsc(state::State) = (state.xc, state.vc, state.qc, state.ωc)
@inline posargsk(state::State; k=1) = (state.xk[k], state.qk[k])
@inline posargssol(state::State) = (state.xsol[2], state.qsol[2])
@inline fullargssol(state::State) = (state.xsol[2], state.vsol[2], state.qsol[2], state.ωsol[2])
@inline posargsnext(state::State, Δt) = (getx3(state, Δt), getq3(state, Δt))

@inline function derivωbar(ω::SVector{3}, Δt)
    msq = -sqrt(4 / Δt^2 - dot(ω, ω))
    return [ω' / msq; I]
end

@inline function ωbar(ω, Δt)
    return UnitQuaternion(sqrt(4 / Δt^2 - dot(ω, ω)), ω, false)
end

@inline function setForce!(state::State, F, τ)
    state.Fk[1] = F
    state.τk[1] = τ
    return
end

@inline function discretizestate!(body::Body{T}, Δt) where T
    state = body.state
    xc = state.xc
    qc = state.qc
    vc = state.vc
    ωc = state.ωc

    state.xk[1] = xc + vc*Δt
    state.qk[1] = qc * ωbar(ωc,Δt) * Δt / 2

    state.Fk[1] = szeros(T,3)
    state.τk[1] = szeros(T,3)

    return
end

@inline function currentasknot!(body::Body)
    state = body.state

    state.xk[1] = state.xc
    state.qk[1] = state.qc

    return
end

@inline function updatestate!(body::Body{T}, Δt) where T
    state = body.state

    state.xc = state.xsol[2]
    state.qc = state.qsol[2]
    state.vc = state.vsol[2]
    state.ωc = state.ωsol[2]

    state.xk[1] = state.xk[1] + state.vsol[2]*Δt
    state.qk[1] = state.qk[1] * ωbar(state.ωsol[2],Δt) * Δt / 2

    state.xsol[2] = state.xk[1]
    state.qsol[2] = state.qk[1]

    state.Fk[1] = szeros(T,3)
    state.τk[1] = szeros(T,3)
    return
end

@inline function setsolution!(body::Body)
    state = body.state
    state.xsol[2] = state.xk[1]
    state.qsol[2] = state.qk[1]
    state.vsol[1] = state.vc
    state.vsol[2] = state.vc
    state.ωsol[1] = state.ωc
    state.ωsol[2] = state.ωc
    return
end

@inline function settempvars!(body::Body{T}, x, v, F, q, ω, τ, d) where T
    state = body.state
    stateold = deepcopy(state)

    state.xc = x
    state.qc = q
    state.vc = v
    state.ωc = ω
    state.Fk[1] = F
    state.τk[1] = τ
    state.d = d

    return stateold
end

function ∂integration(q2::UnitQuaternion{T}, ω2::SVector{3,T}, Δt::T) where {T}
    Δ = Δt * SMatrix{3,3,T,9}(Diagonal(sones(T,3)))
    X = hcat(Δ, szeros(T,3,3))
    Q = hcat(szeros(T,4,3), Lmat(q2)*derivωbar(ω2, Δt)*Δt/2)
    return svcat(X, Q) # 7x6
end
