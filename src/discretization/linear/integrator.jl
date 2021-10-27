## Uses linear interpolation between the knot points. For a = b = 1/2 the integration follows the trapezoid rule.
## This results in the Symplectic Euler method or Stoermer-Verlet, depending on how exactly one sets it up.
# u(t) = a t + b = (xdk+1 - xdk)/Δt t + xdk
# du(t) = a = (xdk+1 - xdk)/Δt
# ωc(t) = 2 V L(qc)ᵀ du(t)
# L(xck,vck) -> Δt (a Ld(u(0),du(0)) + b Ld(u(h),u(h)), where a + b = 1
# L(qck,ωck) -> Δt (a Ld(u(0),2 V qdk† (qdk+1-qdk)/Δt) + b Ld(qdk+1,2 V qdk† (qdk+1-qdk)/Δt), where a + b = 1
# ωckw = sqrt((2/Δt)^2 - ωckᵀωck) - 2/Δt
# Fckᵀxck -> a Fdkᵀxdk + b Fdk+1ᵀxdk, where a + b = 1

METHODORDER = 1 # This refers to the interpolating spline
getGlobalOrder() = (global METHODORDER; return METHODORDER)

# Convenience functions
@inline getx3(state::State, Δt) = state.xk[1] + state.vsol[2]*Δt
@inline getx25(state::State, Δt) = state.xk[1] + state.vsol[2] * Δt/2
@inline getq3(state::State, Δt) = state.qk[1] * ωbar(state.ωsol[2],Δt) * Δt / 2
@inline getq25(state::State, Δt) = state.qk[1] * ωbar(state.ωsol[2]/2,Δt) * Δt / 2

@inline posargsc(state::State) = (state.xc, state.qc)
@inline fullargsc(state::State) = (state.xc, state.vc, state.qc, state.ωc)
@inline posargsk(state::State; k=1) = (state.xk[k], state.qk[k])
@inline posargssol(state::State) = (state.xsol[2], state.qsol[2])
@inline fullargssol(state::State) = (state.xsol[2], state.vsol[2], state.qsol[2], state.ωsol[2])
@inline posargsnext(state::State, Δt) = (getx3(state, Δt), getq3(state, Δt))
@inline posargshalf(state::State, Δt) = (getx25(state, Δt), getq3(state, Δt))



# q0 = UnitQuaternion(rand(4)...)
# ω0 = [20., 0, 0]
# Δt = 0.0999
# q1 = q0 * ωbar(ω0, Δt) * Δt / 2
# qh = q0 * ωbar(ω0/2, Δt) * Δt / 2
# qh1 = q0 * ωbar(ω0, Δt/2) * Δt / 4

# qh_ = UnitQuaternion(midpoint([q0.w, q0.x, q0.y, q0.z], [q1.w, q1.x, q1.y, q1.z])..., false)
# norm([qh.w, qh.x, qh.y, qh.z])
# norm([qh.w, qh.x, qh.y, qh.z] - [qh_.w, qh_.x, qh_.y, qh_.z])
# norm([qh1.w, qh1.x, qh1.y, qh1.z] - [qh_.w, qh_.x, qh_.y, qh_.z])
# δq = (qh_ * qh')
# [δq.w, δq.x, δq.y, δq.z]

# 2π / 0.23

# Lmat(sδq)
# Square root of quaternion: https://www.johndcook.com/blog/2021/01/06/quaternion-square-roots/
function sqrt_quat(q; ϵ=1e-16)
	r = norm(q)
	# Angle
	theta = acos(q[1]/r)
	# Axis
	u = q[2:4]
	u ./= norm(u) + ϵ
	# Half axis-angle rotation
	x = [cos(theta/2);  sin(theta/2)*u]
	# Sqrt on the norm of the quaternion
	x .*= r^0.5 # useless for unit quaternion
	return x
end

function midpoint(q0, q1)
	# Small delta rotation
	qψ = R_multiply(q0)' * q1
	# Divide delta in two
	sqψ = sqrt_quat(qψ)
	# Apply small rotation
	qmid = L_multiply(sqψ) * q0
	return qmid
end

function L_multiply(q)
	s = q[1]
	v = q[2:4]

	SMatrix{4,4}([s -transpose(v);
	              v s * I + skew(v)])
end

function R_multiply(q)
	s = q[1]
	v = q[2:4]

	SMatrix{4,4}([s -transpose(v);
	              v s * I - skew(v)])
end













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
