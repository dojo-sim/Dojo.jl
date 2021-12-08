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

function initialize!(state::State, order)
    state.order = order

    state.x2 = [state.x2[1] for i = 1:order]
    state.q2 = [state.q2[1] for i = 1:order]
    state.F2 = [state.F2[1] for i = 1:order]
    state.τ2 = [state.τ2[1] for i = 1:order]

    return
end

function setState!(mechanism::Mechanism, z::AbstractVector)
    Δt = mechanism.Δt
    off = 0
    for body in mechanism.bodies
        x2, v15, q2, ϕ15 = unpackdata(z[off+1:end]); off += 13
        q2 = UnitQuaternion(q2..., false)
        body.state.v15 = v15
        body.state.ϕ15 = ϕ15
        body.state.x2[1] = x2
        body.state.q2[1] = q2
		discretizestate!(mechanism) #set x1, q1 and zeroes out F2 τ2
    end
	foreach(setsolution!, mechanism.bodies) # warm-start solver
end

function setControl!(mechanism::Mechanism{T}, u::AbstractVector) where {T}
	eqcs = mechanism.eqconstraints
	# set the controls in the equality constraints
	off = 0
	for eqc in eqcs
		nu = controldim(eqc)
		setForce!(mechanism, eqc, SVector{nu,T}(u[off .+ (1:nu)]))
		off += nu
	end
	# apply the controls to each body's state
	foreach(applyFτ!, eqcs, mechanism)
end

function getState(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}) where {T,Nn,Ne,Nb,Ni}
	z = zeros(T,13Nb)
	for (i, body) in enumerate(mechanism.bodies)
		v15 = body.state.v15
		ϕ15 = body.state.ϕ15
		x2 = body.state.x2[1]
		q2 = body.state.q2[1]
		z[13*(i-1) .+ (1:13)] = [x2; v15; vector(q2); ϕ15]
	end
	return z
end

function getNextState(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}) where {T,Nn,Ne,Nb,Ni}
	Δt = mechanism.Δt
	z̄ = zeros(T,13Nb)
	for (i, body) in enumerate(mechanism.bodies)
		v25 = body.state.vsol[2]
		ϕ25 = body.state.ϕsol[2]
		x3 = getx3(body.state, Δt)
		q3 = getq3(body.state, Δt)
		z̄[13*(i-1) .+ (1:13)] = [x3; v25; vector(q3); ϕ25]
	end
	return z̄
end