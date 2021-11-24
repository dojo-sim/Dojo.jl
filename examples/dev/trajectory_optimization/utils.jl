function step!(mech::Mechanism, x, u;
    btol=1.0e-4, undercut=Inf,
    control_inputs = (a, b, c) -> nothing)

    # set data
    data = [x; u]

    off = 0

    for body in mech.bodies
        x2, v15, q2, ω15 = unpackdata(data[off+1:end]); off += 13
        body.state.x1 = x2 - v15 * mech.Δt
        body.state.v15 = v15
        body.state.q1 = UnitQuaternion(q2...) * ωbar(-ω15, mech.Δt) * mech.Δt / 2.0
        body.state.ϕ15 = ω15
    end

    # controller
    controller!(mech, k) = control_inputs(mech, k, u)

    # simulate
    storage = simulate!(mech, mech.Δt,
        controller!, record=true, verbose=false, solver=:mehrotra!, btol=btol, undercut=undercut)

    # next state
    x_next = []

    for body in mech.bodies
        x3 = body.state.x2[1]
        v25 = body.state.vsol[2]
        _q3 = body.state.q2[1]
        q3 = [_q3.w; _q3.x; _q3.y; _q3.z]
        ω25 = body.state.ϕsol[2]
        push!(x_next, [x3; v25; q3; ω25]...)
    end

    return x_next
end

function step_grad_x!(mech::Mechanism, x, u;
    btol=1.0e-3, undercut=1.5,
    control_inputs = (a, b, c) -> nothing)
    m = length(u)
    ∂step∂x = fdjac(w -> step!(mech, w[1:(end-m)], w[(end-m+1):end],
        btol=btol, undercut=undercut, control_inputs=control_inputs),
        [x; u])[:, 1:(end-m)]

    return ∂step∂x
end

function step_grad_u!(mech::Mechanism, x, u;
    btol=1.0e-3, undercut=1.5,
    control_inputs = (a, b, c) -> nothing)
    m = length(u)
    ∂step∂u = fdjac(w -> step!(mech, w[1:(end-m)], w[(end-m+1):end],
        btol=btol, undercut=undercut, control_inputs=control_inputs),
        [x; u])[:, (end-m+1):end]

    return ∂step∂u
end

function generate_storage(mech, x)
    steps = length(x)
    nbodies = length(mech.bodies)
    storage = Storage{Float64}(steps, nbodies)

    for t = 1:steps
        off = 0
        for (i, body) in enumerate(mech.bodies)
            storage.x[i][t] = x[t][off .+ (1:3)]
            storage.v[i][t] = x[t][off .+ (4:6)]
            storage.q[i][t] = UnitQuaternion(x[t][off .+ (7:10)]...)
            storage.ω[i][t] = x[t][off .+ (11:13)]
            off += 13
        end
    end

    return storage
end





function fast_grad_x!(mech::Mechanism{T,Nn,Ne,Nb}, x, u; btol=1.0e-3, undercut=1.5,
        control_inputs = (a, b, c) -> nothing) where {T,Nn,Ne,Nb}

    step!(mech, x, u; btol = btol, undercut = undercut, control_inputs = control_inputs)

    ∂step∂x = full_matrix(mech.system) * full_data_matrix(mech)[:,1:12Nb]
    return ∂step∂x
end


function fast_grad_u!(mech::Mechanism{T,Nn,Ne,Nb}, x, u; btol=1.0e-3, undercut=1.5,
        control_inputs = (a, b, c) -> nothing) where {T,Nn,Ne,Nb}

    step!(mech, x, u; btol = btol, undercut = undercut, control_inputs = control_inputs)

    ∂step∂u = full_matrix(mech.system) * full_data_matrix(mech)[:,1:12Nb]
    return ∂step∂u
end


#
# using BenchmarkTools
# # @benchmark
# full_data_matrix(mech)
# @benchmark step_grad_x!(mech, z[1], u_control)
# @benchmark fast_grad_x!(mech, z[1], u_control)






















function simon_step!(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, x::AbstractVector{T}, u::AbstractVector{T};
		ϵ::T = 1e-6, newtonIter::Int = 100, lineIter::Int = 10, verbose::Bool = true,
		btol::T = ϵ, undercut::T = Inf) where {T,Nn,Ne,Nb,Ni}

	# set the initial conditions: x1, v15, x2...
		# set x2, v15, q2, ϕ15
		# set x1, q1
		# set F2 and τ2 to zero
		# warm-start the solver
	setState!(mechanism, x)

	# set the controls in the equality constraints
		# apply the controls to each body's state
	setControl!(mechanism, u)

	# solve the 1 step simulation problem
	mehrotra!(mechanism, ϵ = ϵ, newtonIter = newtonIter, lineIter = lineIter, verbose = verbose,
		opts=InteriorPointOptions(rtol=ϵ, max_iter=newtonIter, btol=btol, undercut=undercut, verbose=verbose))

	# extract the next state
	x̄ = getNextState(mechanism)
    return x̄
end

function setState!(mechanism::Mechanism, x::AbstractVector)
    Δt = mechanism.Δt
    off = 0
    for body in mechanism.bodies
        x2, v15, q2, ϕ15 = unpackdata(x[off+1:end]); off += 13
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

function getNextState(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}) where {T,Nn,Ne,Nb,Ni}
	Δt = mechanism.Δt
	x̄ = zeros(T,13Nb)
	for (i, body) in enumerate(mechanism.bodies)
		v25 = body.state.vsol[2]
		ϕ25 = body.state.ϕsol[2]
		x3 = getx3(body.state, Δt)
		q3 = getq3(body.state, Δt)
		x̄[13*(i-1) .+ (1:13)] = [x3; v25; vector(q3); ϕ25]
	end
	return x̄
end

function getGradients(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}) where {T,Nn,Ne,Nb,Ni}
	Δt = mechanism.Δt
	nu = controldim(mechanism)
	attjac = false
	nic = attjac ? 12Nb : 13Nb
	neqcs = eqcdim(mechanism)
	datamat = full_data_matrix(mechanism, attjac = attjac)
	solmat = full_matrix(mechanism.system)

	∇data_z = - solmat \ datamat
	∇x_vϕ = ∇data_z[neqcs .+ (1:6Nb),1:nic]
	∇u_vϕ = ∇data_z[neqcs .+ (1:6Nb),nic .+ (1:nu)]
	∇x_x̄ = zeros(13Nb,13Nb)
	∇u_x̄ = zeros(13Nb,nu)
	for (i, body) in enumerate(mechanism.bodies)
		# Fill in gradients of v25, ϕ25
		∇x_x̄[13*(i-1) .+ [4:6; 11:13],:] += ∇x_vϕ[6*(i-1) .+ (1:6),:]
		∇u_x̄[13*(i-1) .+ [4:6; 11:13],:] += ∇u_vϕ[6*(i-1) .+ (1:6),:]

		# Fill in gradients of x3, q3
		q2 = body.state.q2[1]
		ϕ25 = body.state.ϕsol[2]
		∇x_x̄[13*(i-1) .+ (1:3),:] += ∂integrator∂v(Δt) * ∇x_vϕ[6*(i-1) .+ (1:3),:]
		∇x_x̄[13*(i-1) .+ (1:3),13*(i-1) .+ (1:13)] += ∂integrator∂x() * [I(3) zeros(3,10)]
		∇x_x̄[13*(i-1) .+ (7:10),:] += ∂integrator∂ϕ(q2, ϕ25, Δt) * ∇x_vϕ[6*(i-1) .+ (4:6),:]
		∇x_x̄[13*(i-1) .+ (7:10),13*(i-1) .+ (1:13)] += ∂integrator∂q(q2, ϕ25, Δt, attjac = false) * [zeros(4,6) I(4) zeros(4,3)]

		∇u_x̄[13*(i-1) .+ (1:3),:] += ∂integrator∂v(Δt) * ∇u_vϕ[6*(i-1) .+ (1:3),:]
		∇u_x̄[13*(i-1) .+ (7:10),:] += ∂integrator∂ϕ(q2, ϕ25, Δt) * ∇u_vϕ[6*(i-1) .+ (4:6),:]
	end
	return ∇x_x̄, ∇u_x̄
end
q20 = UnitQuaternion(rand(4)...)
ϕ250 = srand(3)
∂integrator∂q(q20, ϕ250, Δt, attjac = false)
[zeros(4,6) I(4) zeros(4,3)]




function getGradients!(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, x::AbstractVector{T}, u::AbstractVector{T};
		ϵ::T = 1e-6, newtonIter::Int = 100, lineIter::Int = 10, verbose::Bool = true,
		btol::T = ϵ, undercut::T = Inf) where {T,Nn,Ne,Nb,Ni}
	simon_step!(mechanism, x, u, ϵ = ϵ, newtonIter = newtonIter, lineIter = lineIter,
		verbose = verbose, btol = btol, undercut = undercut)
	∇x_x̄, ∇u_x̄ = getGradients(mechanism)
	return ∇x_x̄, ∇u_x̄
end
#
# mech = gethopper(Δt = Δt, g = gravity)
# x0 = hopper_initial_state()
# u0 = [0.0; 0.0; mech.g * mech.Δt]
# u0 = Δt*0.1*Vector(1:1.0:7)
# simon_step!(mech, x0, u0, ϵ = 1e-6, btol = 1e-6, undercut = Inf, verbose = true)
# @benchmark gradx0, gradu0 = getGradients(mech)
# gradx0, gradu0 = getGradients!(mech, x0, u0, ϵ = 1e-6, btol = 1e-6, undercut = Inf, verbose = true)
# using BenchmarkTools
# @benchmark gradx1 = FiniteDiff.finite_difference_jacobian(
# 	x -> simon_step!(mech, x, u0, ϵ = 1e-6, btol = 1e-6, undercut = Inf, verbose = false),
# 	x0)
# gradu1 = FiniteDiff.finite_difference_jacobian(
# 	u -> simon_step!(mech, x0, u, ϵ = 1e-6, btol = 1e-6, undercut = Inf, verbose = false),
# 	u0)
# norm(gradx0 - gradx1, Inf) / norm(gradx1, Inf)
# norm(gradu0 - gradu1, Inf) / norm(gradu1, Inf)
#
# plot(Gray.(abs.(gradx0)))
# plot(Gray.(abs.(gradx1)))
# plot(Gray.(abs.(gradu0)))
# plot(Gray.(abs.(gradu1)))
#
#
#
#
# mech = getmechanism(:pendulum, Δt = Δt, g = gravity)
# x0 = [0;0;-0.5; zeros(3); vector(one(UnitQuaternion)); zeros(3)]
# u0 = [0.01]
# simon_step!(mech, x0, u0, ϵ = 1e-6, btol = 1e-6, undercut = Inf, verbose = true)
# gradx0, gradu0 = getGradients(mech)
# gradx0, gradu0 = getGradients!(mech, x0, u0, ϵ = 1e-6, btol = 1e-6, undercut = Inf, verbose = true)
#
# gradx1 = FiniteDiff.finite_difference_jacobian(
# 	x -> simon_step!(mech, x, u0, ϵ = 1e-6, btol = 1e-6, undercut = Inf, verbose = true),
# 	x0)
# gradu1 = FiniteDiff.finite_difference_jacobian(
# 	u -> simon_step!(mech, x0, u, ϵ = 1e-6, btol = 1e-6, undercut = Inf, verbose = true),
# 	u0)
#
# plot(Gray.(abs.(gradx0)))
# plot(Gray.(abs.(gradx1)))
# plot(Gray.(abs.(gradu0)))
# plot(Gray.(abs.(gradu1)))
#
# norm(gradx0 - gradx1) / norm(gradx1)
# norm(gradu0 - gradu1) / norm(gradu1)
#
# round.(gradx1, digits=3)
#
# round.(gradx0, digits=3)
# norm((gradx0 - gradx1)[1:3,:])
# norm((gradx0 - gradx1)[4:6,:])
# norm((gradx0 - gradx1)[7:10,:])
# norm((gradx0 - gradx1)[11:13,:])
#
# norm((gradx0 - gradx1)[1:3,1:3])
# norm((gradx0 - gradx1)[1:3,4:13])
#
# norm((gradx0 - gradx1)[7:10,1:6])
# norm((gradx0 - gradx1)[7:10,7:10])
# norm((gradx0 - gradx1)[7:10,11:13])
#
# gradx0[1:3,1:3]
# gradx1[1:3,1:3]
