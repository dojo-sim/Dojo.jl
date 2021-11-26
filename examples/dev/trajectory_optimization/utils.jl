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

function getGradients!(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, x::AbstractVector{T}, u::AbstractVector{T};
		ϵ::T = 1e-6, newtonIter::Int = 100, lineIter::Int = 10, verbose::Bool = true,
		btol::T = ϵ, undercut::T = Inf) where {T,Nn,Ne,Nb,Ni}
	simon_step!(mechanism, x, u, ϵ = ϵ, newtonIter = newtonIter, lineIter = lineIter,
		verbose = verbose, btol = btol, undercut = undercut)
	∇x_x̄, ∇u_x̄ = getGradients(mechanism)
	return ∇x_x̄, ∇u_x̄
end

function max2min(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, z::AbstractVector{Tz}) where {T,Nn,Ne,Nb,Ni,Tz}
	# z = [[x2, v15, q2, ϕ15]body1  [x2, v15, q2, ϕ15]body2 ....]
	x = []
	for eqc in mechanism.eqconstraints
		c = zeros(Tz,0)
		v = zeros(Tz,0)
		for (i,joint) in enumerate(eqc.constraints)
			ichild = eqc.childids[i] - Ne
			xb, vb, qb, ϕb = unpackMaxState(z, ichild)
			if eqc.parentid != nothing
				iparent = eqc.parentid - Ne
				xa, va, qa, ϕa = getMaxState(z, iparent)
				if typeof(joint) <: Translational
					push!(c, minimalCoordinates(joint, xa, qa, xb, qb)...)
					push!(v, minimalVelocities(joint, va, vb)...)
				elseif typeof(joint) <: Rotational
					push!(c, minimalCoordinates(joint, qa, qb)...)
					push!(v, minimalVelocities(joint, qa, ϕa, qb, ϕb)...)
				end
			else
				if typeof(joint) <: Translational
					push!(c, minimalCoordinates(joint, xb, qb)...)
					push!(v, minimalVelocities(joint, vb)...)
				elseif typeof(joint) <: Rotational
					push!(c, minimalCoordinates(joint, qb)...)
					push!(v, minimalVelocities(joint, qb, ϕb)...)
				end
			end
		end
		# @show c
		# @show v
		push!(x, [c; v]...)
	end
	x = [x...]
	return x
end

function min2max(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, x::AbstractVector{Tx}) where {T,Nn,Ne,Nb,Ni,Tx}
	# z = [[x2, v15, q2, ϕ15]body1  [x2, v15, q2, ϕ15]body2 ....]
	off = 0
	for eqc in mechanism.eqconstraints
		n = controldim(eqc)
		if eqc.parentid != nothing
			c = x[off .+ (1:n)]; off += n
			v = x[off .+ (1:n)]; off += n
			_setPosition!(mechanism, eqc, c)
			setVelocity!(mechanism, eqc, v)
		else
			@assert length(Set(eqc.childids)) == 1 # only one body is linked to the origin
			for ichild in Set(eqc.childids)
				@warn "this is special cased for floating-base"
				cv = x[off .+ (1:13)]; off += 13
				x2 = cv[1:3]
				q2 = UnitQuaternion(cv[4:7]...)
				c = [x2; rotation_vector(q2)]
				v = cv[8:13]
				_setPosition!(mechanism, eqc, c)
				setVelocity!(mechanism, eqc, v)
			end
		end
	end
	z = getMaxState(mechanism)
	return z
end

function getMaxState(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}) where {T,Nn,Ne,Nb,Ni}
	z = zeros(T, 13Nb)
	for (i, body) in enumerate(mechanism.bodies)
		x2 = body.state.x2[1]
		v15 = body.state.v15
		q2 = body.state.q2[1]
		ϕ15 = body.state.ϕ15
		setMaxState!(z, x2, v15, q2, ϕ15, i)
	end
	return z
end

function unpackMaxState(z::AbstractVector, i::Int)
	zi = z[(i-1)*13 .+ (1:13)]
	x2 = zi[1:3]
	v15 = zi[4:6]
	q2 = UnitQuaternion(zi[7:10]...)
	ϕ15 = zi[4:6]
	return x2, v15, q2, ϕ15
end

function setMaxState!(z::AbstractVector, x2::AbstractVector, v15::AbstractVector,
		q2::UnitQuaternion, ϕ15::AbstractVector, i::Int)
	z[(i-1)*13 .+ (1:3)] = x2
	z[(i-1)*13 .+ (4:6)] = v15
	z[(i-1)*13 .+ (7:10)] = vector(q2)
	z[(i-1)*13 .+ (11:13)] = ϕ15
	return nothing
end

function hopper_initial_state()
    # initial state
    x2b1 = [0.0; 0.0; 0.5]
    v15b1 = [0.0; 0.0; 0.0]
    q2b1 = [1.0; 0.0; 0.0; 0.0]
    ϕ15b1 = [0.0; 0.0; 0.0]
    z1b1 = [x2b1; v15b1; q2b1; ϕ15b1]

    x2b2 = [0.0; 0.0; 0.0]
    v15b2 = [0.0; 0.0; 0.0]
    q2b2 = [1.0; 0.0; 0.0; 0.0]
    ϕ15b2 = [0.0; 0.0; 0.0]
    z1b2 = [x2b2; v15b2; q2b2; ϕ15b2]
    z1 = [z1b1; z1b2]
end

# # Open visualizer
# vis = Visualizer()
# open(vis)

Δt0 = 0.01
g0 = -0.81
mech = getmechanism(:hopper, Δt = Δt0, g = g0, spring = 1.0, damper = 10.0)
initialize!(mech, :hopper)
function controller!(mechanism, k,)
    for (i,eqc) in enumerate(collect(mechanism.eqconstraints))
        nu = controldim(eqc)
		@show nu
        setForce!(mechanism, eqc, SVector{nu}(Δt0 * (rand(nu) .- 0.5) ))
    end
    return
end
@elapsed storage = simulate!(mech, 4.00, controller!, record = true, solver = :mehrotra!, verbose = false)
visualize(mech, storage, vis = vis)

Nb = length(storage.x)
Nt = length(storage.x[1])
err = []
for t = 1:Nt
	z0 = zeros(13Nb)
	for i = 1:Nb
		x2 = storage.x[i][t]
		v15 = storage.v[i][t]
		q2 = storage.q[i][t]
		ϕ15 = storage.ω[i][t]
		z0[13*(i-1) .+ (1:13)] = [x2; v15; vector(q2); ϕ15]
	end
	x0 = max2min(mech, z0)
	z1 = min2max(mech, x0)
	push!(err, norm(z0[1:3] - z1[1:3], Inf))
	push!(err, norm(z0[13 .+ (1:3)] - z1[13 .+ (1:3)], Inf))

	push!(err, norm(z0[4:6] - z1[4:6], Inf))
	# push!(err, norm(z0[13 .+ (4:6)] - z1[13 .+ (4:6)], Inf))

	qa0 = UnitQuaternion(z0[7:10]...)
	qa1 = UnitQuaternion(z1[7:10]...)
	δqa = qa0 * inv(qa1)
	qb0 = UnitQuaternion(z0[13 .+ (7:10)]...)
	qb1 = UnitQuaternion(z1[13 .+ (7:10)]...)
	δqb = qa0 * inv(qb1)
	push!(err, norm([δqa.x, δqa.y, δqa.z], Inf))
	push!(err, norm([δqb.x, δqb.y, δqb.z], Inf))

	# push!(err, norm(z0[11:13] - z1[11:13], Inf))
	# push!(err, norm(z0[13 .+ (11:13)] - z1[13 .+ (11:13)], Inf))
end
norm(err, Inf)
plot(err)

z0 = hopper_initial_state()
z0[1:3] = [0.1, 0.0, 0.5]
z0[14:16] = [0.1, 0.0, 0.5]
x0 = max2min(mech, z0)
z1 = min2max(mech, x0)
norm(z1 - z0, Inf)

@show z0[1:6]
@show z0[7:13]
@show z0[14:19]
@show z0[20:26]

@test norm(x0[1:3] - z0[1:3], Inf) < 1e-10
@test norm(x0[4:7] - z0[7:10], Inf) < 1e-10
@test norm(x0[8:10] - z0[4:6], Inf) < 1e-10
@test norm(x0[11:13] - z0[11:13], Inf) < 1e-10
@show x0[14:15]

@show z1[1:6]
@show z1[7:13]
@show z1[14:19]
@show z1[20:26]

@test norm((z0 - z1)[1:6], Inf) < 1e-10
@test norm((z0 - z1)[7:13], Inf) < 1e-10
@test norm((z0 - z1)[14:19], Inf) < 1e-10
@test norm((z0 - z1)[20:26], Inf) < 1e-10



jac = ForwardDiff.jacobian(z -> max2min(mech, z), z0)
plot(Gray.(abs.(jac)))
norm(jac * pinv(jac) - I(18), Inf)

mech
eqc1 = collect(mech.eqconstraints)[1]
xθ = rand(6)
_setPosition!(mech, eqc1, xθ)

eqc1.childids
a = 10
a = 10
a = 10
a = 10
a = 10
a = 10

# ARS demo on hopper (3h)
# fix min2max maxtomin (1h)
# test hopper halfcheetah with minimal coordinates (1h)


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
