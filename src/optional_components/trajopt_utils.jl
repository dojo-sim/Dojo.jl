function generate_storage(mechanism, x)
    steps = length(x)
    nbodies = length(mechanism.bodies)
    storage = Storage{Float64}(steps, nbodies)

    for t = 1:steps
        off = 0
        for (i, body) in enumerate(mechanism.bodies)
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
	# When we set the Δv and Δω in the mechanical graph, we need to start from the root and get down to the leaves.
	# Thus go through the eqconstraints in order, start from joint between robot and origin and go down the tree.
	for id in reverse(mechanism.system.dfs_list)
		(id > Ne) && continue # only treat eqconstraints
		eqc = mechanism.eqconstraints[id]
		c = zeros(Tz,0)
		v = zeros(Tz,0)
		for (i,joint) in enumerate(eqc.constraints)
			ichild = eqc.childids[i] - Ne
			xb, vb, qb, ϕb = unpackMaxState(z, ichild)
			if eqc.parentid != nothing
				iparent = eqc.parentid - Ne
				xa, va, qa, ϕa = unpackMaxState(z, iparent)
				if typeof(joint) <: Translational
					push!(c, minimalCoordinates(joint, xa, qa, xb, qb)...) # Δx in bodya's coordinates projected on jointAB's nullspace
					push!(v, minimalVelocities(joint, xa, qa, va, ϕa, xb, qb, vb, ϕb)...) # Δv in bodya's coordinates projected on jointAB's nullspace
				elseif typeof(joint) <: Rotational
					push!(c, minimalCoordinates(joint, qa, qb)...) # Δq in bodya's coordinates projected on jointAB's nullspace
					push!(v, minimalVelocities(joint, qa, ϕa, qb, ϕb)...) # Δϕ in bodya's coordinates projected on jointAB's nullspace
				end
			else
				# we need a special case: when the first link has free rotation wrt the origin Rotational0
				#TODO @warn "this is special cased for floating-base"
				if typeof(joint) <: Rotational0 #
					push!(c, minimalCoordinates(joint, qb)...) # Δq = q2 of body b
					# ϕw = minimalVelocities(joint, qb, ϕb) # Δϕ = ϕ15 of bodyb in world's coordinates
					# # @show ϕ
					# @show ϕw
					# push!(v, ϕw...) # Δϕ = ϕ15 of bodyb in origin = world's coordinates
					push!(v,  minimalVelocities(joint, qb, ϕb)...) # Δϕ = ϕ15 of bodyb in origin = world's coordinates
				elseif typeof(joint) <: Rotational
					push!(c, minimalCoordinates(joint, qb)...) # Δq = q2 of body b
					push!(v, minimalVelocities(joint, qb, ϕb)...) # Δϕ = ϕ15 of bodyb in bodyb's coordinates
				elseif typeof(joint) <: Translational
					push!(c, minimalCoordinates(joint, xb, qb)...) # Δx = x2 of bodyb in world's coordinates
					push!(v, minimalVelocities(joint, qb, vb, ϕb)...) # Δv = v15 of bodyb in world's coordinates
				end
			end
		end
		push!(x, [c; v]...)
	end
	x = [x...]
	return x
end

function min2max(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, x::AbstractVector{Tx}) where {T,Nn,Ne,Nb,Ni,Tx}
	# x = [[x2, q2, v15, ϕ15]body1  θ_body1-body2 ....]
	# z = [[x2, v15, q2, ϕ15]body1  [x2, v15, q2, ϕ15]body2 ....]
	# When we set the Δv and Δω in the mechanical graph, we need to start from the root and get down to the leaves.
	# Thus go through the eqconstraints in order, start from joint between robot and origin and go down the tree.
	off = 0
	for id in reverse(mechanism.system.dfs_list)
		(id > Ne) && continue # only treat eqconstraints
		eqc = mechanism.eqconstraints[id]
		n = controldim(eqc)
		if eqc.parentid != nothing
			c = x[off .+ (1:n)]; off += n
			v = x[off .+ (1:n)]; off += n#in body1
			_setPosition!(mechanism, eqc, c)
			setVelocity!(mechanism, eqc, v)#in body1
		else
			@assert length(Set(eqc.childids)) == 1 # only one body is linked to the origin
			c = zeros(Tx,0)
			v = zeros(Tx,0)
			# we need a special case: when the first link has free rotation wrt the origin
			q2 = one(UnitQuaternion) # only use for the special case
			for joint in eqc.constraints
				if typeof(joint) <: Rotational0
					nj = 4 # 4 instead of 3 since we use a quaternion for the first link in the body, when there is no rotational constraint
					q2 = UnitQuaternion(x[off .+ (1:nj)]..., false)
					push!(c, rotation_vector(q2)...); off += nj
				else
					nj = controldim(joint)
					push!(c, x[off .+ (1:nj)]...); off += nj
				end
			end

			push!(v, x[off .+ (1:n)]...); off += n
			_setPosition!(mechanism, eqc, c)
			setVelocity!(mechanism, eqc, v) # assumes we provide v and ϕ in body1's coordinates i.e world coordinates
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
	ϕ15 = zi[11:13]
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

function visualizeMaxCoord(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, z::AbstractVector, vis::Visualizer) where {T,Nn,Ne,Nb,Ni}
	storage = Storage(1,Nb)
	for t = 1:1
		for b = 1:Nb
			x2, v15, q2, ϕ15 = unpackMaxState(z, b)
			storage.x[b][t] = x2
			storage.v[b][t] = v15
			storage.q[b][t] = q2
			storage.ω[b][t] = ϕ15
		end
	end
	visualize(mechanism, storage, vis = vis)
end
