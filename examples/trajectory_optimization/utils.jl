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

function getMaxGradients(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}) where {T,Nn,Ne,Nb,Ni}
	Δt = mechanism.Δt
	nu = controldim(mechanism)
	attjac = false
	nic = attjac ? 12Nb : 13Nb
	neqcs = eqcdim(mechanism)
	datamat = full_data_matrix(mechanism, attjac = attjac)
	solmat = full_matrix(mechanism.system)

	∇data_z = - solmat \ datamat
	∇z_vϕ = ∇data_z[neqcs .+ (1:6Nb),1:nic]
	∇u_vϕ = ∇data_z[neqcs .+ (1:6Nb),nic .+ (1:nu)]
	∇z_z̄ = zeros(13Nb,13Nb)
	∇u_z̄ = zeros(13Nb,nu)
	for (i, body) in enumerate(mechanism.bodies)
		# Fill in gradients of v25, ϕ25
		∇z_z̄[13*(i-1) .+ [4:6; 11:13],:] += ∇z_vϕ[6*(i-1) .+ (1:6),:]
		∇u_z̄[13*(i-1) .+ [4:6; 11:13],:] += ∇u_vϕ[6*(i-1) .+ (1:6),:]

		# Fill in gradients of x3, q3
		q2 = body.state.q2[1]
		ϕ25 = body.state.ϕsol[2]
		∇z_z̄[13*(i-1) .+ (1:3),:] += ∂integrator∂v(Δt) * ∇z_vϕ[6*(i-1) .+ (1:3),:]
		∇z_z̄[13*(i-1) .+ (1:3),13*(i-1) .+ (1:13)] += ∂integrator∂x() * [I(3) zeros(3,10)]
		∇z_z̄[13*(i-1) .+ (7:10),:] += ∂integrator∂ϕ(q2, ϕ25, Δt) * ∇z_vϕ[6*(i-1) .+ (4:6),:]
		∇z_z̄[13*(i-1) .+ (7:10),13*(i-1) .+ (1:13)] += ∂integrator∂q(q2, ϕ25, Δt, attjac = false) * [zeros(4,6) I(4) zeros(4,3)]

		∇u_z̄[13*(i-1) .+ (1:3),:] += ∂integrator∂v(Δt) * ∇u_vϕ[6*(i-1) .+ (1:3),:]
		∇u_z̄[13*(i-1) .+ (7:10),:] += ∂integrator∂ϕ(q2, ϕ25, Δt) * ∇u_vϕ[6*(i-1) .+ (4:6),:]
	end
	return ∇z_z̄, ∇u_z̄
end


function getMaxGradients!(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, z::AbstractVector{T}, u::AbstractVector{T};
		opts=InteriorPointOptions()) where {T,Nn,Ne,Nb,Ni}
	step!(mechanism, z, u, opts=opts)
	∇z_z̄, ∇u_z̄ = getMaxGradients(mechanism)
	return ∇z_z̄, ∇u_z̄
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
				if typeof(joint) <: Rotational
					push!(c, minimalCoordinates(joint, qb)...) # Δq = q2 of body b
					push!(v, minimalVelocities(joint, qb, ϕb)...) # Δϕ = ϕ15 of bodyb in origin = world's coordinates
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
				nj = controldim(joint)
				push!(c, x[off .+ (1:nj)]...); off += nj
			end
			push!(v, x[off .+ (1:n)]...); off += n
			_setPosition!(mechanism, eqc, c)
			setVelocity!(mechanism, eqc, v) # assumes we provide v and ϕ in body1's coordinates i.e world coordinates
		end
	end
	z = getMaxState(mechanism)
	return z
end

function ∇min2max(mechanism::Mechanism, x)
	FiniteDiff.finite_difference_jacobian(x -> min2max(mechanism, x), x)
end

function ∇max2min(mechanism::Mechanism, z)
	FiniteDiff.finite_difference_jacobian(z -> max2min(mechanism, z), z)
end

function getMinGradients!(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, z::AbstractVector{T}, u::AbstractVector{T};
		opts=InteriorPointOptions()) where {T,Nn,Ne,Nb,Ni}

	step!(mechanism, z, u, opts=opts)
	z = getState(mechanism)
	z̄ = getNextState(mechanism)
	x = max2min(mechanism, z)
	x̄ = max2min(mechanism, z̄)
	∇z_z̄, ∇u_z̄ = getMaxGradients(mechanism)

	∇x = ∇min2max(mechanism, x)
	∇z̄ = ∇max2min(mechanism, z̄)
	∇x_x̄ = ∇z̄ * ∇z_z̄ * ∇x
	∇u_x̄ = ∇z̄ * ∇u_z̄
	return ∇x_x̄, ∇u_x̄
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

function getMinState(mechanism::Mechanism{T,Nn,Ne,Nb,Ni};
	pos_noise=nothing, vel_noise=nothing,
	pos_noise_range=[-Inf, Inf], vel_noise_range=[-3.9 / mechanism.Δt^2, 3.9 / mechanism.Δt^2]) where {T,Nn,Ne,Nb,Ni}
	# z = [[x2, v15, q2, ϕ15]body1  [x2, v15, q2, ϕ15]body2 ....]
	x = []
	# When we set the Δv and Δω in the mechanical graph, we need to start from the root and get down to the leaves.
	# Thus go through the eqconstraints in order, start from joint between robot and origin and go down the tree.
	for id in reverse(mechanism.system.dfs_list)
		(id > Ne) && continue # only treat eqconstraints
		eqc = mechanism.eqconstraints[id]
		c = zeros(T,0)
		v = zeros(T,0)
		for (i,joint) in enumerate(eqc.constraints)
			cbody = getbody(mechanism, eqc.childids[i])
			pbody = getbody(mechanism, eqc.parentid)
			# push!(c, minimalCoordinates(joint, pbody, cbody)...)
			# push!(v, minimalVelocities(joint, pbody, cbody)...)
			pos = minimalCoordinates(joint, pbody, cbody)
			vel = minimalVelocities(joint, pbody, cbody)
			if pos_noise != nothing
				pos += clamp.(length(pos) == 1 ? rand(pos_noise, length(pos))[1] : rand(pos_noise, length(pos)), pos_noise_range...)
			end
			if vel_noise != nothing
				vel += clamp.(length(vel) == 1 ? rand(vel_noise, length(vel))[1] : rand(vel_noise, length(vel)), vel_noise_range...)
			end
			push!(c, pos...)
			push!(v, vel...)
		end
		push!(x, [c; v]...)
	end
	x = [x...]
	return x
end

function velocity_index(mechanism::Mechanism{T,Nn,Ne}) where {T,Nn,Ne}
    ind = []
    off = 0
    for id in reverse(mechanism.system.dfs_list)
        (id > Ne) && continue # only treat eqconstraints
        eqc = mechanism.eqconstraints[id]
        nu = controldim(eqc)
        push!(ind, Vector(off + nu .+ (1:nu)))
        off += 2nu
    end
    return vcat(ind...)
end

function unpackMaxState(z::AbstractVector, i::Int)
	zi = z[(i-1)*13 .+ (1:13)]
	x2 = zi[1:3]
	v15 = zi[4:6]
	q2 = UnitQuaternion(zi[7:10]..., false)
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

function setSpringOffset!(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, x::AbstractVector) where {T,Nn,Ne,Nb,Ni}
	# When we set the Δv and Δω in the mechanical graph, we need to start from the root and get down to the leaves.
	# Thus go through the eqconstraints in order, start from joint between robot and origin and go down the tree.
	off = 0
	for id in reverse(mechanism.system.dfs_list)
		(id > Ne) && continue # only treat eqconstraints
		eqc = mechanism.eqconstraints[id]
		for (i,joint) in enumerate(eqc.constraints)
			cbody = getbody(mechanism, eqc.childids[i])
			pbody = getbody(mechanism, eqc.parentid)
			N̄ = 3 - length(joint)
			joint.spring_offset = x[off .+ (1:N̄)]
			off += 2N̄
		end
	end
	return nothing
end

function gravity_compensation(mechanism::Mechanism)
    # only works with revolute joints for now
    nu = controldim(mechanism)
    u = zeros(nu)
    off  = 0
    for eqc in mechanism.eqconstraints
        nu = controldim(eqc)
        if eqc.parentid != nothing
            body = getbody(mechanism, eqc.parentid)
            rot = eqc.constraints[2]
            A = Matrix(nullspacemat(rot))
            Fτ = springforce(mechanism, eqc, body)
            F = Fτ[1:3]
            τ = Fτ[4:6]
            u[off .+ (1:nu)] = -A * τ
        else
            @warn "need to treat the joint to origin"
        end
        off += nu
    end
    return u
end

function inverse_control(mechanism::Mechanism, x, x̄; ϵtol = 1e-5)
	nu = controldim(mechanism)
	u = zeros(nu)
	# starting point of the local search
	for k = 1:10
		err = inverse_control_error(mechanism, x, x̄, u, ϵtol = ϵtol)
		norm(err, Inf) < 1e-10 && continue
		∇ = FiniteDiff.finite_difference_jacobian(u -> inverse_control_error(mechanism, x, x̄, u, ϵtol = ϵtol), u)
		u -= ∇ \ err
	end
	return u
end

function inverse_control_error(mechanism, x, x̄, u; ϵtol = 1e-5)
	z = min2max(mechanism, x)
	z̄ = min2max(mechanism, x̄)
	setState!(mechanism, z)
	opts = InteriorPointOptions(rtol=ϵtol, btol=ϵtol, undercut=1.5)
	err = x̄ - max2min(mechanism, step!(mechanism, min2max(mechanism, x), u, opts=opts))
	return err
end
