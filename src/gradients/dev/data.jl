################################################################################
# Dimension
################################################################################
# Mechanism
data_dim(mechanism::Mechanism; attjac::Bool=true) =
	sum(Vector{Int64}(data_dim.(mechanism.eqconstraints))) +
    sum(Vector{Int64}(data_dim.(mechanism.bodies, attjac=attjac))) +
	sum(Vector{Int64}(data_dim.(mechanism.ineqconstraints)))
# Eqconstraints
data_dim(eqc::EqualityConstraint) = 2 + sum(data_dim.(eqc.constraints)) # [utra, urot, spring, damper, tra_spring_offset, rot_spring_offset]
data_dim(joint::Rotational{T,Nλ,Nb,N,Nb½,N̄λ}) where {T,Nλ,Nb,N,Nb½,N̄λ} = 2N̄λ # [u, spring, damper, spring_offset]
data_dim(joint::Translational{T,Nλ,Nb,N,Nb½,N̄λ}) where {T,Nλ,Nb,N,Nb½,N̄λ} = 2N̄λ # [u, spring, damper, spring_offset]
# Body
data_dim(body::Body; attjac::Bool=true) = attjac ? 19 : 20 # 1+6+6+6 or 1+6+6+7 [m,flat(J),x2,v15,q2,ϕ15] with attjac
# Ineqconstraints
data_dim(ineqc::InequalityConstraint) = sum(data_dim.(ineqc.constraints))
data_dim(bound::ContactBound) = 7 # [cf, p, offset]
data_dim(bound::LinearContactBound) = 7 # [cf, p, offset]
data_dim(bound::ImpactBound) = 6 # [p, offset]


################################################################################
# Attitude Jacobian
################################################################################
# Mechanism
function data_attitude_jacobian(mechanism::Mechanism)
	attjacs = [data_attitude_jacobian.(mechanism.eqconstraints);
		data_attitude_jacobian.(mechanism.bodies);
		data_attitude_jacobian.(mechanism.ineqconstraints)]
	attjac = cat(attjacs..., dims=(1,2))
	return attjac
end
# Eqconstraints
function data_attitude_jacobian(eqc::EqualityConstraint)
	return I(data_dim(eqc))
end
# Body
function data_attitude_jacobian(body::Body)
	# [m,flat(J),x1,q1,x2,q2]
	x2, q2 = posargs2(body.state)
	attjac = cat(I(1+6+6), G(vector(q2)), I(3), dims=(1,2))
	return attjac
end
# Ineqconstraints
function data_attitude_jacobian(ineqc::InequalityConstraint)
	return I(data_dim(ineqc))
end


################################################################################
# Get Data
################################################################################
# Mechanism
get_data(mechanism::Mechanism) = vcat([get_data.(mechanism.eqconstraints);
	get_data.(mechanism.bodies); get_data.(mechanism.ineqconstraints)]...)
# Eqconstraints
function get_data(eqc::EqualityConstraint)
	joints = eqc.constraints
	u = vcat(nullspacemat.(joints) .* getfield.(joints, :Fτ)...)
	spring = joints[1].spring # assumes we have the same spring and dampers for translational and rotational joint.
	damper = joints[1].damper # assumes we have the same spring and dampers for translational and rotational joint.
	spring_offset = vcat(getfield.(joints, :spring_offset)...)
	return [u; spring; damper; spring_offset]
end
# Body
function get_data(body::Body)
	m = body.m
	j = flatten_inertia(body.J)
	x2, q2 = posargs2(body.state)
	v15 = body.state.v15
	ϕ15 = body.state.ϕ15
	return [m; j; x2; v15; vector(q2); ϕ15]
end
# Ineqconstraints
get_data(bound::ContactBound) = [bound.cf; bound.offset; bound.p]
get_data(bound::LinearContactBound) = [bound.cf; bound.offset; bound.p]
get_data(bound::ImpactBound) = [bound.offset; bound.p]
get_data(ineqc::InequalityConstraint) = vcat(get_data.(ineqc.constraints)...)


################################################################################
# Set Data
################################################################################
# Mechanism
function set_data!(mechanism::Mechanism, data::AbstractVector)
	# It's important to treat bodies before eqcs
	# set_data!(body) will erase state.F2[1] and state.τ2[1]
	# set_data!(eqc) using applyFτ!, will write in state.F2[1] and state.τ2[1]
	c = 0
	for eqc in mechanism.eqconstraints
		Nd = data_dim(eqc)
		set_data!(mechanism, eqc, data[c .+ (1:Nd)]); c += Nd
	end
	for body in mechanism.bodies
		Nd = data_dim(body, attjac=false)
		set_data!(body, data[c .+ (1:Nd)], mechanism.Δt); c += Nd
	end
	for ineqc in mechanism.ineqconstraints
		Nd = data_dim(ineqc)
		set_data!(ineqc, data[c .+ (1:Nd)]); c += Nd
	end
	foreach(applyFτ!, mechanism.eqconstraints, mechanism, false)
	return nothing
end
 # Eqconstraints
function set_data!(mechanism::Mechanism, eqc::EqualityConstraint, data::AbstractVector)
	nu = controldim(eqc)
	u = data[SUnitRange(1,nu)]
	spring = data[nu+1]
	damper = data[nu+2]
	spring_offset = data[nu+2 .+ (1:nu)]

	setForce!(eqc, u)
	c = 0
	for joint in eqc.constraints
		nu = controldim(joint)
		joint.spring = spring
		joint.damper = damper
		joint.spring_offset = spring_offset[SUnitRange(c+1,c+nu)]; c += nu
	end
	# applyFτ!(eqc, mechanism, false)
	return nothing
end
# Body
function set_data!(body::Body, data::AbstractVector, Δt)
	# [m,flat(J),x2,v15,q2,ϕ15]
	m = data[1]
	J = lift_inertia(data[SUnitRange(2,7)])
	x2 = data[SUnitRange(8,10)]
	v15 = data[SUnitRange(11,13)]
	q2 = UnitQuaternion(data[14:17]..., false)
	ϕ15 = data[SUnitRange(18,20)]
	x1 = getx3(x2, -v15, Δt)
	q1 = getq3(q2, -ϕ15, Δt)

	body.m = m
	body.J = J
	body.state.x1 = x1
	body.state.v15 = v15
	body.state.q1 = q1
	body.state.ϕ15 = ϕ15
	body.state.x2[1] = x2
	body.state.q2[1] = q2
	body.state.F2[1] = SVector{3}(0,0,0.)
	body.state.τ2[1] = SVector{3}(0,0,0.)
	return nothing
end
# Ineqconstraints
function set_data!(bound::ContactBound, data::AbstractVector)
	bound.cf = data[1]
    bound.offset = data[SVector{3,Int}(2:4)]
    bound.p = data[SVector{3,Int}(5:7)]
    return nothing
end
function set_data!(bound::LinearContactBound, data::AbstractVector)
	bound.cf = data[1]
    bound.offset = data[SVector{3,Int}(2:4)]
    bound.p = data[SVector{3,Int}(5:7)]
    return nothing
end
function set_data!(bound::ImpactBound, data::AbstractVector)
    bound.offset = data[SVector{3,Int}(1:3)]
    bound.p = data[SVector{3,Int}(4:6)]
    return nothing
end
function set_data!(ineqc::InequalityConstraint, data::AbstractVector)
    c = 0
	for bound in ineqc.constraints
		N = data_dim(bound)
        set_data!(bound, data[c .+ (1:N)]); c += N
    end
    return nothing
end
