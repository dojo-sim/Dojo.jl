@inline function get_position_delta(joint::Translational, body1::Node, body2::Node, x::SVector)
    Δx = zerodimstaticadjoint(nullspace_mask(joint)) * x # in body1 frame
    return Δx
end

@inline function get_velocity_delta(joint::Translational, body1::Node, body2::Node, v::SVector)
    Δv = zerodimstaticadjoint(nullspace_mask(joint)) * v # in body1 frame
    return Δv
end

@inline function position_error(joint::Translational, xa::AbstractVector,
		qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion; rotate::Bool = true)
	# TODO remove rotate
    vertices = joint.vertices
    d = xb + vrotate(vertices[2], qb) - (xa + vrotate(vertices[1], qa)) # in the world frame
    rotate && (d = vrotate(d, inv(qa))) # in the a frame
    return d
end

@inline function minimal_coordinates(joint::Translational, body1::Node, body2::Node)
    statea = body1.state
    stateb = body2.state
    return minimal_coordinates(joint, statea.x2[1], statea.q2[1], stateb.x2[1], stateb.q2[1])
end

@inline function minimal_coordinates(joint::Translational, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
    return nullspace_mask(joint) * position_error(joint, xa, qa, xb, qb)
end

@inline function minimal_velocities(joint::Translational, body1::Node, body2::Node)
    statea = body1.state
    stateb = body2.state
    return minimal_velocities(joint, statea.x2[1], statea.v15, statea.q2[1], statea.ϕ15, stateb.x2[1], stateb.v15, stateb.q2[1], stateb.ϕ15)
end

@inline function minimal_velocities(joint::Translational, xa::AbstractVector,
        va::AbstractVector,  qa::UnitQuaternion, ωa::AbstractVector,
        xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ωb::AbstractVector)
	vertices = joint.vertices
    pbcb_w = vrotate(-vertices[2], qb)
    pbca_w = xa - (xb + vrotate(vertices[2], qb))
    # Δvw = V(pb,B/A)w - V(pa,A/A)w
    Δvw = vb + skew(pbcb_w) * vrotate(ωb, qb) - (va + skew(pbca_w) * vrotate(ωa, qa)) # in world frame
    Δv = vrotate(Δvw, inv(qa)) # in the a frame
    return nullspace_mask(joint) * Δv
end

function minimal_coordinates_jacobian_configuration(jacobian_relative::Symbol, joint::Translational,
	xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
	jacobian_relative==:parent && (return [FiniteDiff.finite_difference_jacobian(x -> minimal_coordinates(joint, x, qa, xb, qb), xa) FiniteDiff.finite_difference_jacobian(q -> minimal_coordinates(joint, xa, UnitQuaternion(q..., false), xb, qb), vector(qa))]) # * LVᵀmat(qa)])
	jacobian_relative==:child && (return [FiniteDiff.finite_difference_jacobian(x -> minimal_coordinates(joint, xa, qa, x, qb), xb) FiniteDiff.finite_difference_jacobian(q -> minimal_coordinates(joint, xa, qa, xb, UnitQuaternion(q..., false)), vector(qb))]) # * LVᵀmat(qb)])
end

function minimal_velocities_jacobian_configuration(jacobian_relative::Symbol,
	joint::Translational, xa::AbstractVector, va::AbstractVector,
	qa::UnitQuaternion, ωa::AbstractVector, xb::AbstractVector,
	vb::AbstractVector, qb::UnitQuaternion, ωb::AbstractVector)

	(jacobian_relative == :parent) && (return [FiniteDiff.finite_difference_jacobian(x -> minimal_velocities(joint, x, va, qa, ωa, xb, vb, qb, ωb), xa) FiniteDiff.finite_difference_jacobian(q -> minimal_velocities(joint, xa, va, UnitQuaternion(q..., false), ωa, xb, vb, qb, ωb), vector(qa))]) # * LVᵀmat(qa)])
	(jacobian_relative == :child) && (return [FiniteDiff.finite_difference_jacobian(x -> minimal_velocities(joint, xa, va, qb, ωa, x, vb, qb, ωb), xb) FiniteDiff.finite_difference_jacobian(q -> minimal_velocities(joint, xa, va, qa, ωa, xb, vb, UnitQuaternion(q..., false), ωb), vector(qb))]) # * LVᵀmat(qb)])
	return
end

function minimal_velocities_jacobian_velocity(jacobian_relative::Symbol,
	joint::Translational, xa::AbstractVector, va::AbstractVector,
	qa::UnitQuaternion, ωa::AbstractVector, xb::AbstractVector,
	vb::AbstractVector, qb::UnitQuaternion, ωb::AbstractVector)
	(jacobian_relative == :parent) && (return [FiniteDiff.finite_difference_jacobian(v -> minimal_velocities(joint, xa, v, qa, ωa, xb, vb, qb, ωb), va) FiniteDiff.finite_difference_jacobian(ω -> minimal_velocities(joint, xa, va, qa, ω, xb, vb, qb, ωb), ωa)])
	(jacobian_relative == :child) && (return [FiniteDiff.finite_difference_jacobian(v -> minimal_velocities(joint, xa, va, qa, ωa, xb, v, qb, ωb), vb) FiniteDiff.finite_difference_jacobian(ω -> minimal_velocities(joint, xa, va, qa, ωa, xb, vb, qb, ω), ωb)])
	return
end






mech = get_humanoid()
joint = mech.joints[1].constraints[1]
qa = UnitQuaternion(rand(4)...)
qb = UnitQuaternion(rand(4)...)
xa = rand(3)
va = rand(3)
ωa = rand(3)
xb = rand(3)
vb = rand(3)
ωb = rand(3)
minimal_velocities(joint, xa, va, qa, ωa, xb, vb, qb, ωb)

∇0 = minimal_velocities_jacobian_configuration_parent(joint, xa, va, qa, ωa, xb, vb, qb, ωb)
∇1 = FiniteDiff.finite_difference_jacobian(
	xq -> minimal_velocities(joint, xq[1:3], va, UnitQuaternion(xq[4:7]..., false), ωa, xb, vb, qb, ωb),
	[xa; vector(qa)]) * cat(I(3), LVᵀmat(qa), dims=(1,2))
@test norm(∇0 - ∇1, Inf) < 1e-6

∇0 = minimal_velocities_jacobian_configuration_child(joint, xa, va, qa, ωa, xb, vb, qb, ωb)
∇1 = FiniteDiff.finite_difference_jacobian(
	xq -> minimal_velocities(joint, xa, va, qa, ωa, xq[1:3], vb, UnitQuaternion(xq[4:7]..., false), ωb),
	[xb; vector(qb)]) * cat(I(3), LVᵀmat(qb), dims=(1,2))
@test norm(∇0 - ∇1, Inf) < 1e-6

∇0 = minimal_velocities_jacobian_velocity_parent(joint, xa, va, qa, ωa, xb, vb, qb, ωb)
∇1 = FiniteDiff.finite_difference_jacobian(
	vϕ -> minimal_velocities(joint, xa, vϕ[1:3], qa, vϕ[4:6], xb, vb, qb, ωb),
	[va; ωa])
@test norm(∇0 - ∇1, Inf) < 1e-6

∇0 = minimal_velocities_jacobian_velocity_child(joint, xa, va, qa, ωa, xb, vb, qb, ωb)
∇1 = FiniteDiff.finite_difference_jacobian(
	vϕ -> minimal_velocities(joint, xa, va, qa, ωa, xb, vϕ[1:3], qb, vϕ[4:6]),
	[vb; ωb])
@test norm(∇0 - ∇1, Inf) < 1e-6


@testset "minimal velocity jacobian" begin
    mech = getpendulum()
    joint = mech.joints[1].constraints[2]
    qa = UnitQuaternion(rand(4)...)
    qb = UnitQuaternion(rand(4)...)
    xa = rand(3)
    va = rand(3)
    ωa = rand(3)
    xb = rand(3)
    vb = rand(3)
    ωb = rand(3)
    minimal_velocities(joint, xa, va, qa, ωa, xb, vb, qb, ωb)

    ∇0 = minimal_velocities_jacobian_configuration_parent(joint, xa, va, qa, ωa, xb, vb, qb, ωb)
    ∇1 = FiniteDiff.finite_difference_jacobian(
        xq -> minimal_velocities(joint, xq[1:3], va, UnitQuaternion(xq[4:7]..., false), ωa, xb, vb, qb, ωb),
        [xa; vector(qa)]) * cat(I(3), LVᵀmat(qa), dims=(1,2))
    @test norm(∇0 - ∇1, Inf) < 1e-6

    ∇0 = minimal_velocities_jacobian_configuration_child(joint, xa, va, qa, ωa, xb, vb, qb, ωb)
    ∇1 = FiniteDiff.finite_difference_jacobian(
        xq -> minimal_velocities(joint, xa, va, qa, ωa, xq[1:3], vb, UnitQuaternion(xq[4:7]..., false), ωb),
        [xb; vector(qb)]) * cat(I(3), LVᵀmat(qb), dims=(1,2))
    @test norm(∇0 - ∇1, Inf) < 1e-6

    ∇0 = minimal_velocities_jacobian_velocity_parent(joint, xa, va, qa, ωa, xb, vb, qb, ωb)
    ∇1 = FiniteDiff.finite_difference_jacobian(
        vϕ -> minimal_velocities(joint, xa, vϕ[1:3], qa, vϕ[4:6], xb, vb, qb, ωb),
        [va; ωa])
    @test norm(∇0 - ∇1, Inf) < 1e-6

    ∇0 = minimal_velocities_jacobian_velocity_child(joint, xa, va, qa, ωa, xb, vb, qb, ωb)
    ∇1 = FiniteDiff.finite_difference_jacobian(
        vϕ -> minimal_velocities(joint, xa, va, qa, ωa, xb, vϕ[1:3], qb, vϕ[4:6]),
        [vb; ωb])
    @test norm(∇0 - ∇1, Inf) < 1e-6
end

@testset "minimal coordinates jacobian" begin
    mech = getpendulum()
    joint = mech.joints[1].constraints[2]
    qa = UnitQuaternion(rand(4)...)
    qb = UnitQuaternion(rand(4)...)
    xa = rand(3)
    va = rand(3)
    ωa = rand(3)
    xb = rand(3)
    vb = rand(3)
    ωb = rand(3)
    minimal_coordinates(joint, xa, qa, xb, qb)

    ∇0 = minimal_coordinates_jacobian_configuration(:parent, joint, xa, qa, xb, qb)
    ∇1 = FiniteDiff.finite_difference_jacobian(
        xq -> minimal_coordinates(joint, xq[1:3], UnitQuaternion(xq[4:7]..., false), xb, qb),
        [xa; vector(qa)]) * cat(I(3), LVᵀmat(qa), dims=(1,2))
    @test norm(∇0 - ∇1, Inf) < 1e-6

    ∇0 = minimal_coordinates_jacobian_configuration(:child, joint, xa, qa, xb, qb)
    ∇1 = FiniteDiff.finite_difference_jacobian(
        xq -> minimal_coordinates(joint, xa, qa, xq[1:3], UnitQuaternion(xq[4:7]..., false)),
        [xb; vector(qb)]) * cat(I(3), LVᵀmat(qb), dims=(1,2))
    @test norm(∇0 - ∇1, Inf) < 1e-6
end
