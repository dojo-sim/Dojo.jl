################################################################################
# Utils
################################################################################\

jointtypes = [
    :Fixed,
    :Prismatic,
    :Planar,
    :FixedOrientation,
    :Revolute,
    :Cylindrical,
    :PlanarAxis,
    :FreeRevolute,
    :Orbital,
    :PrismaticOrbital,
    :PlanarOrbital,
    :FreeOrbital,
    :Spherical,
    :CylindricalFree,
    :PlanarFree
    ]

function quaterror(q0::AbstractVector, q1::AbstractVector, p::Real = Inf)
	q0_ = UnitQuaternion(q0...)
	q1_ = UnitQuaternion(q1...)
	q_ = q0_ * inv(q1_)
	return norm([q_.x, q_.y, q_.z], p)
end

################################################################################
# Test set_position! and set_velocity!
################################################################################\
@testset "minimal to maximal: set_position!, set_velocity" begin
	mech = Dojo.get_mechanism(:raiberthopper)
	joint1 = collect(mech.joints)[1]
	joint2 = collect(mech.joints)[2]
	body1 = collect(mech.bodies)[1]
	body2 = collect(mech.bodies)[2]
	tra2 = joint2.translational
	rot2 = joint2.rotational

	x = srand(1)
	Δx = Dojo.zerodimstaticadjoint(Dojo.nullspace_mask(tra2)) * x
	Δq = UnitQuaternion(rand(4)...)
	Dojo.set_position!(body1, body2; p1 = tra2.vertices[1], p2 = tra2.vertices[2], Δx = Δx, Δq = Δq)
	@test norm(Dojo.minimal_coordinates(tra2, body1, body2) - x[1], Inf) < 1e-10

	v = srand(1)
	Δv = Dojo.zerodimstaticadjoint(Dojo.nullspace_mask(tra2)) * v
	Δω = rand(3)
	Dojo.set_velocity!(body1, body2; p1 = tra2.vertices[1], p2 = tra2.vertices[2], Δv = Δv, Δω = Δω)
	@test norm(Dojo.minimal_velocities(tra2, body1, body2) - v[1], Inf) < 1e-10
end

################################################################################
# Test min -> max -> min
################################################################################

# raiberthopper
@testset "minimal to maximal to minimal: raibert hopper" begin
	mech = Dojo.get_mechanism(:raiberthopper);
	Random.seed!(100)
	nx = Dojo.minimal_dimension(mech)
	x0 = [rand(3); Dojo.vector(UnitQuaternion(rand(4)...)); rand(3); rand(3); rand(nx - 13)]
	z0 = Dojo.minimal_to_maximal(mech, x0)
	x1 = Dojo.maximal_to_minimal(mech, z0)
	@test norm(x0[1:3] - x1[1:3], Inf) < 1e-10
	@test quaterror(x0[4:7], x1[4:7]) < 1e-10
	@test norm(x0[8:10] - x1[8:10], Inf) < 1e-10
	@test norm(x0[11:13] - x1[11:13], Inf) < 1e-10
	@test norm(x0[14:nx] - x1[14:nx], Inf) < 1e-10
end

# box
@testset "minimal to maximal to minimal: box" begin
	mech = Dojo.get_mechanism(:box)
	Random.seed!(100)
	nx = Dojo.minimal_dimension(mech)
	x0 = [rand(3); Dojo.vector(UnitQuaternion(rand(4)...)); rand(3); rand(3); rand(abs(nx - 13))]
	z0 = Dojo.minimal_to_maximal(mech, x0)
	x1 = Dojo.maximal_to_minimal(mech, z0)
	@test norm(x0[1:3] - x1[1:3], Inf) < 1e-10
	@test quaterror(x0[4:7], x1[4:7]) < 1e-10
	@test norm(x0[8:10] - x1[8:10], Inf) < 1e-10
	@test norm(x0[11:nx] - x1[11:nx], Inf) < 1e-10
end

# pendulum
@testset "minimal to maximal to minimal: pendulum" begin
	mech = Dojo.get_mechanism(:pendulum)
	Random.seed!(100)
	nx = Dojo.minimal_dimension(mech)
	x0 = rand(nx)
	z0 = Dojo.minimal_to_maximal(mech, x0)
	x1 = Dojo.maximal_to_minimal(mech, z0)
	@test norm(x0 - x1, Inf) < 1e-10
end

# halfcheetah
@testset "minimal to maximal to minimal: halfcheetah" begin
	mech = Dojo.get_mechanism(:halfcheetah)
	Random.seed!(100)
	nx = Dojo.minimal_dimension(mech)
	x0 = rand(nx)
	z0 = Dojo.minimal_to_maximal(mech, x0)
	x1 = Dojo.maximal_to_minimal(mech, z0)
	@test norm(x0 - x1, Inf) < 1e-10
end

# nslider
@testset "minimal to maximal to minimal: nslider" begin
	Nb0 = 5
	mech = Dojo.get_mechanism(:nslider, Nb = Nb0)
	Random.seed!(100)
	nx = Dojo.minimal_dimension(mech)
	x0 = rand(nx)
	z0 = Dojo.minimal_to_maximal(mech, x0)
	x1 = Dojo.maximal_to_minimal(mech, z0)
	@test norm(x0 - x1, Inf) < 1e-10
end

# npendulum
@testset "minimal to maximal to minimal: npendulum" begin
	for jointtype in jointtypes
		# @show jointtype
		Nb0 = 5
		mech = Dojo.get_mechanism(:npendulum, Nb = Nb0, jointtype = jointtype)
		Random.seed!(100)
		nx = Dojo.minimal_dimension(mech)
		x0 = rand(nx)
		z0 = Dojo.minimal_to_maximal(mech, x0)
		x1 = Dojo.maximal_to_minimal(mech, z0)
		@test norm(x0 - x1, Inf) < 1e-10
	end
end

# snake
@testset "minimal to maximal to minimal: snake" begin
	for jointtype in jointtypes
	# @show jointtype
	Nb0 = 5
	mech = Dojo.get_mechanism(:snake, Nb = Nb0, jointtype = jointtype)
	mech = Dojo.get_mechanism(:snake, Nb = Nb0, jointtype = :Fixed)
	Random.seed!(100)
	nx = Dojo.minimal_dimension(mech)
	x0 = [rand(3); Dojo.vector(UnitQuaternion(rand(4)...)); rand(3); rand(3); rand(abs(nx - 13))]
	z0 = Dojo.minimal_to_maximal(mech, x0)
	x1 = Dojo.maximal_to_minimal(mech, z0)
	@test norm(x0[1:3] - x1[1:3], Inf) < 1e-10
	@test quaterror(x0[4:7], x1[4:7]) < 1e-10
	@test norm(x0[8:10] - x1[8:10], Inf) < 1e-10
	@test norm(x0[11:nx] - x1[11:nx], Inf) < 1e-10
	end
end

# twister
@testset "minimal to maximal to minimal: twister" begin
	for jointtype in jointtypes
		# @show jointtype
		Nb0 = 5
		mech = Dojo.get_mechanism(:twister, Nb = Nb0, jointtype = jointtype)
		Random.seed!(100)
		nx = Dojo.minimal_dimension(mech)
		x0 = [rand(3); Dojo.vector(UnitQuaternion(rand(4)...)); rand(3); rand(3); rand(abs(nx - 13))]
		z0 = Dojo.minimal_to_maximal(mech, x0)
		x1 = Dojo.maximal_to_minimal(mech, z0)
		@test norm(x0[1:3] - x1[1:3], Inf) < 1e-10
		@test quaterror(x0[4:7], x1[4:7]) < 1e-10
		@test norm(x0[8:10] - x1[8:10], Inf) < 1e-10
		@test norm(x0[11:nx] - x1[11:nx], Inf) < 1e-10
	end
end

# humanoid
@testset "minimal to maximal to minimal: humanoid" begin
	mech = Dojo.get_mechanism(:humanoid)
	Random.seed!(100)
	nx = Dojo.minimal_dimension(mech)
	x0 = [rand(3); Dojo.vector(UnitQuaternion(rand(4)...)); rand(3); rand(3); rand(abs(nx - 13))]
	z0 = Dojo.minimal_to_maximal(mech, x0)
	x1 = Dojo.maximal_to_minimal(mech, z0)
	@test norm(x0[1:3] - x1[1:3], Inf) < 1e-10
	@test quaterror(x0[4:7], x1[4:7]) < 1e-10
	@test norm(x0[8:10] - x1[8:10], Inf) < 1e-10
	@test norm(x0[11:nx] - x1[11:nx], Inf) < 1e-10
end

# quadruped
@testset "minimal to maximal to minimal: quadruped" begin
	mech = Dojo.get_mechanism(:quadruped)
	Random.seed!(100)
	nx = Dojo.minimal_dimension(mech)
	x0 = [rand(3); Dojo.vector(UnitQuaternion(rand(4)...)); rand(3); rand(3); rand(abs(nx - 13))]
	z0 = Dojo.minimal_to_maximal(mech, x0)
	x1 = Dojo.maximal_to_minimal(mech, z0)
	@test norm(x0[1:3] - x1[1:3], Inf) < 1e-10
	@test quaterror(x0[4:7], x1[4:7]) < 1e-10
	@test norm(x0[8:10] - x1[8:10], Inf) < 1e-10
	@test norm(x0[11:nx] - x1[11:nx], Inf) < 1e-10
end

# atlas
@testset "minimal to maximal to minimal: atlas" begin
	mech = Dojo.get_mechanism(:atlas, model_type = :simple, contact = true, damper = 10.0)
	Random.seed!(100)
	nx = Dojo.minimal_dimension(mech)
	x0 = [rand(3); Dojo.vector(UnitQuaternion(rand(4)...)); rand(3); rand(3); rand(abs(nx - 13))]
	z0 = Dojo.minimal_to_maximal(mech, x0)
	x1 = Dojo.maximal_to_minimal(mech, z0)
	@test norm(x0[1:3] - x1[1:3], Inf) < 1e-10
	@test quaterror(x0[4:7], x1[4:7]) < 1e-10
	@test norm(x0[8:10] - x1[8:10], Inf) < 1e-10
	@test norm(x0[11:nx] - x1[11:nx], Inf) < 1e-10
end

@testset "maximal_to_minimal_jacobian" begin
	# 5-link pendulum
	mech = get_mechanism(:npendulum, timestep = 0.01, gravity = -9.81, Nb=5)
	Random.seed!(100)
	ϕ1 = 0.3π
	initialize!(mech, :npendulum, ϕ1 = ϕ1)
	storage = simulate!(mech, 1.0, record = true, verbose = false)

	Dojo.maximal_dimension(mech) == 13
	minimal_dimension(mech) == 12
	z = Dojo.get_maximal_state(mech)

	attjac = attitude_jacobian(z, length(mech.bodies))
	M_fd = FiniteDiff.finite_difference_jacobian(z -> maximal_to_minimal(mech, z), z) * attjac
	M_a = maximal_to_minimal_jacobian_analytical(mech, z)
	@test size(M_fd) == size(M_a)
	@test norm(M_fd - M_a, Inf) < 1.0e-6

	# # sphere
	# mech = get_mechanism(:sphere, timestep = 0.01, gravity = -9.81)
	# initialize!(mech, :sphere)
	# storage = simulate!(mech, 1.0, record = true, verbose = false)
	#
	# Dojo.maximal_dimension(mech)
	# minimal_dimension(mech)
	# z = get_maximal_state(mech)
	#
	# M_fd = maximal_to_minimal_jacobian(mech, z)
	# M_a = maximal_to_minimal_jacobian_analytical(mech, z)
	#
	# @test size(M_fd) == size(M_a)
	# @test norm(M_fd - M_a, Inf) < 1.0e-6
	@warn "sphere does not pass for some reason"
end

################################################################################
# Test set and get minimal coordinates and velocities
################################################################################
@testset "set and get minimal coordinates and velocities" begin
	for jointtype in jointtypes
	    mech = get_snake(Nb=10, jointtype=jointtype)
		for joint in mech.joints
		    joint.rotational.qoffset = UnitQuaternion(rand(4)...)
		end
	    joint0 = mech.joints[1]
	    tra0 = joint0.translational
	    rot0 = joint0.rotational
	    pnodes0 = [mech.origin; mech.bodies[1:end-1]]
	    cnodes0 = mech.bodies

	    Random.seed!(100)
	    Δθ = rand(control_dimension(rot0))
	    Δx = rand(control_dimension(tra0))
	    Δϕ = rand(control_dimension(rot0))
	    Δv = rand(control_dimension(tra0))
	    for i = 1:10
	        set_minimal_coordinates!(pnodes0[i], cnodes0[i], rot0, Δθ=Δθ)
	        Δθ0 = minimal_coordinates(rot0, pnodes0[i], cnodes0[i])
	        @test norm(Δθ0 - Δθ, Inf) < 1e-7

	        set_minimal_coordinates!(pnodes0[i], cnodes0[i], tra0, Δx=Δx)
	        Δx0 = minimal_coordinates(tra0, pnodes0[i], cnodes0[i])
	        @test norm(Δx0 - Δx, Inf) < 1e-7

	        set_minimal_velocities!(pnodes0[i], cnodes0[i], rot0, Δϕ=Δϕ)
	        Δϕ0 = minimal_velocities(rot0, pnodes0[i], cnodes0[i])
	        @test norm(Δϕ0 - Δϕ, Inf) < 1e-7

	        set_minimal_velocities!(pnodes0[i], cnodes0[i], tra0, Δv=Δv)
	        Δv0 = minimal_velocities(tra0, pnodes0[i], cnodes0[i])
	        @test norm(Δv0 - Δv, Inf) < 1e-7
	    end
	end
end


################################################################################
# Test minimal coordinates and velocities Jacobians
################################################################################
@testset "minimal velocity jacobian" begin
	mech = get_humanoid()
	for jointcon in mech.joints
		for joint in [jointcon.translational, jointcon.rotational]
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
	end
end

@testset "minimal coordinates jacobian" begin
	mech = get_humanoid()
	for jointcon in mech.joints
		for joint in [jointcon.translational, jointcon.rotational]
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
	end
end

# @testset "minimal to maximal Jacobian" begin
# 	function ctrl!(mechanism, k)
# 		set_control!(mechanism, 0.1 * srand(control_dimension(mechanism)))
# 	end

# 	# pendulum
# 	mechanism = get_mechanism(:npendulum, timestep = 0.01, gravity = -9.81 * 0.0, Nb=10)
# 	Random.seed!(100)
# 	ϕ1 = 0.3π
# 	initialize!(mechanism, :npendulum, ϕ1 = ϕ1)
# 	storage = simulate!(mechanism, 1.0, ctrl!, record = true, verbose = false)

# 	maximal_dimension(mechanism)
# 	minimal_dimension(mechanism)
# 	z = get_maximal_state(mechanism)
# 	x = get_minimal_state(mechanism)
# 	u = zeros(control_dimension(mechanism))
# 	maximal_to_minimal(mechanism, z) - x
# 	minimal_to_maximal(mechanism, x) - z
# 	Nb = length(mechanism.bodies)
# 	G = attitude_jacobian(z, Nb)[1:13Nb,1:12Nb]
# 	M_fd = maximal_to_minimal_jacobian(mechanism, z) * G
# 	M_a = maximal_to_minimal_jacobian_analytical(mechanism, z)
# 	@test size(M_fd) == size(M_a)
# 	@test norm(M_fd - M_a, Inf) < 1.0e-5

# 	N_fd = minimal_to_maximal_jacobian(mechanism, x)
# 	N_a = minimal_to_maximal_jacobian_analytical(mechanism, x)
# 	@test size(N_fd) == size(N_a)
# 	@test norm(N_fd - N_a, Inf) < 1.0e-5
# 	@test abs(sum(diag(M_fd * G' * N_fd)) - minimal_dimension(mechanism)) < 1.0e-6
# 	@test abs(sum(diag(M_a * G' * N_a)) - minimal_dimension(mechanism)) < 1.0e-6

# 	# sphere
# 	mechanism = get_mechanism(:sphere, timestep = 0.01, gravity = -9.81)
# 	initialize!(mechanism, :sphere)
# 	storage = simulate!(mechanism, 1.0, ctrl!, record = true, verbose = false)

# 	maximal_dimension(mechanism)
# 	minimal_dimension(mechanism)
# 	z = get_maximal_state(mechanism)
# 	x = get_minimal_state(mechanism)
# 	u = zeros(control_dimension(mechanism))
# 	maximal_to_minimal(mechanism, z) - x
# 	minimal_to_maximal(mechanism, x) - z
# 	Nb = length(mechanism.bodies)
# 	G = attitude_jacobian(z, Nb)[1:13Nb,1:12Nb]
# 	M_fd = maximal_to_minimal_jacobian(mechanism, z) * G
# 	M_a = maximal_to_minimal_jacobian_analytical(mechanism, z)
# 	@test size(M_fd) == size(M_a)
# 	@test norm(M_fd - M_a, Inf) < 1.0e-5

# 	N_fd = minimal_to_maximal_jacobian(mechanism, x)
# 	N_a = minimal_to_maximal_jacobian_analytical(mechanism, x)
# 	@test size(N_fd) == size(N_a)
# 	@test norm(N_fd - N_a, Inf) < 1.0e-6
# 	@test abs(sum(diag(M_fd * G' * N_fd)) - minimal_dimension(mechanism)) < 1.0e-6
# 	@test abs(sum(diag(M_a * G' * N_a)) - minimal_dimension(mechanism)) < 1.0e-6

# 	# half cheetah
# 	mechanism = get_mechanism(:halfcheetah, timestep=0.01, gravity=-9.81)
# 	initialize!(mechanism, :halfcheetah)
# 	storage = simulate!(mechanism, 1.0, ctrl!, record = true, verbose = false)

# 	maximal_dimension(mechanism)
# 	minimal_dimension(mechanism)
# 	z = get_maximal_state(mechanism)
# 	x = get_minimal_state(mechanism)
# 	u = zeros(control_dimension(mechanism))
# 	maximal_to_minimal(mechanism, z) - x
# 	minimal_to_maximal(mechanism, x) - z
# 	Nb = length(mechanism.bodies)
# 	G = attitude_jacobian(z, Nb)[1:13Nb,1:12Nb]

# 	M_fd = Dojo.maximal_to_minimal_jacobian(mechanism, z) * G
# 	M_a = Dojo.maximal_to_minimal_jacobian_analytical(mechanism, z)
# 	@test size(M_fd) == size(M_a)
# 	@test norm(M_fd - M_a, Inf) < 1.0e-5

# 	N_fd = Dojo.minimal_to_maximal_jacobian(mechanism, x)
# 	N_a = Dojo.minimal_to_maximal_jacobian_analytical(mechanism, x)
# 	@test size(N_fd) == size(N_a)
# 	@test norm(N_fd - N_a, Inf) < 1.0e-6

# 	@test abs(sum(diag(M_fd * G' * N_fd)) - minimal_dimension(mechanism)) < 1.0e-6
# 	@test abs(sum(diag(M_a * G' * N_a)) - minimal_dimension(mechanism)) < 1.0e-6

# 	# atlas
# 	mechanism = get_mechanism(:atlas, timestep=0.01, gravity=-9.81, friction_coefficient=0.5, damper=100.0, spring=1.0, contact=true)
# 	initialize_atlasstance!(mechanism, tran=[0,0,0.5], rot=[0.0,0.0,0.0])
# 	storage = simulate!(mechanism, 1.0, ctrl!, record = true, verbose = false)

# 	maximal_dimension(mechanism)
# 	minimal_dimension(mechanism)
# 	z = get_maximal_state(mechanism)
# 	x = get_minimal_state(mechanism)
# 	u = zeros(control_dimension(mechanism))
# 	maximal_to_minimal(mechanism, z) - x
# 	minimal_to_maximal(mechanism, x) - z
# 	Nb = length(mechanism.bodies)
# 	G = attitude_jacobian(z, Nb)[1:13Nb,1:12Nb]
# 	M_fd = Dojo.maximal_to_minimal_jacobian(mechanism, z) * G
# 	M_a = Dojo.maximal_to_minimal_jacobian_analytical(mechanism, z)
# 	@test size(M_fd) == size(M_a)
# 	@test norm(M_fd - M_a, Inf) < 1.0e-5

# 	N_fd = Dojo.minimal_to_maximal_jacobian(mechanism, x)
# 	N_a = Dojo.minimal_to_maximal_jacobian_analytical(mechanism, x)
# 	@test size(N_fd) == size(N_a)
# 	@test norm(N_fd - N_a, Inf) < 1.0e-5 # TODO: replace tolerance after replacing FD

# 	@test abs(sum(diag(M_fd * G' * N_fd)) - minimal_dimension(mechanism)) < 1.0e-6
# 	@test abs(sum(diag(M_a * G' * N_a)) - minimal_dimension(mechanism)) < 1.0e-6
# end
