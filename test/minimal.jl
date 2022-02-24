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

################################################################################
# Test set_position! and set_maximal_velocity!
################################################################################\
@testset "minimal to maximal: set_position!, set_velocity" begin
	mech = Dojo.get_mechanism(:raiberthopper)
	timestep = mech.timestep
	joint1 = mech.joints[1]
	joint2 = mech.joints[2]
	pbody = mech.bodies[1]
	cbody = mech.bodies[2]
	tra2 = joint2.translational
	rot2 = joint2.rotational

	x = srand(1)
	Δx = Dojo.zerodimstaticadjoint(Dojo.nullspace_mask(tra2)) * x
	Δq = UnitQuaternion(rand(4)...)
	Dojo.set_maximal_configuration!(pbody, cbody; parent_vertex = tra2.vertices[1], child_vertex = tra2.vertices[2], Δx = Δx, Δq = Δq)
	@test norm(Dojo.minimal_coordinates(tra2, pbody, cbody) - x[1], Inf) < 1e-10

	v = srand(1)
	Δv = Dojo.zerodimstaticadjoint(Dojo.nullspace_mask(tra2)) * v
	Δω = rand(3)
	Dojo.set_maximal_velocity!(pbody, cbody; parent_vertex = tra2.vertices[1], child_vertex = tra2.vertices[2], Δv = Δv, Δω = Δω)
	@test norm(Dojo.minimal_velocities(tra2, pbody, cbody, timestep) - v[1], Inf) < 1e-10
end

################################################################################
# Test min -> max -> min
################################################################################

# raiberthopper
@testset "minimal to maximal to minimal: raibert hopper" begin
	mech = Dojo.get_mechanism(:raiberthopper);
	Random.seed!(100)
	nx = Dojo.minimal_dimension(mech)
	x0 = rand(nx)
	z0 = Dojo.minimal_to_maximal(mech, x0)
	x1 = Dojo.maximal_to_minimal(mech, z0)

	@test norm(x0[1:3] - x1[1:3], Inf) < 1e-10
	@test norm(x0[4:6] - x1[4:6]) < 1e-10
	@test norm(x0[7:9] - x1[7:9], Inf) < 1e-10
	@test norm(x0[11:12] - x1[11:12], Inf) < 1e-10
	@test norm(x0[13] - x1[13], Inf) < 1e-10
	@test norm(x0[14] - x1[14], Inf) < 1e-10
end

# box
@testset "minimal to maximal to minimal: box" begin
	mech = Dojo.get_mechanism(:box)
	Random.seed!(100)
	nx = Dojo.minimal_dimension(mech)
	x0 = rand(nx)
	z0 = Dojo.minimal_to_maximal(mech, x0)
	x1 = Dojo.maximal_to_minimal(mech, z0)
	@test norm(x0 - x1, Inf) < 1.0e-8
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
		x0 = rand(nx)
		z0 = Dojo.minimal_to_maximal(mech, x0)
		x1 = Dojo.maximal_to_minimal(mech, z0)
		@test norm(x0 - x1, Inf) < 1e-10
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
		x0 = rand(nx)
		z0 = Dojo.minimal_to_maximal(mech, x0)
		x1 = Dojo.maximal_to_minimal(mech, z0)
		@test norm(x0 - x1, Inf) < 1e-10
	end
end

# humanoid
@testset "minimal to maximal to minimal: humanoid" begin
	mech = Dojo.get_mechanism(:humanoid)
	Random.seed!(100)
	nx = Dojo.minimal_dimension(mech)
	x0 = rand(nx)
	z0 = Dojo.minimal_to_maximal(mech, x0)
	x1 = Dojo.maximal_to_minimal(mech, z0)
	@test norm(x0 - x1, Inf) < 1e-10
end

# quadruped
@testset "minimal to maximal to minimal: quadruped" begin
	mech = Dojo.get_mechanism(:quadruped)
	Random.seed!(100)
	nx = Dojo.minimal_dimension(mech)
	x0 = rand(nx)
	z0 = Dojo.minimal_to_maximal(mech, x0)
	x1 = Dojo.maximal_to_minimal(mech, z0)
	@test norm(x0 - x1, Inf) < 1e-10
end

# atlas
@testset "minimal to maximal to minimal: atlas" begin
	mech = Dojo.get_mechanism(:atlas, model_type = :simple, contact = true, damper = 10.0)
	Random.seed!(100)
	nx = Dojo.minimal_dimension(mech)
	x0 = rand(nx)
	z0 = Dojo.minimal_to_maximal(mech, x0)
	x1 = Dojo.maximal_to_minimal(mech, z0)
	@test norm(x0 - x1, Inf) < 1e-10
end

################################################################################
# Test set and get minimal coordinates and velocities
################################################################################
@testset "set and get minimal coordinates and velocities" begin
	for jointtype in jointtypes
	    mech = Dojo.get_snake(Nb=10, jointtype=jointtype)
		timestep = mech.timestep
		for joint in mech.joints
		    joint.rotational.axis_offset = UnitQuaternion(rand(4)...)
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
	        Dojo.set_minimal_coordinates!(rot0, pnodes0[i], cnodes0[i],  timestep, Δθ=Δθ)
	        Δθ0 = Dojo.minimal_coordinates(rot0, pnodes0[i], cnodes0[i])
	        @test norm(Δθ0 - Δθ, Inf) < 1e-7

	        Dojo.set_minimal_coordinates!(tra0, pnodes0[i], cnodes0[i], timestep, Δx=Δx)
	        Δx0 = Dojo.minimal_coordinates(tra0, pnodes0[i], cnodes0[i])
	        @test norm(Δx0 - Δx, Inf) < 1e-7

	        Dojo.set_minimal_velocities!(joint0, pnodes0[i], cnodes0[i], timestep, Δv=Δv, Δϕ=Δϕ)
	        Δϕ0 = Dojo.minimal_velocities(rot0, pnodes0[i], cnodes0[i], timestep)
			Δv0 = Dojo.minimal_velocities(tra0, pnodes0[i], cnodes0[i], timestep)
	        @test norm(Δϕ0 - Δϕ, Inf) < 1e-7
			@test norm(Δv0 - Δv, Inf) < 1e-7
	    end
	end
end


################################################################################
# Test minimal coordinates and velocities Jacobians
################################################################################
@testset "minimal velocity jacobian" begin
	mech = Dojo.get_humanoid()
	timestep = mech.timestep
	for jointcon in mech.joints
		for joint in [jointcon.translational, jointcon.rotational]
			# generate random configuration in minimal space
			x = rand(minimal_dimension(mech))

			# convert to maximal
			z = Dojo.minimal_to_maximal(mech, x)

			# extract body states
			Ne = Dojo.length(mech.joints)
			if Dojo.get_body(mech, jointcon.parent_id).name == :origin 
				zp = [mech.origin.state.x2[1]; mech.origin.state.v15; Dojo.vector(mech.origin.state.q2[1]); mech.origin.state.ϕ15]
			else
				zp = z[(jointcon.parent_id - Ne - 1) * 13 .+ (1:13)]
			end
			zc = z[(jointcon.child_id - Ne - 1) * 13 .+ (1:13)]

			xa = SVector{3}(zp[1:3])
			va = SVector{3}(zp[3 .+ (1:3)])
			qa = UnitQuaternion(zp[6 .+ (1:4)]..., false)
			ωa = SVector{3}(zp[10 .+ (1:3)])

			xb = SVector{3}(zc[1:3])
			vb = SVector{3}(zc[3 .+ (1:3)])
			qb = UnitQuaternion(zc[6 .+ (1:4)]..., false)
			ωb = SVector{3}(zc[10 .+ (1:3)])

			Dojo.minimal_velocities(joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep)

			# Jacobians
			∇0 = Dojo.minimal_velocities_jacobian_configuration(:parent, joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep)
			∇1 = FiniteDiff.finite_difference_jacobian(
				xq -> Dojo.minimal_velocities(joint, xq[Dojo.SUnitRange(1,3)], va, UnitQuaternion(xq[4:7]..., false), ωa, xb, vb, qb, ωb, timestep),
				[xa; Dojo.vector(qa)]) * cat(I(3), Dojo.LVᵀmat(qa), dims=(1,2))
			@test norm(∇0 - ∇1, Inf) < 1e-5

			∇0 = Dojo.minimal_velocities_jacobian_configuration(:child, joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep)
			∇1 = FiniteDiff.finite_difference_jacobian(
				xq -> Dojo.minimal_velocities(joint, xa, va, qa, ωa, xq[Dojo.SUnitRange(1,3)], vb, UnitQuaternion(xq[4:7]..., false), ωb, timestep),
				[xb; Dojo.vector(qb)]) * cat(I(3), Dojo.LVᵀmat(qb), dims=(1,2))
			@test norm(∇0 - ∇1, Inf) < 1e-5

			∇0 = Dojo.minimal_velocities_jacobian_velocity(:parent, joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep)
			∇1 = FiniteDiff.finite_difference_jacobian(
				vϕ -> Dojo.minimal_velocities(joint, xa, vϕ[Dojo.SUnitRange(1,3)], qa, vϕ[Dojo.SUnitRange(4,6)], xb, vb, qb, ωb, timestep),
				[va; ωa])
			@test norm(∇0 - ∇1, Inf) < 1e-5

			∇0 = Dojo.minimal_velocities_jacobian_velocity(:child, joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep)
			∇1 = FiniteDiff.finite_difference_jacobian(
				vϕ -> Dojo.minimal_velocities(joint, xa, va, qa, ωa, xb, vϕ[Dojo.SUnitRange(1,3)], qb, vϕ[Dojo.SUnitRange(4,6)], timestep),
				[vb; ωb])
			@test norm(∇0 - ∇1, Inf) < 1e-5
		end
	end
end

@testset "minimal coordinates jacobian" begin
	mech = Dojo.get_humanoid()
	for jointcon in mech.joints
		for joint in [jointcon.translational, jointcon.rotational]
			# generate random configuration in minimal space
			x = rand(Dojo.minimal_dimension(mech))

			# convert to maximal
			z = Dojo.minimal_to_maximal(mech, x)

			# extract body states
			Ne = Dojo.length(mech.joints)
			if Dojo.get_body(mech, jointcon.parent_id).name == :origin 
				zp = [mech.origin.state.x2[1]; mech.origin.state.v15; Dojo.vector(mech.origin.state.q2[1]); mech.origin.state.ϕ15]
			else
				zp = z[(jointcon.parent_id - Ne - 1) * 13 .+ (1:13)]
			end
			zc = z[(jointcon.child_id - Ne - 1) * 13 .+ (1:13)]

			xa = SVector{3}(zp[1:3])
			# va = SVector{3}(zp[3 .+ (1:3)])
			qa = UnitQuaternion(zp[6 .+ (1:4)]..., false)
			# ωa = SVector{3}(zp[10 .+ (1:3)])

			xb = SVector{3}(zc[1:3])
			# vb = SVector{3}(zc[3 .+ (1:3)])
			qb = UnitQuaternion(zc[6 .+ (1:4)]..., false)
			# ωb = SVector{3}(zc[10 .+ (1:3)])

			Dojo.minimal_coordinates(joint, xa, qa, xb, qb)

			∇0 = Dojo.minimal_coordinates_jacobian_configuration(:parent, joint, xa, qa, xb, qb)
			∇1 = FiniteDiff.finite_difference_jacobian(
				xq -> Dojo.minimal_coordinates(joint, xq[1:3], UnitQuaternion(xq[4:7]..., false), xb, qb),
				[xa; Dojo.vector(qa)]) * cat(I(3), Dojo.LVᵀmat(qa), dims=(1,2))
			@test norm(∇0 - ∇1, Inf) < 1e-6

			∇0 = Dojo.minimal_coordinates_jacobian_configuration(:child, joint, xa, qa, xb, qb)
			∇1 = FiniteDiff.finite_difference_jacobian(
				xq -> Dojo.minimal_coordinates(joint, xa, qa, xq[1:3], UnitQuaternion(xq[4:7]..., false)),
				[xb; Dojo.vector(qb)]) * cat(I(3), Dojo.LVᵀmat(qb), dims=(1,2))
			@test norm(∇0 - ∇1, Inf) < 1e-6
		end
	end
end

@testset "minimal-maximal Jacobians" begin
	function maximal_to_minimal_jacobian_fd(mechanism::Mechanism, z)
		J = FiniteDiff.finite_difference_jacobian(y -> Dojo.maximal_to_minimal(mechanism, y), z)
		G = attitude_jacobian(z, length(mechanism.bodies))
		return J * G
	end

	function minimal_to_maximal_jacobian_fd(mechanism::Mechanism, x)
		J = FiniteDiff.finite_difference_jacobian(y -> Dojo.minimal_to_maximal(mechanism, y), x)
		z = minimal_to_maximal(mechanism, x)
		G = attitude_jacobian(z, length(mechanism.bodies))
		return G'*J
	end

	function ctrl!(mechanism, k)
		Dojo.set_control!(mechanism, 0.1 * srand(Dojo.control_dimension(mechanism)))
	end

	# n-pendulum
	mechanism = Dojo.get_mechanism(:npendulum, timestep = 0.1, gravity = -9.81, Nb=1)
	ϕ1 = 0.3 * π
	Dojo.initialize!(mechanism, :npendulum, ϕ1 = ϕ1)
	storage = Dojo.simulate!(mechanism, 1.0, record = true, verbose = false)

	x = Dojo.get_minimal_state(mechanism)
	z = Dojo.get_maximal_state(mechanism)
	u = zeros(Dojo.control_dimension(mechanism))

	@test norm(minimal_to_maximal(mechanism, x) - z) < 1.0e-5 # NOTE: this won't necessarily pass
	@test norm(Dojo.maximal_to_minimal(mechanism, z) - x) < 1.0e-5

	M_fd = maximal_to_minimal_jacobian_fd(mechanism, z)
	M_a = Dojo.maximal_to_minimal_jacobian(mechanism, z)
	@test size(M_fd) == size(M_a)
	@test norm(M_fd - M_a, Inf) < 1.0e-5

	N_fd = minimal_to_maximal_jacobian_fd(mechanism, Dojo.maximal_to_minimal(mechanism, z))
	N_a = Dojo.minimal_to_maximal_jacobian(mechanism, Dojo.maximal_to_minimal(mechanism, z))
	@test size(N_fd) == size(N_a)
	@test norm(N_fd - N_a, Inf) < 1.0e-5
	@test abs(sum(diag(M_fd * N_fd)) - Dojo.minimal_dimension(mechanism)) < 1.0e-5
	@test abs(sum(diag(M_a * N_a)) - Dojo.minimal_dimension(mechanism)) < 1.0e-5

	# n-pendulum
	mechanism = Dojo.get_mechanism(:npendulum, timestep = 0.1, gravity = -9.81, Nb=2)
	ϕ1 = 0.3 * π
	Dojo.initialize!(mechanism, :npendulum, ϕ1 = ϕ1)
	storage = Dojo.simulate!(mechanism, 1.0, record = true, verbose = false)

	x = Dojo.get_minimal_state(mechanism)
	z = Dojo.get_maximal_state(mechanism)
	u = zeros(Dojo.control_dimension(mechanism))

	@test norm(minimal_to_maximal(mechanism, x) - z) < 1.0e-5 # NOTE: this won't necessarily pass
	@test norm(Dojo.maximal_to_minimal(mechanism, z) - x) < 1.0e-5

	M_fd = maximal_to_minimal_jacobian_fd(mechanism, z)
	M_a = Dojo.maximal_to_minimal_jacobian(mechanism, z)
	@test size(M_fd) == size(M_a)
	@test norm(M_fd - M_a, Inf) < 1.0e-5

	N_fd = minimal_to_maximal_jacobian_fd(mechanism, Dojo.maximal_to_minimal(mechanism, z))
	N_a = Dojo.minimal_to_maximal_jacobian(mechanism, Dojo.maximal_to_minimal(mechanism, z))
	@test size(N_fd) == size(N_a)
	@test norm(N_fd - N_a, Inf) < 1.0e-5
	@test abs(sum(diag(M_fd * N_fd)) - Dojo.minimal_dimension(mechanism)) < 1.0e-5
	@test abs(sum(diag(M_a * N_a)) - Dojo.minimal_dimension(mechanism)) < 1.0e-5

	# sphere
	mechanism = Dojo.get_mechanism(:sphere, timestep = 0.01, gravity = -9.81)
	Dojo.initialize!(mechanism, :sphere)
	storage = Dojo.simulate!(mechanism, 1.0, record = true, verbose = false)

	z = Dojo.get_maximal_state(mechanism)
	x = Dojo.get_minimal_state(mechanism)
	u = zeros(Dojo.control_dimension(mechanism))

	@test norm(minimal_to_maximal(mechanism, x) - z) < 1.0e-5 # NOTE: this won't necessarily pass
	@test norm(Dojo.maximal_to_minimal(mechanism, z) - x) < 1.0e-5

	M_fd = maximal_to_minimal_jacobian_fd(mechanism, z)
	M_a = Dojo.maximal_to_minimal_jacobian(mechanism, z)
	@test size(M_fd) == size(M_a)
	@test norm(M_fd - M_a, Inf) < 1.0e-5

	N_fd = minimal_to_maximal_jacobian_fd(mechanism, Dojo.maximal_to_minimal(mechanism, z))
	N_a = Dojo.minimal_to_maximal_jacobian(mechanism, Dojo.maximal_to_minimal(mechanism, z))
	@test size(N_fd) == size(N_a)
	@test norm(N_fd - N_a, Inf) < 1.0e-6
	@test abs(sum(diag(M_fd * N_fd)) - Dojo.minimal_dimension(mechanism)) < 1.0e-5
	@test abs(sum(diag(M_a * N_a)) - Dojo.minimal_dimension(mechanism)) < 1.0e-5

	# half cheetah
	mechanism = Dojo.get_mechanism(:halfcheetah, timestep=0.01, gravity=-9.81)
	Dojo.initialize!(mechanism, :halfcheetah)
	storage = Dojo.simulate!(mechanism, 1.0, ctrl!, record = true, verbose = false)

	z = Dojo.get_maximal_state(mechanism)
	x = Dojo.get_minimal_state(mechanism)
	u = zeros(Dojo.control_dimension(mechanism))

	@test norm(minimal_to_maximal(mechanism, x) - z) < 1.0e-5 # NOTE: this won't necessarily pass
	@test norm(Dojo.maximal_to_minimal(mechanism, z) - x) < 1.0e-5

	M_fd = maximal_to_minimal_jacobian_fd(mechanism, z)
	M_a = Dojo.maximal_to_minimal_jacobian(mechanism, z)
	@test size(M_fd) == size(M_a)
	@test norm(M_fd - M_a, Inf) < 1.0e-5

	N_fd = minimal_to_maximal_jacobian_fd(mechanism, Dojo.maximal_to_minimal(mechanism, z))
	N_a = Dojo.minimal_to_maximal_jacobian(mechanism, Dojo.maximal_to_minimal(mechanism, z))
	@test size(N_fd) == size(N_a)
	@test norm(N_fd - N_a, Inf) < 1.0e-5

	@test abs(sum(diag(M_fd * N_fd)) - Dojo.minimal_dimension(mechanism)) < 1.0e-5
	@test abs(sum(diag(M_a * N_a)) - Dojo.minimal_dimension(mechanism)) < 1.0e-5
	@test abs(sum(diag(M_a * N_fd)) - Dojo.minimal_dimension(mechanism)) < 1.0e-5

	# atlas
	mechanism = Dojo.get_mechanism(:atlas, timestep=0.01, gravity=-9.81, friction_coefficient=0.5, damper=100.0, spring=1.0, contact=true)
	Dojo.initialize_atlasstance!(mechanism, tran=[0,0,0.5], rot=[0.0,0.0,0.0])
	storage = Dojo.simulate!(mechanism, 1.0, ctrl!, record = true, verbose = false)

	z = Dojo.get_maximal_state(mechanism)
	x = Dojo.get_minimal_state(mechanism)
	u = zeros(Dojo.control_dimension(mechanism))

	@test norm(minimal_to_maximal(mechanism, x) - z) < 1.0e-5 # NOTE: this won't necessarily pass
	@test norm(Dojo.maximal_to_minimal(mechanism, z) - x) < 1.0e-5

	M_fd = maximal_to_minimal_jacobian_fd(mechanism, z)
	M_a = Dojo.maximal_to_minimal_jacobian(mechanism, z)
	@test size(M_fd) == size(M_a)
	@test norm(M_fd - M_a, Inf) < 1.0e-5

	N_fd = minimal_to_maximal_jacobian_fd(mechanism, Dojo.maximal_to_minimal(mechanism, z))
	N_a = Dojo.minimal_to_maximal_jacobian(mechanism, Dojo.maximal_to_minimal(mechanism, z))
	@test size(N_fd) == size(N_a)
	@test norm(N_fd - N_a, Inf) < 5.0e-5

	@test abs(sum(diag(M_fd * N_fd)) - Dojo.minimal_dimension(mechanism)) < 1.0e-5
	@test abs(sum(diag(M_a * N_a)) - Dojo.minimal_dimension(mechanism)) < 1.0e-5
end

@testset "maximal_to_minimal_jacobian" begin
	function maximal_to_minimal_jacobian_fd(mechanism::Mechanism, z)
		J = FiniteDiff.finite_difference_jacobian(y -> maximal_to_minimal(mechanism, y), z)
		G = attitude_jacobian(z, length(mechanism.bodies))
		return J * G
	end

	# 5-link pendulum
	mech = Dojo.get_mechanism(:npendulum, timestep = 0.01, gravity = -9.81, Nb=5)
	Random.seed!(100)
	ϕ1 = 0.3π
	Dojo.initialize!(mech, :npendulum, ϕ1 = ϕ1)
	storage = Dojo.simulate!(mech, 1.0, record = true, verbose = false)

	Dojo.maximal_dimension(mech) == 13
	Dojo.minimal_dimension(mech) == 12
	z = Dojo.get_maximal_state(mech)

	attjac = Dojo.attitude_jacobian(z, length(mech.bodies))
	M_fd = maximal_to_minimal_jacobian_fd(mech, z)
	M_a = Dojo.maximal_to_minimal_jacobian(mech, z)
	@test size(M_fd) == size(M_a)
	@test norm(M_fd - M_a, Inf) < 1.0e-6

	# # sphere
	mech = Dojo.get_mechanism(:sphere, timestep = 0.01, gravity = -9.81)
	Dojo.initialize!(mech, :sphere)
	storage = Dojo.simulate!(mech, 1.0, record = true, verbose = false)

	Dojo.maximal_dimension(mech)
	Dojo.minimal_dimension(mech)
	z = Dojo.get_maximal_state(mech)

	attjac = Dojo.attitude_jacobian(z, length(mech.bodies))
	M_fd = maximal_to_minimal_jacobian_fd(mech, z)
	M_a = Dojo.maximal_to_minimal_jacobian(mech, z)

	@test size(M_fd) == size(M_a)
	@test norm(M_fd - M_a, Inf) < 1.0e-6
end
