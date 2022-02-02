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
	tra2 = joint2.constraints[1]
	rot2 = joint2.constraints[1]

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

	maximal_dimension(mech) == 13
	minimal_dimension(mech) == 12
	z = get_maximal_state(mech)

	M_fd = FiniteDiff.finite_difference_jacobian(z -> maximal_to_minimal(mech, z), z)
	M_a = maximal_to_minimal_jacobian_analytical(mech, z)
	@test size(M_fd) == size(M_a)
	@test norm(M_fd - M_a, Inf) < 1.0e-6

	# sphere
	mech = get_mechanism(:sphere, timestep = 0.01, gravity = -9.81)
	initialize!(mech, :sphere)
	storage = simulate!(mech, 1.0, record = true, verbose = false)

	maximal_dimension(mech)
	minimal_dimension(mech)
	z = get_maximal_state(mech)

	M_fd = maximal_to_minimal_jacobian(mech, z)
	M_a = maximal_to_minimal_jacobian_analytical(mech, z)

	@test size(M_fd) == size(M_a)
	@test norm(M_fd - M_a, Inf) < 1.0e-6
end