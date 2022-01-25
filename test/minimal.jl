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
# Test set_position and set_velocity!
################################################################################\
@testset "minMaxCoord: setPos!, setVel!" begin
	mech = Dojo.getmechanism(:raiberthopper)
	eqc1 = collect(mech.eqconstraints)[1]
	eqc2 = collect(mech.eqconstraints)[2]
	body1 = collect(mech.bodies)[1]
	body2 = collect(mech.bodies)[2]
	tra2 = eqc2.constraints[1]
	rot2 = eqc2.constraints[1]

	x = srand(1)
	Δx = Dojo.zerodimstaticadjoint(Dojo.nullspace_mask(tra2)) * x
	Δq = UnitQuaternion(rand(4)...)
	Dojo.set_position(body1, body2; p1 = tra2.vertices[1], p2 = tra2.vertices[2], Δx = Δx, Δq = Δq)
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
@testset "min -> max -> min: raiberthopper" begin
	mech = Dojo.getmechanism(:raiberthopper);
	Random.seed!(100)
	nx = Dojo.minimal_dimension(mech)
	x0 = [rand(3); Dojo.vector(UnitQuaternion(rand(4)...)); rand(3); rand(3); rand(nx - 13)]
	z0 = Dojo.min2max(mech, x0)
	x1 = Dojo.max2min(mech, z0)
	@test norm(x0[1:3] - x1[1:3], Inf) < 1e-10
	@test quaterror(x0[4:7], x1[4:7]) < 1e-10
	@test norm(x0[8:10] - x1[8:10], Inf) < 1e-10
	@test norm(x0[11:13] - x1[11:13], Inf) < 1e-10
	@test norm(x0[14:nx] - x1[14:nx], Inf) < 1e-10
end

# box
@testset "min -> max -> min: box" begin
	mech = Dojo.getmechanism(:box)
	Random.seed!(100)
	nx = Dojo.minimal_dimension(mech)
	x0 = [rand(3); Dojo.vector(UnitQuaternion(rand(4)...)); rand(3); rand(3); rand(abs(nx - 13))]
	z0 = Dojo.min2max(mech, x0)
	x1 = Dojo.max2min(mech, z0)
	@test norm(x0[1:3] - x1[1:3], Inf) < 1e-10
	@test quaterror(x0[4:7], x1[4:7]) < 1e-10
	@test norm(x0[8:10] - x1[8:10], Inf) < 1e-10
	@test norm(x0[11:nx] - x1[11:nx], Inf) < 1e-10
end

# pendulum
@testset "min -> max -> min: pendulum" begin
	mech = Dojo.getmechanism(:pendulum)
	Random.seed!(100)
	nx = Dojo.minimal_dimension(mech)
	x0 = rand(nx)
	z0 = Dojo.min2max(mech, x0)
	x1 = Dojo.max2min(mech, z0)
	@test norm(x0 - x1, Inf) < 1e-10
end

# halfcheetah
@testset "min -> max -> min: halfcheetah" begin
	mech = Dojo.getmechanism(:halfcheetah)
	Random.seed!(100)
	nx = Dojo.minimal_dimension(mech)
	x0 = rand(nx)
	z0 = Dojo.min2max(mech, x0)
	x1 = Dojo.max2min(mech, z0)
	@test norm(x0 - x1, Inf) < 1e-10
end

# nslider
@testset "min -> max -> min: nslider" begin
	Nb0 = 5
	mech = Dojo.getmechanism(:nslider, Nb = Nb0)
	Random.seed!(100)
	nx = Dojo.minimal_dimension(mech)
	x0 = rand(nx)
	z0 = Dojo.min2max(mech, x0)
	x1 = Dojo.max2min(mech, z0)
	@test norm(x0 - x1, Inf) < 1e-10
end

# npendulum
@testset "min -> max -> min: npendulum" begin
	for jointtype in jointtypes
		# @show jointtype
		Nb0 = 5
		mech = Dojo.getmechanism(:npendulum, Nb = Nb0, jointtype = jointtype)
		Random.seed!(100)
		nx = Dojo.minimal_dimension(mech)
		x0 = rand(nx)
		z0 = Dojo.min2max(mech, x0)
		x1 = Dojo.max2min(mech, z0)
		@test norm(x0 - x1, Inf) < 1e-10
	end
end

# snake
@testset "min -> max -> min: snake" begin
	for jointtype in jointtypes
	# @show jointtype
	Nb0 = 5
	mech = Dojo.getmechanism(:snake, Nb = Nb0, jointtype = jointtype)
	mech = Dojo.getmechanism(:snake, Nb = Nb0, jointtype = :Fixed)
	Random.seed!(100)
	nx = Dojo.minimal_dimension(mech)
	x0 = [rand(3); Dojo.vector(UnitQuaternion(rand(4)...)); rand(3); rand(3); rand(abs(nx - 13))]
	z0 = Dojo.min2max(mech, x0)
	x1 = Dojo.max2min(mech, z0)
	@test norm(x0[1:3] - x1[1:3], Inf) < 1e-10
	@test quaterror(x0[4:7], x1[4:7]) < 1e-10
	@test norm(x0[8:10] - x1[8:10], Inf) < 1e-10
	@test norm(x0[11:nx] - x1[11:nx], Inf) < 1e-10
	end
end

# twister
@testset "min -> max -> min: twister" begin
	for jointtype in jointtypes
		# @show jointtype
		Nb0 = 5
		mech = Dojo.getmechanism(:twister, Nb = Nb0, jointtype = jointtype)
		Random.seed!(100)
		nx = Dojo.minimal_dimension(mech)
		x0 = [rand(3); Dojo.vector(UnitQuaternion(rand(4)...)); rand(3); rand(3); rand(abs(nx - 13))]
		z0 = Dojo.min2max(mech, x0)
		x1 = Dojo.max2min(mech, z0)
		@test norm(x0[1:3] - x1[1:3], Inf) < 1e-10
		@test quaterror(x0[4:7], x1[4:7]) < 1e-10
		@test norm(x0[8:10] - x1[8:10], Inf) < 1e-10
		@test norm(x0[11:nx] - x1[11:nx], Inf) < 1e-10
	end
end

# humanoid
@testset "min -> max -> min: humanoid" begin
	mech = Dojo.getmechanism(:humanoid)
	Random.seed!(100)
	nx = Dojo.minimal_dimension(mech)
	x0 = [rand(3); Dojo.vector(UnitQuaternion(rand(4)...)); rand(3); rand(3); rand(abs(nx - 13))]
	z0 = Dojo.min2max(mech, x0)
	x1 = Dojo.max2min(mech, z0)
	@test norm(x0[1:3] - x1[1:3], Inf) < 1e-10
	@test quaterror(x0[4:7], x1[4:7]) < 1e-10
	@test norm(x0[8:10] - x1[8:10], Inf) < 1e-10
	@test norm(x0[11:nx] - x1[11:nx], Inf) < 1e-10
end

# quadruped
@testset "min -> max -> min: quadruped" begin
	mech = Dojo.getmechanism(:quadruped)
	Random.seed!(100)
	nx = Dojo.minimal_dimension(mech)
	x0 = [rand(3); Dojo.vector(UnitQuaternion(rand(4)...)); rand(3); rand(3); rand(abs(nx - 13))]
	z0 = Dojo.min2max(mech, x0)
	x1 = Dojo.max2min(mech, z0)
	@test norm(x0[1:3] - x1[1:3], Inf) < 1e-10
	@test quaterror(x0[4:7], x1[4:7]) < 1e-10
	@test norm(x0[8:10] - x1[8:10], Inf) < 1e-10
	@test norm(x0[11:nx] - x1[11:nx], Inf) < 1e-10
end

# atlas
@testset "min -> max -> min: atlas" begin
	mech = Dojo.getmechanism(:atlas, model_type = :simple, contact = true, damper = 10.0)
	Random.seed!(100)
	nx = Dojo.minimal_dimension(mech)
	x0 = [rand(3); Dojo.vector(UnitQuaternion(rand(4)...)); rand(3); rand(3); rand(abs(nx - 13))]
	z0 = Dojo.min2max(mech, x0)
	x1 = Dojo.max2min(mech, z0)
	@test norm(x0[1:3] - x1[1:3], Inf) < 1e-10
	@test quaterror(x0[4:7], x1[4:7]) < 1e-10
	@test norm(x0[8:10] - x1[8:10], Inf) < 1e-10
	@test norm(x0[11:nx] - x1[11:nx], Inf) < 1e-10
end
