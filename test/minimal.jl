################################################################################
# Utils
################################################################################\

joint_types = [
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

# TODO this is necessary, because some conversions do not return a vector for scalar values
# See TODO below for removing finite diff
function force_to_jacobian_finite_diff(f,x)
	if typeof(f(x)) <: AbstractVector
		return FiniteDiff.finite_difference_jacobian(f,x)
	else
		return FiniteDiff.finite_difference_jacobian(x -> [f(x)],x)
	end
end


function force_to_jacobian_forward_diff(f,x)
	if typeof(f(x)) <: AbstractVector
		return ForwardDiff.jacobian(f,x)
	else
		return ForwardDiff.jacobian(x -> [f(x)],x)
	end
end

################################################################################
# Test get and set position and velocities
################################################################################
@testset "Get and set position and velocities" begin
	@testset "Maximal coordinates" begin
		mech = DojoEnvironments.get_mechanism(:raiberthopper)
		timestep= mech.timestep
		joint1 = mech.joints[1]
		joint2 = mech.joints[2]
		pbody = mech.bodies[1]
		cbody = mech.bodies[2]
		tra2 = joint2.translational
		rot2 = joint2.rotational

		x = srand(1)
		Δx = Dojo.zerodimstaticadjoint(Dojo.nullspace_mask(tra2)) * x
		Δq = rand(QuatRotation).q
		Dojo.set_maximal_configurations!(pbody, cbody;
			parent_vertex=tra2.vertices[1],
			child_vertex=tra2.vertices[2],
			Δx=Δx,
			Δq=Δq)
		@test norm(Dojo.minimal_coordinates(tra2, pbody, cbody) - x[1], Inf) < 1.0e-8

		v = srand(1)
		Δv = Dojo.zerodimstaticadjoint(Dojo.nullspace_mask(tra2)) * v
		Δω = rand(3)
		Dojo.set_maximal_velocities!(pbody, cbody;
			parent_vertex=tra2.vertices[1],
			child_vertex=tra2.vertices[2],
			Δv=Δv,
			Δω=Δω)
		@test norm(Dojo.minimal_velocities(tra2, pbody, cbody, timestep) - v[1], Inf) < 1.0e-8
	end

	@testset "Minimal coordinates" begin
		for joint_type in joint_types
			mech = DojoEnvironments.get_snake(
					num_bodies=10,
					joint_type=joint_type)

			timestep = mech.timestep
			for joint in mech.joints
				joint.rotational.axis_offset = rand(QuatRotation).q
			end
			joint0 = mech.joints[1]
			tra0 = joint0.translational
			rot0 = joint0.rotational
			pnodes0 = [mech.origin; mech.bodies[1:end-1]]
			cnodes0 = mech.bodies

			Random.seed!(100)
			Δθ = rand(input_dimension(rot0))
			Δx = rand(input_dimension(tra0))
			Δϕ = rand(input_dimension(rot0))
			Δv = rand(input_dimension(tra0))
			for i = 1:10
				Dojo.set_minimal_coordinates!(rot0, pnodes0[i], cnodes0[i],  timestep,
					Δθ=Δθ)
				Δθ0 = Dojo.minimal_coordinates(rot0, pnodes0[i], cnodes0[i])
				@test norm(Δθ0 - Δθ, Inf) < 1.0e-8

				Dojo.set_minimal_coordinates!(tra0, pnodes0[i], cnodes0[i], timestep,
					Δx=Δx)
				Δx0 = Dojo.minimal_coordinates(tra0, pnodes0[i], cnodes0[i])
				@test norm(Δx0 - Δx, Inf) < 1.0e-8

				Dojo.set_minimal_velocities!(joint0, pnodes0[i], cnodes0[i], timestep,
					Δv=Δv,
					Δϕ=Δϕ)
				Δϕ0 = Dojo.minimal_velocities(rot0, pnodes0[i], cnodes0[i], timestep)
				Δv0 = Dojo.minimal_velocities(tra0, pnodes0[i], cnodes0[i], timestep)
				@test norm(Δϕ0 - Δϕ, Inf) < 1.0e-8
				@test norm(Δv0 - Δv, Inf) < 1.0e-8
			end
		end
	end
end

################################################################################
# Test min -> max -> min
################################################################################
@testset "Minimal to maximal to minimal" begin
	# raiberthopper
	@testset "Raibert hopper" begin
		mech = DojoEnvironments.get_mechanism(:raiberthopper);
		Random.seed!(100)
		nx = Dojo.minimal_dimension(mech)
		x0 = rand(nx)
		z0 = Dojo.minimal_to_maximal(mech, x0)
		x1 = Dojo.maximal_to_minimal(mech, z0)

		@test norm(x0[1:3] - x1[1:3], Inf) < 1.0e-8
		@test norm(x0[4:6] - x1[4:6]) < 1.0e-8
		@test norm(x0[7:9] - x1[7:9], Inf) < 1.0e-8
		@test norm(x0[11:12] - x1[11:12], Inf) < 1.0e-8
		@test norm(x0[13] - x1[13], Inf) < 1.0e-8
		@test norm(x0[14] - x1[14], Inf) < 1.0e-8
	end

	# box
	@testset "Box" begin
		mech = DojoEnvironments.get_mechanism(:block)
		Random.seed!(100)
		nx = Dojo.minimal_dimension(mech)
		x0 = rand(nx)
		z0 = Dojo.minimal_to_maximal(mech, x0)
		x1 = Dojo.maximal_to_minimal(mech, z0)
		@test norm(x0 - x1, Inf) < 1.0e-8
	end

	# pendulum
	@testset "Pendulum" begin
		mech = DojoEnvironments.get_mechanism(:pendulum)
		Random.seed!(100)
		nx = Dojo.minimal_dimension(mech)
		x0 = rand(nx)
		z0 = Dojo.minimal_to_maximal(mech, x0)
		x1 = Dojo.maximal_to_minimal(mech, z0)
		@test norm(x0 - x1, Inf) < 1.0e-8
	end

	# halfcheetah
	@testset "Halfcheetah" begin
		mech = DojoEnvironments.get_mechanism(:halfcheetah)
		Random.seed!(100)
		nx = Dojo.minimal_dimension(mech)
		x0 = rand(nx)
		z0 = Dojo.minimal_to_maximal(mech, x0)
		x1 = Dojo.maximal_to_minimal(mech, z0)
		@test norm(x0 - x1, Inf) < 1.0e-8
	end

	# nslider
	@testset "Nslider" begin
		Nb0 = 5
		mech = DojoEnvironments.get_mechanism(:nslider,
			num_bodies=Nb0)
		Random.seed!(100)
		nx = Dojo.minimal_dimension(mech)
		x0 = rand(nx)
		z0 = Dojo.minimal_to_maximal(mech, x0)
		x1 = Dojo.maximal_to_minimal(mech, z0)
		@test norm(x0 - x1, Inf) < 1.0e-8
	end

	# npendulum
	@testset "Npendulum" begin
		for joint_type in joint_types
			# @show joint_type
			Nb0 = 5
			mech = DojoEnvironments.get_mechanism(:npendulum,
				num_bodies=Nb0,
				joint_type=joint_type)
			Random.seed!(100)
			nx = Dojo.minimal_dimension(mech)
			x0 = rand(nx)
			z0 = Dojo.minimal_to_maximal(mech, x0)
			x1 = Dojo.maximal_to_minimal(mech, z0)
			@test norm(x0 - x1, Inf) < 1.0e-8
		end
	end

	# snake
	@testset "Snake" begin
		for joint_type in joint_types
			# @show joint_type
			Nb0 = 5
			mech = DojoEnvironments.get_mechanism(:snake,
				num_bodies=Nb0,
				joint_type=joint_type)
			mech = DojoEnvironments.get_mechanism(:snake,
				num_bodies=Nb0,
				joint_type=:Fixed)
			Random.seed!(100)
			nx = Dojo.minimal_dimension(mech)
			x0 = rand(nx)
			z0 = Dojo.minimal_to_maximal(mech, x0)
			x1 = Dojo.maximal_to_minimal(mech, z0)
			@test norm(x0 - x1, Inf) < 1.0e-8
		end
	end

	# twister
	@testset "Twister" begin
		for joint_type in joint_types
			# @show joint_type
			Nb0 = 5
			mech = DojoEnvironments.get_mechanism(:twister,
				num_bodies=Nb0,
				joint_type=joint_type)
			Random.seed!(100)
			nx = Dojo.minimal_dimension(mech)
			x0 = rand(nx)
			z0 = Dojo.minimal_to_maximal(mech, x0)
			x1 = Dojo.maximal_to_minimal(mech, z0)
			@test norm(x0 - x1, Inf) < 1.0e-8
		end
	end

	# humanoid
	@testset "Humanoid" begin
		mech = DojoEnvironments.get_mechanism(:humanoid)
		Random.seed!(100)
		nx = Dojo.minimal_dimension(mech)
		x0 = rand(nx)
		z0 = Dojo.minimal_to_maximal(mech, x0)
		x1 = Dojo.maximal_to_minimal(mech, z0)
		@test norm(x0 - x1, Inf) < 1.0e-8
	end

	# quadruped
	@testset "Quadruped" begin
		mech = DojoEnvironments.get_mechanism(:quadruped)
		Random.seed!(100)
		nx = Dojo.minimal_dimension(mech)
		x0 = rand(nx)
		z0 = Dojo.minimal_to_maximal(mech, x0)
		x1 = Dojo.maximal_to_minimal(mech, z0)
		@test norm(x0 - x1, Inf) < 1.0e-8
	end

	# atlas
	@testset "Atlas" begin
		mech = DojoEnvironments.get_mechanism(:atlas,
			model_type=:simple,
			contact_feet=true,
			damper=10.0)

		Random.seed!(100)
		nx = Dojo.minimal_dimension(mech)
		x0 = rand(nx)
		z0 = Dojo.minimal_to_maximal(mech, x0)
		x1 = Dojo.maximal_to_minimal(mech, z0)
		@test norm(x0 - x1, Inf) < 1.0e-8
	end
end

################################################################################
#Test minimal coordinates and velocities Jacobians
################################################################################
@testset "Jacobians" begin
	@testset "Minimal velocity Jacobian" begin
		mech = DojoEnvironments.get_humanoid()
		timestep= mech.timestep
		for jointcon in mech.joints
			for joint in [jointcon.translational, jointcon.rotational]
				# generate random configuration in minimal space
				x = rand(minimal_dimension(mech))

				# convert to maximal
				z = Dojo.minimal_to_maximal(mech, x)

				# extract body states
				Ne = Dojo.length(mech.joints)
				if Dojo.get_body(mech, jointcon.parent_id).name == :origin
					zp = [mech.origin.state.x2; mech.origin.state.v15; Dojo.vector(mech.origin.state.q2); mech.origin.state.ϕ15]
				else
					zp = z[(jointcon.parent_id - Ne - 1) * 13 .+ (1:13)]
				end
				zc = z[(jointcon.child_id - Ne - 1) * 13 .+ (1:13)]

				xa = SVector{3}(zp[1:3])
				va = SVector{3}(zp[3 .+ (1:3)])
				qa = Quaternion(zp[6 .+ (1:4)]...)
				ωa = SVector{3}(zp[10 .+ (1:3)])

				xb = SVector{3}(zc[1:3])
				vb = SVector{3}(zc[3 .+ (1:3)])
				qb = Quaternion(zc[6 .+ (1:4)]...)
				ωb = SVector{3}(zc[10 .+ (1:3)])

				Dojo.minimal_velocities(joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep)

				# Jacobians
				∇0 = Dojo.minimal_velocities_jacobian_configuration(:parent, joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep)
				∇1 = force_to_jacobian_forward_diff(
					xq -> Dojo.minimal_velocities(joint, xq[Dojo.SUnitRange(1,3)], va, Quaternion(xq[4:7]...), ωa, xb, vb, qb, ωb, timestep),
					[xa; Dojo.vector(qa)]) * cat(I(3), Dojo.LVᵀmat(qa), dims=(1,2))
				@test norm(∇0 - ∇1, Inf) < 1.0e-8

				∇0 = Dojo.minimal_velocities_jacobian_configuration(:child, joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep)
				∇1 = force_to_jacobian_forward_diff(
					xq -> Dojo.minimal_velocities(joint, xa, va, qa, ωa, xq[Dojo.SUnitRange(1,3)], vb, Quaternion(xq[4:7]...), ωb, timestep),
					[xb; Dojo.vector(qb)]) * cat(I(3), Dojo.LVᵀmat(qb), dims=(1,2))
				@test norm(∇0 - ∇1, Inf) < 1.0e-8

				∇0 = Dojo.minimal_velocities_jacobian_velocity(:parent, joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep)
				∇1 = force_to_jacobian_forward_diff(
					vϕ -> Dojo.minimal_velocities(joint, xa, vϕ[Dojo.SUnitRange(1,3)], qa, vϕ[Dojo.SUnitRange(4,6)], xb, vb, qb, ωb, timestep),
					[va; ωa])
				@test norm(∇0 - ∇1, Inf) < 1.0e-8

				∇0 = Dojo.minimal_velocities_jacobian_velocity(:child, joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep)
				∇1 = force_to_jacobian_forward_diff(
					vϕ -> Dojo.minimal_velocities(joint, xa, va, qa, ωa, xb, vϕ[Dojo.SUnitRange(1,3)], qb, vϕ[Dojo.SUnitRange(4,6)], timestep),
					[vb; ωb])
				@test norm(∇0 - ∇1, Inf) < 1.0e-8
			end
		end
	end

	@testset "Minimal coordinate Jacobian" begin
		mech = DojoEnvironments.get_humanoid()
		for jointcon in mech.joints
			for joint in [jointcon.translational, jointcon.rotational]
				# generate random configuration in minimal space
				x = rand(Dojo.minimal_dimension(mech))

				# convert to maximal
				z = Dojo.minimal_to_maximal(mech, x)

				# extract body states
				Ne = Dojo.length(mech.joints)
				if Dojo.get_body(mech, jointcon.parent_id).name == :origin
					zp = [mech.origin.state.x2; mech.origin.state.v15; Dojo.vector(mech.origin.state.q2); mech.origin.state.ϕ15]
				else
					zp = z[(jointcon.parent_id - Ne - 1) * 13 .+ (1:13)]
				end
				zc = z[(jointcon.child_id - Ne - 1) * 13 .+ (1:13)]

				xa = SVector{3}(zp[1:3])
				# va = SVector{3}(zp[3 .+ (1:3)])
				qa = Quaternion(zp[6 .+ (1:4)]...)
				# ωa = SVector{3}(zp[10 .+ (1:3)])

				xb = SVector{3}(zc[1:3])
				# vb = SVector{3}(zc[3 .+ (1:3)])
				qb = Quaternion(zc[6 .+ (1:4)]...)
				# ωb = SVector{3}(zc[10 .+ (1:3)])

				Dojo.minimal_coordinates(joint, xa, qa, xb, qb)

				∇0 = Dojo.minimal_coordinates_jacobian_configuration(:parent, joint, xa, qa, xb, qb)
				∇1 = force_to_jacobian_forward_diff(
					xq -> Dojo.minimal_coordinates(joint, xq[1:3], Quaternion(xq[4:7]...), xb, qb),
					[xa; Dojo.vector(qa)]) * cat(I(3), Dojo.LVᵀmat(qa), dims=(1,2))
				@test norm(∇0 - ∇1, Inf) < 1.0e-8

				∇0 = Dojo.minimal_coordinates_jacobian_configuration(:child, joint, xa, qa, xb, qb)
				∇1 = force_to_jacobian_forward_diff(
					xq -> Dojo.minimal_coordinates(joint, xa, qa, xq[1:3], Quaternion(xq[4:7]...)),
					[xb; Dojo.vector(qb)]) * cat(I(3), Dojo.LVᵀmat(qb), dims=(1,2))
				@test norm(∇0 - ∇1, Inf) < 1.0e-8
			end
		end
	end

	@testset "Minimal to maximal Jacobian" begin
		function maximal_to_minimal_jacobian_fd(mechanism::Mechanism, z)
			J = force_to_jacobian_forward_diff(y -> Dojo.maximal_to_minimal(mechanism, y), z)
			G = attitude_jacobian(z, length(mechanism.bodies))
			return J * G
		end

		# TODO switch to ForwardDiff once it works
		function maximal_to_minimal_jacobian_fd_finite_diff(mechanism::Mechanism, z)
			J = force_to_jacobian_finite_diff(y -> Dojo.maximal_to_minimal(mechanism, y), z)
			G = attitude_jacobian(z, length(mechanism.bodies))
			return J * G
		end

		# TODO switch to ForwardDiff once it works
		function minimal_to_maximal_jacobian_fd(mechanism::Mechanism, x)
			J = force_to_jacobian_finite_diff(y -> Dojo.minimal_to_maximal(mechanism, y), x)
			z = minimal_to_maximal(mechanism, x)
			G = attitude_jacobian(z, length(mechanism.bodies))
			return G' * J
		end

		function ctrl!(mechanism, k)
			Dojo.set_input!(mechanism, 0.1 * srand(Dojo.input_dimension(mechanism)))
		end

		# n-pendulum
		mechanism = DojoEnvironments.get_mechanism(:npendulum,
			timestep=0.1,
			gravity=-9.81,
			num_bodies=1)

		base_angle = 0.3 * π
		Dojo.initialize!(mechanism, :npendulum,
			base_angle=base_angle)
		storage = Dojo.simulate!(mechanism, 1.0,
			record=true,
			verbose=false)

		x = Dojo.get_minimal_state(mechanism)
		z = Dojo.get_maximal_state(mechanism)
		u = zeros(Dojo.input_dimension(mechanism))

		@test norm(minimal_to_maximal(mechanism, x) - z) < 1.0e-5
		@test norm(Dojo.maximal_to_minimal(mechanism, z) - x) < 1.0e-8

		M_fd = maximal_to_minimal_jacobian_fd(mechanism, z)
		M_a = Dojo.maximal_to_minimal_jacobian(mechanism, z)
		@test size(M_fd) == size(M_a)
		@test norm(M_fd - M_a, Inf) < 1.0e-8

		N_fd = minimal_to_maximal_jacobian_fd(mechanism, Dojo.maximal_to_minimal(mechanism, z))
		N_a = Dojo.minimal_to_maximal_jacobian(mechanism, Dojo.maximal_to_minimal(mechanism, z))
		@test size(N_fd) == size(N_a)
		@test norm(N_fd - N_a, Inf) < 1.0e-5
		@test norm(diag(M_fd * N_fd) .- 1.0, Inf) < 1.0e-5
		@test norm(diag(M_a * N_a) .- 1.0, Inf) < 1.0e-8

		# n-pendulum
		mechanism = DojoEnvironments.get_mechanism(:npendulum,
			timestep=0.1,
			gravity=-9.81,
			num_bodies=2)
		base_angle = 0.3 * π
		Dojo.initialize!(mechanism, :npendulum,
			base_angle=base_angle)
		storage = Dojo.simulate!(mechanism, 1.0,
			record=true,
			verbose=false)

		x = Dojo.get_minimal_state(mechanism)
		z = Dojo.get_maximal_state(mechanism)
		u = zeros(Dojo.input_dimension(mechanism))

		@test norm(minimal_to_maximal(mechanism, x) - z) < 1.0e-5
		@test norm(Dojo.maximal_to_minimal(mechanism, z) - x) < 1.0e-8

		M_fd = maximal_to_minimal_jacobian_fd(mechanism, z)
		M_a = Dojo.maximal_to_minimal_jacobian(mechanism, z)
		@test size(M_fd) == size(M_a)
		@test norm(M_fd - M_a, Inf) < 1.0e-8

		N_fd = minimal_to_maximal_jacobian_fd(mechanism, Dojo.maximal_to_minimal(mechanism, z))
		N_a = Dojo.minimal_to_maximal_jacobian(mechanism, Dojo.maximal_to_minimal(mechanism, z))
		@test size(N_fd) == size(N_a)
		@test norm(N_fd - N_a, Inf) < 1.0e-5
		@test norm(diag(M_fd * N_fd) .- 1.0, Inf) < 1.0e-5
		@test norm(diag(M_a * N_a) .- 1.0, Inf) < 1.0e-5

		# sphere
		mechanism = DojoEnvironments.get_mechanism(:sphere,
			timestep=0.01,
			gravity=-9.81)
		Dojo.initialize!(mechanism, :sphere)
		storage = Dojo.simulate!(mechanism, 1.0,
			record=true,
			verbose=false)

		z = Dojo.get_maximal_state(mechanism)
		x = Dojo.get_minimal_state(mechanism)
		u = zeros(Dojo.input_dimension(mechanism))

		@test norm(minimal_to_maximal(mechanism, x) - z) < 1.0e-8
		@test norm(Dojo.maximal_to_minimal(mechanism, z) - x) < 1.0e-8

		M_fd = maximal_to_minimal_jacobian_fd_finite_diff(mechanism, z)
		M_a = Dojo.maximal_to_minimal_jacobian(mechanism, z)
		@test size(M_fd) == size(M_a)
		@test norm(M_fd - M_a, Inf) < 1.0e-5

		N_fd = minimal_to_maximal_jacobian_fd(mechanism, Dojo.maximal_to_minimal(mechanism, z))
		N_a = Dojo.minimal_to_maximal_jacobian(mechanism, Dojo.maximal_to_minimal(mechanism, z))
		@test size(N_fd) == size(N_a)
		@test norm(N_fd - N_a, Inf) < 1.0e-5
		@test norm(diag(M_fd * N_fd) .- 1.0, Inf) < 1.0e-5
		@test norm(diag(M_a * N_a) .- 1.0, Inf) < 1.0e-8

		# half cheetah
		mechanism = DojoEnvironments.get_mechanism(:halfcheetah,
			timestep=0.01,
			gravity=-9.81)
		Dojo.initialize!(mechanism, :halfcheetah)
		storage = Dojo.simulate!(mechanism, 1.0, ctrl!,
			record=true,
			verbose=false)

		z = Dojo.get_maximal_state(mechanism)
		x = Dojo.get_minimal_state(mechanism)
		u = zeros(Dojo.input_dimension(mechanism))

		@test norm(minimal_to_maximal(mechanism, x) - z) < 1.0e-5
		@test norm(Dojo.maximal_to_minimal(mechanism, z) - x) < 1.0e-8

		M_fd = maximal_to_minimal_jacobian_fd(mechanism, z)
		M_a = Dojo.maximal_to_minimal_jacobian(mechanism, z)
		@test size(M_fd) == size(M_a)
		@test norm(M_fd - M_a, Inf) < 1.0e-8

		N_fd = minimal_to_maximal_jacobian_fd(mechanism, Dojo.maximal_to_minimal(mechanism, z))
		N_a = Dojo.minimal_to_maximal_jacobian(mechanism, Dojo.maximal_to_minimal(mechanism, z))
		@test size(N_fd) == size(N_a)
		@test norm(N_fd - N_a, Inf) < 1.0e-5

		@test norm(diag(M_fd * N_fd) .- 1.0, Inf) < 1.0e-5
		@test norm(diag(M_a * N_a) .- 1.0, Inf) < 1.0e-5
		@test norm(diag(M_a * N_fd) .- 1.0, Inf) < 1.0e-5

		# atlas
		mechanism = DojoEnvironments.get_mechanism(:atlas,
			timestep=0.01,
			gravity=-9.81,
			friction_coefficient=0.5,
			damper=100.0,
			spring=1.0,
			contact_feet=true)
		Dojo.initialize!(mechanism, :atlas_stance,
			body_position=[0,0,0.5],
			body_orientation=[0.0,0.0,0.0])
		storage = Dojo.simulate!(mechanism, 1.0, ctrl!,
			record=true,
			verbose=false)

		z = Dojo.get_maximal_state(mechanism)
		x = Dojo.get_minimal_state(mechanism)
		u = zeros(Dojo.input_dimension(mechanism))

		@test norm(minimal_to_maximal(mechanism, x) - z) < 1.0e-5
		@test norm(Dojo.maximal_to_minimal(mechanism, z) - x) < 1.0e-8

		M_fd = maximal_to_minimal_jacobian_fd(mechanism, z)
		M_a = Dojo.maximal_to_minimal_jacobian(mechanism, z)
		@test size(M_fd) == size(M_a)
		@test norm(M_fd - M_a, Inf) < 1.0e-8

		N_fd = minimal_to_maximal_jacobian_fd(mechanism, Dojo.maximal_to_minimal(mechanism, z))
		N_a = Dojo.minimal_to_maximal_jacobian(mechanism, Dojo.maximal_to_minimal(mechanism, z))
		@test size(N_fd) == size(N_a)
		@test norm(N_fd - N_a, Inf) < 5.0e-5

		@test norm(diag(M_fd * N_fd) .- 1.0, Inf) < 1.0e-5
		@test norm(diag(M_a * N_a) .- 1.0, Inf) < 1.0e-5
	end

	@testset "Maximal to minimal Jacobian" begin

		function maximal_to_minimal_jacobian_fd_finite_diff(mechanism::Mechanism, z)
			J = force_to_jacobian_finite_diff(y -> maximal_to_minimal(mechanism, y), z)
			G = attitude_jacobian(z, length(mechanism.bodies))
			return J * G
		end

		function maximal_to_minimal_jacobian_fd(mechanism::Mechanism, z)
			J = force_to_jacobian_forward_diff(y -> maximal_to_minimal(mechanism, y), z)
			G = attitude_jacobian(z, length(mechanism.bodies))
			return J * G
		end

		# 5-link pendulum
		mech = DojoEnvironments.get_mechanism(:npendulum,
			timestep=0.01,
			gravity=-9.81,
			num_bodies=5)

		Random.seed!(100)
		base_angle = 0.3π
		Dojo.initialize!(mech, :npendulum,
			base_angle=base_angle)
		storage = Dojo.simulate!(mech, 1.0,
			record=true,
			verbose=false)

		Dojo.maximal_dimension(mech) == 13
		Dojo.minimal_dimension(mech) == 12
		z = Dojo.get_maximal_state(mech)

		attjac = Dojo.attitude_jacobian(z, length(mech.bodies))
		M_fd = maximal_to_minimal_jacobian_fd(mech, z)
		M_a = Dojo.maximal_to_minimal_jacobian(mech, z)
		@test size(M_fd) == size(M_a)
		@test norm(M_fd - M_a, Inf) < 1.0e-8

		# # sphere
		mech = DojoEnvironments.get_mechanism(:sphere,
			timestep=0.01,
			gravity=-9.81)
		Dojo.initialize!(mech, :sphere)
		storage = Dojo.simulate!(mech, 1.0,
			record=true,
			verbose=false)

		Dojo.maximal_dimension(mech)
		Dojo.minimal_dimension(mech)
		z = Dojo.get_maximal_state(mech)

		attjac = Dojo.attitude_jacobian(z, length(mech.bodies))
		M_fd = maximal_to_minimal_jacobian_fd_finite_diff(mech, z)
		M_a = Dojo.maximal_to_minimal_jacobian(mech, z)

		@test size(M_fd) == size(M_a)
		@test norm(M_fd - M_a, Inf) < 1.0e-6
	end
end



# ################################################################################
# # Test minimal velocities
# ################################################################################
# mech = DojoEnvironments.get_snake(gravity=0.00, num_bodies=2, damper=0.3, spring=0.2, joint_type=:Revolute)
# DojoEnvironments.initialize_snake!(mech)
# function ctrl!(m,k)
#     set_input!(m, 0.01*m.timestep*ones(minimal_dimension(m)))
# end
# storage = Dojo.simulate!(mech, 1.0, ctrl!)
# Dojo.visualize(mech, storage, vis=vis)
#
# mech.joints[2]
# rot0 = mech.joints[2].rotational
#
# timestep0 = 0.01
# xa0 = srand(3)
# qa0 = rand(QuatRotation).q
# va0 = srand(3)
# ϕa0 = srand(3)
# xb0 = srand(3)
# qb0 = rand(QuatRotation).q
# vb0 = srand(3)
# ϕb0 = srand(3)
#
# J0 = minimal_velocities_jacobian_configuration(:parent,
#     rot0, xa0, va0, qa0, ϕa0, xb0, vb0, qb0, ϕb0, timestep0)
# J1 = FiniteDiff.finite_difference_jacobian(
#     xq -> Dojo.minimal_velocities(rot0, xq[SUnitRange(1,3)], va0, Quaternion(xq[SUnitRange(4,7)]...,true), ϕa0,
#     xb0, vb0, qb0, ϕb0, timestep0),
#     [xa0; vector(qa0)]) * cat(I(3), LVᵀmat(qa0), dims=(1,2))
# norm(J0 - J1, Inf)
# norm(J0 - J1, Inf) < 1e-4
#
# J0 = minimal_velocities_jacobian_configuration(:child,
#     rot0, xa0, va0, qa0, ϕa0, xb0, vb0, qb0, ϕb0, timestep0)
# J1 = FiniteDiff.finite_difference_jacobian(
#     xq -> Dojo.minimal_velocities(rot0, xa0, va0, qa0, ϕa0,
#     xq[SUnitRange(1,3)], vb0, Quaternion(xq[SUnitRange(4,7)]...,true), ϕb0, timestep0),
#     [xb0; vector(qb0)]) * cat(I(3), LVᵀmat(qb0), dims=(1,2))
# norm(J0 - J1, Inf)
# norm(J0 - J1, Inf) < 1e-4
#
#
# J0 = minimal_velocities_jacobian_velocity(:parent,
#     rot0, xa0, va0, qa0, ϕa0, xb0, vb0, qb0, ϕb0, timestep0)
# J1 = FiniteDiff.finite_difference_jacobian(
#     vϕ -> Dojo.minimal_velocities(rot0, xa0, vϕ[SUnitRange(1,3)], qa0, vϕ[SUnitRange(4,6)], xb0, vb0, qb0, ϕb0, timestep0),
#     [va0; ϕa0])
# norm(J0 - J1, Inf)
# norm(J0 - J1, Inf) < 1e-4
#
# J0 = minimal_velocities_jacobian_velocity(:child,
#     rot0, xa0, va0, qa0, ϕa0, xb0, vb0, qb0, ϕb0, timestep0)
# J1 = FiniteDiff.finite_difference_jacobian(
#     vϕ -> Dojo.minimal_velocities(rot0, xa0, va0, qa0, ϕa0,
#         xb0, vϕ[SUnitRange(1,3)], qb0, vϕ[SUnitRange(4,6)], timestep0),
#     [vb0; ϕb0])
# norm(J0 - J1, Inf)
# norm(J0 - J1, Inf) < 1e-4
