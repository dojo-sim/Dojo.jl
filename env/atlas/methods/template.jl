
################################################################################
# Quadruped Trajectory Generator
################################################################################

# Trajectory template
function low_foot_trajectory(r::T, N::Int) where T
	# r = radius of the foot traj
	# low foot trajectory decomposed into N steps (N+1 waypoints)
	X = Vector(range(-r, r, length = N+1))
	low_traj = [[x, 0, 0.0] for x in X]
	return low_traj
end

function high_foot_trajectory(r::T, N::Int) where T
	# high foot trajectory decomposed into N steps (N+1 waypoints)
	α = range(0, π, length = N+1)
	θ = π * (1 .- cos.(α)) / 2
	X = [r*cos(θi) for θi in θ]
	Z = [r*sin(θi)/2 for θi in θ] # the facto 2 flatten the foo trajectory
	high_traj = [[X[i], 0.0, Z[i]] for i = 1:N+1]
	return high_traj
end

function foot_trajectory(r::T, N::Int) where T
	low_traj = low_foot_trajectory(r, N)
	high_traj = high_foot_trajectory(r, N)
	traj = [low_traj[1:end-1]; high_traj[1:end-1]]
	return traj
end

# Inverse kinematics: given the position of the first link, finds the knee and hip angles that will put the foot
# at the p_foot location
function IKatlas(mechanism::Mechanism, p_base, p_foot; leg::Symbol = :r)
	# starting point of the local search
	θ = [0.5, 1.0] # θhip, θknee
	for k = 1:10
		err = IKerror(mechanism, p_base, p_foot, θ; leg = leg)
		norm(err, Inf) < 1e-10 && continue
		∇ = FiniteDiff.finite_difference_jacobian(θ -> IKerror(mechanism, p_base, p_foot, θ; leg = leg), θ)
		θ -= ∇ \ err
	end
	return θ
end

function IKerror(mechanism::Mechanism, p_base, p_foot, θ; leg::Symbol = :r)
	setPosition!(mechanism, geteqconstraint(mechanism, "auto_generated_floating_joint"), [p_base; zeros(3)])
	setPosition!(mechanism, geteqconstraint(mechanism, String(leg)*"_leg_hpxyz"), [0.0, -θ[1], 0.0])
	setPosition!(mechanism, geteqconstraint(mechanism, String(leg)*"_leg_kny"), [θ[2]])
	setPosition!(mechanism, geteqconstraint(mechanism, String(leg)*"_leg_akxy"), [θ[1]-θ[2], 0.0])

	foot = getbody(mechanism, String(leg)*"_foot")
	ineqcs = collect(mechanism.ineqconstraints)
	foot_ineqcs = ineqcs[findall(x -> x.parentid == foot.id, ineqcs)]
	p = mean([contact_location(ineqc, foot) for ineqc in foot_ineqcs]) # average of all contact locations for one foot
	err = p - p_foot
	return err[[1,3]]
end

function atlas_trajectory(mechanism::Mechanism{T}; Δt = 0.05, r = 0.10, z = 0.85, N = 12, Ncycles = 1) where T
	pL = [0, + 0.1145, 0]
	pR = [0, - 0.1145, 0]

	t = reverse(foot_trajectory(r, N))
	t_dephased = [t[N+1:2N]; t[1:N]]
	# Foot positions
	tL = [pL + ti for ti in t_dephased]
	tR = [pR + ti for ti in t]

	# Leg angles
	p_base = [0,0,z]
	θL = [IKatlas(mechanism, p_base, tL[i], leg = :l) for i = 1:2N]
	θR = [IKatlas(mechanism, p_base, tR[i], leg = :r) for i = 1:2N]

	# adding angular velocities
	ωL = [(θL[i%(2N)+1] - θL[i]) / Δt for i = 1:2N]
	ωR = [(θR[i%(2N)+1] - θR[i]) / Δt for i = 1:2N]

	# Minimal Coordinates
	nx = minCoordDim(mechanism)
	X = [zeros(nx) for i = 1:2N]
	for i = 1:2N
		setPosition!(mechanism, geteqconstraint(mechanism, "auto_generated_floating_joint"), [p_base; zeros(3)])
		setPosition!(mechanism, geteqconstraint(mechanism, "l_leg_hpxyz"), [0.0, -θL[i][1], 0.0])
		setPosition!(mechanism, geteqconstraint(mechanism, "l_leg_kny"), [θL[i][2]])
		setPosition!(mechanism, geteqconstraint(mechanism, "l_leg_akxy"), [θL[i][1]-θL[i][2], 0.0])

		setPosition!(mechanism, geteqconstraint(mechanism, "r_leg_hpxyz"), [0.0, -θR[i][1], 0.0])
		setPosition!(mechanism, geteqconstraint(mechanism, "r_leg_kny"), [θR[i][2]])
		setPosition!(mechanism, geteqconstraint(mechanism, "r_leg_akxy"), [θR[i][1]-θR[i][2], 0.0])

		setVelocity!(mechanism, geteqconstraint(mechanism, "auto_generated_floating_joint"), [zeros(3); zeros(3)])
		setVelocity!(mechanism, geteqconstraint(mechanism, "l_leg_hpxyz"), [0.0, -ωL[i][1], 0.0])
		setVelocity!(mechanism, geteqconstraint(mechanism, "l_leg_kny"), [ωL[i][2]])
		setVelocity!(mechanism, geteqconstraint(mechanism, "l_leg_akxy"), [ωL[i][1]-ωL[i][2], 0.0])

		setVelocity!(mechanism, geteqconstraint(mechanism, "r_leg_hpxyz"), [0.0, -ωR[i][1], 0.0])
		setVelocity!(mechanism, geteqconstraint(mechanism, "r_leg_kny"), [ωR[i][2]])
		setVelocity!(mechanism, geteqconstraint(mechanism, "r_leg_akxy"), [ωR[i][1]-ωR[i][2], 0.0])
		X[i] .= getMinState(mechanism)
	end

	X = vcat([deepcopy(X) for i = 1:Ncycles]...)
	# adding forward displacement and velocity
	for i = 1:2N*Ncycles
		X[i][1] += (i-1) * 2r / N
		X[i][7] += 2r / N / Δt
	end
	return X
end


# ################################################################################
# # Compute trajectory
# ################################################################################
#
# mech = getmechanism(:atlas, Δt = 0.05, model_type = :armless, damper = 1000.0)
# initialize!(mech, :atlas, tran = [1,0,0.0], rot = [0,0,0.], αhip = 0.0, αknee = 0.0)
#
# @elapsed storage = simulate!(mech, 0.55, record = true, undercut = Inf,
#     solver = :mehrotra!, verbose = true)
# visualize(mech, storage, vis = vis)
#
# X = atlas_trajectory(mech; Δt = 0.05, r = 0.10, z = 0.89, N = 12, Ncycles = 5)
# zref = [min2max(mech, x) for x in X]
# storage = generate_storage(mech, zref)
# visualize(mech, storage, vis = vis)
#
#
# bodies = collect(mech.bodies)
# eqcs = collect(mech.eqconstraints)
# ineqcs = collect(mech.ineqconstraints)
# getfield.(ineqcs, :parentid)
# getbody(mech, 15)
# nx = minCoordDim(mech)
#
#
# p_base = [0, 0, 0.9385]
# p_base = [0, 0, 0.88]
# p_foot = [0, 0, 0.00]
# θ = [0.4, 0.8]
# IKerror(mech, p_base, p_foot, θ; leg = :r)
# IKatlas(mech, p_base, p_foot; leg = :r)
