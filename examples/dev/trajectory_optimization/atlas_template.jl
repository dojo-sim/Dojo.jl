
################################################################################
# Quadruped Trajectory Generator
################################################################################

# Trajectory template
function trunk_trajectory(z::T, r::T, N::Int) where T
	# z = vertical displacement
	# r = radius of the foot traj
	# trunk trajectory decomposed into N steps (N+1 waypoints)
	X = Vector(range(-r, r, length = N))
	trunk_traj = [[x, 0.0, z] for x in X]
	return trunk_traj
end

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
function inverse_kinematics(mechanism::Mechanism, p_base, p_foot; leg::Symbol = :r)
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

function altas_trajectory(mechanism::Mechanism{T}; Δt = 0.05, r = 0.10, z = 0.85, N = 12, Ncycles = 1) where T
	pL = [0, + 0.1145, 0]
	pR = [0, - 0.1145, 0]

	t = reverse(foot_trajectory(r, N))
	t_dephased = [t[N+1:2N]; t[1:N]]
	# Foot positions
	tL = [pL + ti for ti in t_dephased]
	tR = [pR + ti for ti in t]

	# Leg angles
	p_base = [0,0,z]
	θL = [inverse_kinematics(mechanism, p_base, tL[i], leg = :l) for i = 1:2N]
	θR = [inverse_kinematics(mechanism, p_base, tR[i], leg = :r) for i = 1:2N]

	# adding angular velocities
	ωL = [(θL[i%10+1] - θL[i]) / Δt for i = 1:2N]
	ωR = [(θR[i%10+1] - θR[i]) / Δt for i = 1:2N]

	# Minimal Coordinates
	X = [Vector{T}([
		p_base; zeros(3); zeros(3); zeros(3); # u_torso
		zeros(3); zeros(3); # pelvis
		[0, -θL[i][1], 0]; [0, -ωL[i][1], 0]; # l_uleg
		θL[i][2]; ωL[i][2]; # l_lleg
		[θL[i][1] - θL[i][2], 0]; [ωL[i][1] - ωL[i][2], 0]; # l_foot
		[0, -θR[i][1], 0]; [0, -ωR[i][1], 0]; # R_uleg
		θR[i][2]; ωR[i][2]; # R_lleg
		[θR[i][1] - θR[i][2], 0]; [ωR[i][1] - ωR[i][2], 0]; # R_foot
		]) for i = 1:2N]

	X = vcat([deepcopy(X) for i = 1:Ncycles]...)
	# adding forward displacement and velocity
	for i = 1:2N*Ncycles
		X[i][1] += (i-1) * 2r / N
		X[i][7] += 2r / N / Δt
	end
	return X
end


################################################################################
# Compute trajectory
################################################################################

mech = getmechanism(:atlas, Δt = 0.05, model_type = :armless, damper = 1000.0)
initialize!(mech, :atlas, tran = [1,0,0.0], rot = [0,0,0.], αhip = 0.0, αknee = 0.0)

@elapsed storage = simulate!(mech, 0.55, record = true, undercut = Inf,
    solver = :mehrotra!, verbose = true)

visualize(mech, storage, vis = vis)

atlas_trajectory(mech; Δt = 0.05, r = 0.10, z = 0.85, N = 12, Ncycles = 1)




x0 = [[0,0,0.9385]; zeros(9);
	zeros(6);
	zeros(6);
	zeros(6);
	zeros(18)]
z0 = min2max(mech, x0)
visualizeMaxCoord(mech, z0, vis)

center_of_mass(mech, storage, 1)
contact_location(mech)[1:4]

bodies = collect(mech.bodies)
eqcs = collect(mech.eqconstraints)
ineqcs = collect(mech.ineqconstraints)
getfield.(ineqcs, :parentid)
getbody(mech, 15)
nx = minCoordDim(mech)


p_base = [0, 0, 0.9385]
p_base = [0, 0, 0.88]
p_foot = [0, 0, 0.00]
θ = [0.4, 0.8]
IKerror(mech, p_base, p_foot, θ; leg = :r)
inverse_kinematics(mech, p_base, p_foot; leg = :r)


# discretization
N = 8 # half period of the 1-step loop

# trunk trajectory
z = 0.29 # trunk height

# Half disk foot trajectory
r = 0.10 # foot traj radius

# System
mech = getmechanism(:quadruped, Δt = 0.05)
initialize!(mech, :quadruped, tran = [0,0,0.], v = [0,0,0.])

X = quadruped_trajectory(mech, 0.0, r = 0.10, z = 0.29; Δt = 0.05, N = N, Ncycles = 10)
storage = generate_storage(mech, [min2max(mech, x) for x in X])
visualize(mech, storage, vis = vis)

collect(mech.ineqconstraints)
p_trunk = [0,0,0.31]
p_foot = [0.2,0,0]
inverse_kinematics(mech, p_trunk, p_foot)

N = 6
low_traj = low_foot_trajectory(r, 0.0, y, Nfoot)
high_traj = high_foot_trajectory(r, 0.0, y, Nfoot)
traj = foot_trajectory(r, 0.0, y, Nfoot)
trunk_traj = trunk_trajectory(0.0, z, r, Nfoot)

plot()
scatter!([p[1] for p in low_traj], [p[3] for p in low_traj], )
scatter!([p[1] for p in high_traj], [p[3] for p in high_traj], )
scatter!([p[1] for p in traj], [p[3] for p in traj], )

plot()
scatter!([p[1] for p in tFR], [p[3] for p in tFR])
scatter!([p[1] for p in tFL], [p[3] for p in tFL])
scatter!([p[1] for p in tRR], [p[3] for p in tRR])
scatter!([p[1] for p in tRR], [p[3] for p in tRR])


# contact_location(mech, ineqcs[1])
# contact_location(mech, ineqcs[2])
# contact_location(mech, ineqcs[3])
# contact_location(mech, ineqcs[4])




#
#
#
# ## state space
# Nb = length(mech.bodies)
# n = minCoordDim(mech)
# m = 12
#
# function potato_dynamics(x, u, Δt, m, g)
# 	# x = [x,y,z,ẋ,ẏ,ż]
# 	gv = [0, 0, g]
# 	ẋ = [x[4:6]; u ./ (m*Δt) + gv]
# 	x̄ = x + ẋ * Δt
# 	return x̄
# end
#
# trunk = getbody(mech, "trunk")
# x2_trunk = trunk.state.x2[1]
# v15_trunk = trunk.state.v15
#
#
# U_potato = [fill([0,0,8], 6); fill(zeros(3), 15)]
# x_potato = [x2_trunk; v15_trunk]
# X_potato = []
# for t = 1:21
# 	mass = sum(getfield.(mech.bodies, :m))
# 	alt = x_potato[3] - 0.48
# 	if t > 13
# 		u_potato = -[0, 0, mech.Δt * mass * mech.g + 20*alt + 10*x_potato[6]]
# 	else
# 		u_potato = U_potato[t]
# 	end
# 	push!(X_potato, x_potato)
# 	x_potato = potato_dynamics(x_potato, u_potato, mech.Δt, mass, mech.g)
# end
# plot()
# plot!([x[1] for x in X_potato], linewidth = 5.0)
# plot!([x[2] for x in X_potato], linewidth = 5.0)
# plot!([x[3] for x in X_potato], linewidth = 5.0)
#
