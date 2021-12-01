
################################################################################
# Quadruped Trajectory Generator
################################################################################

# Trajectory template
function low_foot_trajectory(r::T, N::Int) where T
	# r = radius of the foot traj
	# low foot trajectory decomposed into N-1 steps (N waypoints)
	X = Vector(range(-r, r, length = N+1))
	low_traj = [[x, 0, 0.0] for x in X]
	return low_traj
end

function high_foot_trajectory(r::T, N::Int) where T
	# high foot trajectory decomposed into N+1 steps (N+2 waypoints)
	α = range(0, π, length = N+1)
	θ = π * (1 .- cos.(α)) / 2
	X = [r*cos(θi) for θi in θ]
	Z = [r*sin(θi)/2 for θi in θ] # the factor 2 flattens the foot trajectory
	high_traj = [[X[i], 0.0, Z[i]] for i = 1:N+1]
	return high_traj
end

function foot_trajectory(r::T, N::Int) where T
	low_traj = low_foot_trajectory(r, N)
	high_traj = high_foot_trajectory(r, N)
	traj = [low_traj[1:end-1]; high_traj[1:end-1]]
	return traj
end

# Inverse kinematics: given the position of the trunk, finds the knee and hip angles that will put the foot
# at the p_foot location
function IKquadruped(mechanism::Mechanism, p_trunk, p_foot; leg::Symbol = :FR)
	# starting point of the local search
	θ = [0.95, -1.5*0.95] # θhip, θknee
	for k = 1:10
		err = IKerror(mechanism, p_trunk, p_foot, θ; leg = leg)
		norm(err, Inf) < 1e-10 && continue
		∇ = FiniteDiff.finite_difference_jacobian(θ -> IKerror(mechanism, p_trunk, p_foot, θ; leg = leg), θ)
		θ -= ∇ \ err
	end
	return θ
end

function IKerror(mechanism::Mechanism, p_trunk, p_foot, θ; leg::Symbol = :FR)
	setPosition!(mechanism, geteqconstraint(mechanism, "auto_generated_floating_joint"), [p_trunk; zeros(3)])
	setPosition!(mechanism, geteqconstraint(mechanism, String(leg)*"_thigh_joint"), [θ[1]])
	setPosition!(mechanism, geteqconstraint(mechanism, String(leg)*"_calf_joint"), [θ[2]])

	foot = getbody(mechanism, String(leg)*"_calf")
	ineqcs = collect(mechanism.ineqconstraints)
	ineqc = ineqcs[findfirst(x -> x.parentid == foot.id, ineqcs)]
	p = contact_location(ineqc, foot)
	err = p - p_foot
	return err[[1,3]]
end

function quadruped_trajectory(mechanism::Mechanism{T}; Δt = 0.05, Δx = -0.04, r = 0.10, z = 0.29, N = 8, Ncycles = 1) where T
	pFR = [+ 0.13 + Δx, - 0.13205, 0.]
	pFL = [+ 0.13 + Δx, + 0.13205, 0.]
	pRR = [- 0.23 + Δx, - 0.13205, 0.]
	pRL = [- 0.23 + Δx, + 0.13205, 0.]

	t = reverse(foot_trajectory(r, N))
	t_dephased = [t[N+1:2N]; t[1:N]]
	# Foot positions
	tFR = [pFR + ti for ti in t]
	tFL = [pFL + ti for ti in t_dephased]
	tRR = [pRR + ti for ti in t_dephased]
	tRL = [pRL + ti for ti in t]

	# Leg angles
	p_trunk = [0,0,z]
	θFR = [IKquadruped(mechanism, p_trunk, tFR[i], leg = :FR) for i = 1:2N]
	θFL = [IKquadruped(mechanism, p_trunk, tFL[i], leg = :FL) for i = 1:2N]
	θRR = [IKquadruped(mechanism, p_trunk, tRR[i], leg = :RR) for i = 1:2N]
	θRL = [IKquadruped(mechanism, p_trunk, tRL[i], leg = :RL) for i = 1:2N]

	# adding angular velocities
	ωFR = [(θFR[i%10+1] - θFR[i]) / Δt for i = 1:2N]
	ωFL = [(θFL[i%10+1] - θFL[i]) / Δt for i = 1:2N]
	ωRR = [(θRR[i%10+1] - θRR[i]) / Δt for i = 1:2N]
	ωRL = [(θRL[i%10+1] - θRL[i]) / Δt for i = 1:2N]

	# Minimal Coordinates
	X = [Vector{T}([p_trunk; zeros(3); zeros(3); zeros(3);
		zeros(2); θFR[i][1]; ωFR[i][1]; θFR[i][2]; ωFR[i][2];
		zeros(2); θFL[i][1]; ωFL[i][1]; θFL[i][2]; ωFL[i][2];
		zeros(2); θRR[i][1]; ωRR[i][1]; θRR[i][2]; ωRR[i][2];
		zeros(2); θRL[i][1]; ωRL[i][1]; θRL[i][2]; ωRL[i][2];
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

# discretization
N = 9 # half period of the 1-step loop

# trunk trajectory
z = 0.29 # trunk height

# Half disk foot trajectory
r = 0.10 # foot traj radius

# System
mech = getmechanism(:quadruped, Δt = 0.05)
initialize!(mech, :quadruped, tran = [0,0,0.], v = [0,0,0.])

X = quadruped_trajectory(mech, r = 0.08, z = 0.27; Δt = 0.05, Δx = -0.03, N = 9, Ncycles = 10)
storage = generate_storage(mech, [min2max(mech, x) for x in X])
visualize(mech, storage, vis = vis)

collect(mech.ineqconstraints)
p_trunk = [0,0,0.31]
p_foot = [0.2,.0,0]
IKquadruped(mech, p_trunk, p_foot)

N = 6
low_traj = low_foot_trajectory(r, N)
high_traj = high_foot_trajectory(r, N)
traj = foot_trajectory(r, N)

plot()
scatter!([p[1] for p in low_traj], [p[3] for p in low_traj], )
scatter!([p[1] for p in high_traj], [p[3] for p in high_traj], )
scatter!([p[1] for p in traj], [p[3] for p in traj], )

N = 6
t = reverse(foot_trajectory(r, N))
t_dephased = [t[N+1:2N]; t[1:N]]
plt = plot(legend = false)
for i = 1:2N
	scatter!(plt, t[i][1:1], t[i][3:3], markersize = i+3)
	scatter!(plt, t_dephased[i][1:1], t_dephased[i][3:3] .+ 0.10, markersize = i+3)
	display(plt)
end
display(plt)

a = 10
a = 10
a = 10
a = 10
a = 10
a = 10

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
