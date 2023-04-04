"""
    Quadruped <: Environment

    four-legged dog-like robot 
    https://www.unitree.com/products/a1/
"""
struct Quadruped end

function quadruped(; 
    representation=:minimal, 
    timestep=0.05, 
    gravity=[0.0; 0.0; -9.81], 
    friction_coefficient=0.8,
    damper=0.0, 
    spring=0.0, 
    parse_damper=true,
    info=nothing,
    seed=1, 
    contact_feet=true,
    contact_body=true,
    vis=Visualizer(), 
    name=:robot,
    infeasible_control=false,
    opts_step=SolverOptions(rtol=3.0e-4, btol=3.0e-4, undercut=1.5),
    opts_grad=SolverOptions(rtol=3.0e-4, btol=3.0e-4, undercut=1.5),
    T=Float64)

    mechanism = get_mechanism(:quadruped;
        timestep, 
        contact_feet,
        contact_body, 
        gravity, 
        friction_coefficient, 
        parse_damper,
        damper, 
        spring)

    initialize!(mechanism, :quadruped)

    if representation == :minimal
        nx = minimal_dimension(mechanism)
    elseif representation == :maximal
        nx = maximal_dimension(mechanism)
    end

    nu_inf = input_dimension(mechanism)
    nu = infeasible_control ? nu_inf : nu_inf - 6 # remove first 6 controls
    no = nx

    aspace = BoxSpace(nu, 
        low=(-1.0e8 * ones(nu)), 
        high=(1.0e8 * ones(nu)))
    ospace = BoxSpace(no, 
        low=(-Inf * ones(no)), 
        high=(Inf * ones(no)))

    rng = MersenneTwister(seed)
    z = get_maximal_state(mechanism)
    x = representation == :minimal ? maximal_to_minimal(mechanism, z) : z
    fx = zeros(nx, nx)
    fu = zeros(nx, nu)

    u_prev = zeros(nu)
    control_mask = infeasible_control ? I(nu) : [zeros(nu, 6) I(nu)]
    control_scaling = Diagonal(ones(nu))

    build_robot(mechanism, 
        vis=vis, name=name)

    TYPES = [Quadruped, T, typeof(mechanism), typeof(aspace), typeof(ospace), typeof(info)]
    Environment{TYPES...}(mechanism, representation, aspace, ospace,
        x, fx, fu,
        u_prev, control_mask' * control_scaling,
        nx, nu, no,
        info,
        [rng], vis,
        opts_step, opts_grad)
end






################################################################################
# Quadruped Trajectory Generator
################################################################################

# Trajectory template
function low_foot_trajectory(r::T, width_scale::T, N::Int) where T
	# r = radius of the foot traj
	# low foot trajectory decomposed into N-1 steps (N waypoints)
	X = Vector(range(-width_scale*r, width_scale*r, length = N+1))
	low_traj = [[x, 0, 0.0] for x in X]
	return low_traj
end

function high_foot_trajectory(r::T, width_scale::T, height_scale::T, N::Int) where T
	# high foot trajectory decomposed into N+1 steps (N+2 waypoints)
	α = range(0, π, length = N+1)
	θ = π * (1 .- cos.(α)) / 2
	X = [width_scale*r*cos(θi) for θi in θ]
	Z = [height_scale*r*sin(θi) for θi in θ] # the factor α flattens the foot trajectory
	high_traj = [[X[i], 0.0, Z[i]] for i = 1:N+1]
	return high_traj
end

function foot_trajectory(r::T, width_scale::T, height_scale::T, N::Int) where T
	low_traj = low_foot_trajectory(r, width_scale, N)
	high_traj = high_foot_trajectory(r, width_scale, height_scale, N)
	traj = [low_traj[1:end-1]; high_traj[1:end-1]]
	return traj
end

# Inverse kinematics: given the position of the trunk, finds the knee and hip angles that will put the foot
# at the p_foot location
function IKquadruped(mechanism::Mechanism, p_trunk, p_foot;
	leg::Symbol=:FR)

	# starting point of the local search
	θ = [0.95, -1.5*0.95] # θhip, θknee
	for k = 1:10
		err = QuadrupedIKerror(mechanism, p_trunk, p_foot, θ; leg = leg)
		norm(err, Inf) < 1e-10 && continue
		∇ = FiniteDiff.finite_difference_jacobian(θ -> QuadrupedIKerror(mechanism, p_trunk, p_foot, θ; leg = leg), θ)
		θ -= ∇ \ err
	end
	return θ
end

function QuadrupedIKerror(mechanism::Mechanism, p_trunk, p_foot, θ;
	leg::Symbol=:FR)

	Dojo.set_minimal_coordinates!(mechanism, Dojo.get_joint(mechanism, :floating_base), [p_trunk; zeros(3)])
	Dojo.set_minimal_coordinates!(mechanism, Dojo.get_joint(mechanism, Symbol(String(leg)*"_thigh_joint")), [θ[1]])
	Dojo.set_minimal_coordinates!(mechanism, Dojo.get_joint(mechanism, Symbol(String(leg)*"_calf_joint")), [θ[2]])

	foot = Dojo.get_body(mechanism, Symbol(String(leg)*"_calf"))
	contacts = collect(mechanism.contacts)
	contact = contacts[findfirst(x -> x.parent_id == foot.id, contacts)]
	p = Dojo.contact_location(contact, foot)
	err = p - p_foot
	return err[[1,3]]
end

function quadruped_trajectory(mechanism::Mechanism{T};
	timestep=0.05,
	Δx = -0.04,
	Δfront = 0.05,
	width_scale = 1.0,
	height_scale = 0.5,
	r = 0.10,
	z = 0.29,
	N = 8,
	Ncycles = 1) where T

	# pFR = [+ 0.13 + Δx, - 0.13205, 0.]
	# pFL = [+ 0.13 + Δx, + 0.13205, 0.]
	# pRR = [- 0.23 + Δx, - 0.13205, 0.]
	# pRL = [- 0.23 + Δx, + 0.13205, 0.]
	pRL = [- 0.23 + Δx, + 0.13205, 0.]
	pRR = [- 0.23 + Δx, - 0.13205, 0.]
	pFL = [+ 0.13 + Δx + Δfront, + 0.13205, 0.]
	pFR = [+ 0.13 + Δx + Δfront, - 0.13205, 0.]

	t = reverse(foot_trajectory(r, width_scale, height_scale, N))
	t_dephased = [t[N+1:2N]; t[1:N]]
	# Foot positions
	tRL = [pRL + ti for ti in t]
	tRR = [pRR + ti for ti in t_dephased]
	tFL = [pFL + ti for ti in t_dephased]
	tFR = [pFR + ti for ti in t]

	# Leg angles
	p_trunk = [0,0,z]
	θRL = [IKquadruped(mechanism, p_trunk, tRL[i], leg = :RL) for i = 1:2N]
	θRR = [IKquadruped(mechanism, p_trunk, tRR[i], leg = :RR) for i = 1:2N]
	θFL = [IKquadruped(mechanism, p_trunk, tFL[i], leg = :FL) for i = 1:2N]
	θFR = [IKquadruped(mechanism, p_trunk, tFR[i], leg = :FR) for i = 1:2N]

	# adding angular velocities
	ωRL = [(θRL[i%(2N)+1] - θRL[i]) / timestep for i = 1:2N]
	ωRR = [(θRR[i%(2N)+1] - θRR[i]) / timestep for i = 1:2N]
	ωFL = [(θFL[i%(2N)+1] - θFL[i]) / timestep for i = 1:2N]
	ωFR = [(θFR[i%(2N)+1] - θFR[i]) / timestep for i = 1:2N]

	# Minimal Coordinates
	X = [Vector{T}([p_trunk; zeros(3); zeros(3); zeros(3);
		zeros(2); θRL[i][1]; ωRL[i][1]; θRL[i][2]; ωRL[i][2];
		zeros(2); θRR[i][1]; ωRR[i][1]; θRR[i][2]; ωRR[i][2];
		zeros(2); θFL[i][1]; ωFL[i][1]; θFL[i][2]; ωFL[i][2];
		zeros(2); θFR[i][1]; ωFR[i][1]; θFR[i][2]; ωFR[i][2];
		]) for i = 1:2N]

	X = vcat([deepcopy(X) for i = 1:Ncycles]...)
	# adding forward displacement and velocity
	for i = 1:2N*Ncycles
		X[i][1] += (i-1) * 2 * width_scale * r / N
		X[i][7] += 2 * width_scale * r / N / timestep
	end
	return X
end

function gravity_compensation(mechanism::Mechanism)
    # only works with revolute joints for now
    nu = Dojo.input_dimension(mechanism)
    u = zeros(nu)
    off  = 0
    for joint in mechanism.joints
        nu = Dojo.input_dimension(joint)
        if joint.parent_id != 0
            body = Dojo.get_body(mechanism, joint.parent_id)
            rot = joint.rotational
            A = Matrix(nullspace_mask(rot))
            input = Dojo.spring_impulses(mechanism, joint, body)
            F = input[1:3]
            τ = input[4:6]
            u[off .+ (1:nu)] = -A * τ
        else
            @warn "need to treat the joint to origin"
        end
        off += nu
    end
    return u
end


# ################################################################################
# # Compute trajectory
# ################################################################################
#
# # discretization
# N = 9 # half period of the 1-step loop
#
# # trunk trajectory
# z = 0.29 # trunk height
#
# # Half disk foot trajectory
# r = 0.10 # foot traj radius
#
# # System
# mech = get_mechanism(:quadruped; timestep=0.05)
# initialize!(mech, :quadruped, tran = [0,0,0.], v = [0,0,0.])
#
# X = quadruped_trajectory(mech, r = 0.08, z = 0.27; timestep=0.05, Δx = -0.03, N = 9, Ncycles = 10)
# storage = generate_storage(mech, [minimal_to_maximal(mech, x) for x in X])
# visualize(mech, storage, vis=vis)
#
# collect(mech.contacts)
# p_trunk = [0,0,0.31]
# p_foot = [0.2,.0,0]
# IKquadruped(mech, p_trunk, p_foot)
#
# N = 6
# low_traj = low_foot_trajectory(r, N)
# high_traj = high_foot_trajectory(r, N)
# traj = foot_trajectory(r, N)
#
# plot()
# scatter!([p[1] for p in low_traj], [p[3] for p in low_traj], )
# scatter!([p[1] for p in high_traj], [p[3] for p in high_traj], )
# scatter!([p[1] for p in traj], [p[3] for p in traj], )
#
# N = 6
# t = reverse(foot_trajectory(r, N))
# t_dephased = [t[N+1:2N]; t[1:N]]
# plt = plot(legend = false)
# for i = 1:2N
# 	scatter!(plt, t[i][1:1], t[i][3:3], markersize = i+3)
# 	scatter!(plt, t_dephased[i][1:1], t_dephased[i][3:3] .+ 0.10, markersize = i+3)
# 	display(plt)
# end
# display(plt)
#
# a = 10
# a = 10
# a = 10
# a = 10
# a = 10
# a = 10
#
# plot()
# scatter!([p[1] for p in tFR], [p[3] for p in tFR])
# scatter!([p[1] for p in tFL], [p[3] for p in tFL])
# scatter!([p[1] for p in tRR], [p[3] for p in tRR])
# scatter!([p[1] for p in tRR], [p[3] for p in tRR])


# contact_location(mech, contacts[1])
# contact_location(mech, contacts[2])
# contact_location(mech, contacts[3])
# contact_location(mech, contacts[4])




#
#
#
# ## dimensions
# Nb = length(mech.bodies)
# n = minimal_dimension(mech)
# m = 12
#
# function potato_dynamics(x, u, timestep, m, g)
# 	# x = [x,y,z,ẋ,ẏ,ż]
# 	gv = [0, 0, g]
# 	ẋ = [x[4:6]; u ./ (m*timestep) + gv]
# 	x̄ = x + ẋ * timestep
# 	return x̄
# end
#
# trunk = get_body(mech, "trunk")
# x2_trunk = trunk.state.x2
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
# 		u_potato = -[0, 0, mass * mechanism.gravity + 20*alt + 10*x_potato[6]]
# 	else
# 		u_potato = U_potato[t]
# 	end
# 	push!(X_potato, x_potato)
# 	x_potato = potato_dynamics(x_potato, u_potato, mech.timestep, mass, mechanism.gravity)
# end
# plot()
# plot!([x[1] for x in X_potato], linewidth = 5.0)
# plot!([x[2] for x in X_potato], linewidth = 5.0)
# plot!([x[3] for x in X_potato], linewidth = 5.0)
#
