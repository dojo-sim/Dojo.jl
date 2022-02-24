mech = get_slider()
initialize_slider!(mech, z1=0.1)

mech = get_pendulum()
initialize_pendulum!(mech, ϕ1=0.25 * π, ω1=0.1)

mech = get_halfcheetah()
initialize_halfcheetah!(mech)

mech = get_humanoid()
initialize_humanoid!(mech)

simulate!(mech, 0.5)

# mech = mechanism
joint = mech.joints[2]
bodya = get_body(mech, joint.parent_id) 
bodyb = get_body(mech, joint.child_id)
xa, va, qa, ϕa = current_configuration_velocity(bodya.state)
za = [xa; va; vector(qa); ϕa]
xb, vb, qb, ϕb = current_configuration_velocity(bodyb.state)
timestep = mech.timestep 

# nu = control_dimension(joint)
# nu_tra = control_dimension(joint.translational)
# nu_rot = control_dimension(joint.rotational)

# Arot = zerodimstaticadjoint(nullspace_mask(joint.rotational))
# Atra = zerodimstaticadjoint(nullspace_mask(joint.translational))

# x = minimal_coordinates_velocities(joint, bodya, bodyb, timestep)
# Δx = x[SUnitRange(joint.minimal_index[1]...)]
# Δθ = x[SUnitRange(joint.minimal_index[2]...)]
# Δv = x[nu .+ SUnitRange(joint.minimal_index[1]...)]
# Δϕ = x[nu .+ SUnitRange(joint.minimal_index[2]...)]
# xmin = [Δx; Δθ; Δv; Δϕ]

minimal_velocities(joint.rotational,
		xa, va, qa, ϕa,
		xb, vb, qb, ϕb,
		timestep)

J_fd = FiniteDiff.finite_difference_jacobian(w -> minimal_velocities(joint.rotational,
    SVector{3}(w[1:3]), va, UnitQuaternion(w[3 .+ (1:4)]..., false), ϕa,
    xb, vb, qb, ϕb,
    timestep),
    [xa; vector(qa)]) * cat(1.0 * I(3), LVᵀmat(qa), dims=(1,2))

J_a = _minimal_velocities_jacobian_configuration(:parent, joint.rotational,
    xa, va, qa, ϕa,
    xb, vb, qb, ϕb,
    timestep)

norm(J_fd - J_a, Inf)

J_fd = FiniteDiff.finite_difference_jacobian(w -> minimal_velocities(joint.rotational,
    xa, va, qa, ϕa,
    SVector{3}(w[1:3]), vb, UnitQuaternion(w[3 .+ (1:4)]..., false), ϕb,
    timestep),
    [xb; vector(qb)]) * cat(1.0 * I(3), LVᵀmat(qb), dims=(1,2))

J_a = _minimal_velocities_jacobian_configuration(:child, joint.rotational,
    xa, va, qa, ϕa,
    xb, vb, qb, ϕb,
    timestep)

norm(J_fd - J_a, Inf)

J_fd = FiniteDiff.finite_difference_jacobian(w -> minimal_velocities(joint.rotational,
    xa, SVector{3}(w[1:3]), qa, SVector{3}(w[3 .+ (1:3)]),
    xb, vb, qb, ϕb,
    timestep),
    [va; ϕa]) 

J_a = _minimal_velocities_jacobian_velocity(:parent, joint.rotational,
    xa, va, qa, ϕa,
    xb, vb, qb, ϕb,
    timestep)

norm(J_fd - J_a, Inf)


minimal_velocities(joint.translational,
		xa, va, qa, ϕa,
		xb, vb, qb, ϕb,
		timestep)

J_fd = FiniteDiff.finite_difference_jacobian(w -> minimal_velocities(joint.translational,
    SVector{3}(w[1:3]), va, UnitQuaternion(w[3 .+ (1:4)]..., false), ϕa,
    xb, vb, qb, ϕb,
    timestep),
    [xa; vector(qa)]) * cat(1.0 * I(3), LVᵀmat(qa), dims=(1,2))

J_a = _minimal_velocities_jacobian_configuration(:parent, joint.translational,
    xa, va, qa, ϕa,
    xb, vb, qb, ϕb,
    timestep)

norm(J_fd - J_a, Inf)

J_fd = FiniteDiff.finite_difference_jacobian(w -> minimal_velocities(joint.translational,
    xa, va, qa, ϕa,
    SVector{3}(w[1:3]), vb, UnitQuaternion(w[3 .+ (1:4)]..., false), ϕb,
    timestep),
    [xb; vector(qb)]) * cat(1.0 * I(3), LVᵀmat(qb), dims=(1,2))

J_a = _minimal_velocities_jacobian_configuration(:child, joint.translational,
    xa, va, qa, ϕa,
    xb, vb, qb, ϕb,
    timestep)

norm(J_fd - J_a, Inf)

J_fd = FiniteDiff.finite_difference_jacobian(w -> minimal_velocities(joint.translational,
    xa, SVector{3}(w[1:3]), qa, SVector{3}(w[3 .+ (1:3)]),
    xb, vb, qb, ϕb,
    timestep),
    [va; ϕa])

J_a = _minimal_velocities_jacobian_velocity(:parent, joint.translational,
    xa, va, qa, ϕa,
    xb, vb, qb, ϕb,
    timestep)

norm(J_fd - J_a, Inf)

J_fd = FiniteDiff.finite_difference_jacobian(w -> minimal_velocities(joint.translational,
    xa, va, qa, ϕa,
    xb, SVector{3}(w[1:3]), qb, SVector{3}(w[3 .+ (1:3)]),
    timestep),
    [vb; ϕb])

J_a = _minimal_velocities_jacobian_velocity(:child, joint.translational,
    xa, va, qa, ϕa,
    xb, vb, qb, ϕb,
    timestep)

norm(J_fd - J_a, Inf)

# x = rand(minimal_dimension(mech))
# z = minimal_to_maximal(mech, x)

# Ne = length(mech.joints)
# zp = z[(joint.parent_id - Ne - 1) * 13 .+ (1:13)]
# zc = z[(joint.child_id - Ne - 1) * 13 .+ (1:13)]

# xp = zp[1:3] 
# vp = zp[3 .+ (1:3)]
# qp = UnitQuaternion(zp[6 .+ (1:4)]..., false)
# ϕp = zp[10 .+ (1:3)]

# xc = zc[1:3] 
# vc = zc[3 .+ (1:3)]
# qc = UnitQuaternion(zc[6 .+ (1:4)]..., false)
# ϕc = zc[10 .+ (1:3)]


a = 1
