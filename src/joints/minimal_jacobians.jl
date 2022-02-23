mech = get_slider()
initialize_slider!(mech, z1=0.1)

mech = get_pendulum()
initialize_pendulum!(mech, ϕ1=0.25 * π, ω1=0.1)

mech = get_halfcheetah()
initialize_halfcheetah!(mech)

simulate!(mech, 0.5)

mech = mechanism
joint = mech.joints[6]
bodya = get_body(mech, joint.parent_id) 
bodyb = get_body(mech, joint.child_id)
xa, va, qa, ϕa = current_configuration_velocity(bodya.state)
za = [xa; va; vector(qa); ϕa]
xb, qb = current_configuration(bodyb.state)
timestep = mech.timestep 

nu = control_dimension(joint)
nu_tra = control_dimension(joint.translational)
nu_rot = control_dimension(joint.rotational)

Arot = zerodimstaticadjoint(nullspace_mask(joint.rotational))
Atra = zerodimstaticadjoint(nullspace_mask(joint.translational))

x = minimal_coordinates_velocities(joint, bodya, bodyb, timestep)
Δx = x[SUnitRange(joint.minimal_index[1]...)]
Δθ = x[SUnitRange(joint.minimal_index[2]...)]
Δv = x[nu .+ SUnitRange(joint.minimal_index[1]...)]
Δϕ = x[nu .+ SUnitRange(joint.minimal_index[2]...)]
xmin = [Δx; Δθ; Δv; Δϕ]

# set_minimal_coordinates_velocities!(bodya, bodyb, joint, timestep)

# child wrt minimal

J_fd = FiniteDiff.finite_difference_jacobian(w -> set_minimal_coordinates_velocities!(bodya, bodyb, joint, timestep, xmin=w, zp=za), xmin)

J_a = minimal_coordinates_velocities_jacobian_minimal(joint,
    xa, va, qa, ϕa,
    timestep;
    Δx=Δx,
    Δθ=Δθ,
    Δv=Δv,
    Δϕ=Δϕ)

norm(Array(J_fd) - Array(J_a), Inf)

# position wrt maximal
J_fd = FiniteDiff.finite_difference_jacobian(w -> set_minimal_coordinates_velocities!(bodya, bodyb, joint, timestep, xmin=xmin, zp=w), za)

J_a = minimal_coordinates_velocities_jacobian_parent(joint,
    xa, va, qa, ϕa,
    xb, qb, timestep;
    Δx=Δx,
    Δθ=Δθ,
    Δv=Δv,
    Δϕ=Δϕ)

norm(Array(J_fd) - Array(J_a), Inf)

