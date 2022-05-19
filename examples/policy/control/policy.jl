




maximum(abs.(vcat(u_sol...)))

dynamics_jacobian_state(dx, env, x, u, w)
dynamics_jacobian_input(du, env, x, u, w)
A = [zeros(n,n) for i = 1:T-1]
B = [zeros(n,m) for i = 1:T-1]
for i = 1:T-1
    dynamics_jacobian_state(A[i], env, x_sol[i], u_sol[i], zeros(0))
    dynamics_jacobian_input(B[i], env, x_sol[i], u_sol[i], zeros(0))
end
qt = [0.3; 0.05; 0.05;
    5e-1 * ones(3);
    1e-6 * ones(3);
    1e-6 * ones(3);
    fill([2, 1e-6], 12)...]

rt = timestep * 100.0 * ones(m)

Q = [Diagonal(qt) for i = 1:T]
R = [Diagonal(rt) for i = 1:T-1]

K_tv, P_tv = tvlqr(A,B,Q,R)

################################################################################
# Simulate accurately the TVLQR policy
################################################################################
# divide timestep
S = 5
mech_sim = get_mechanism(:quadruped,
    contact_body=true,
    timestep=timestep/S,
    gravity=gravity,
    friction_coefficient=friction_coefficient,
    damper=damper,
    spring=spring)

# controller
function ctrl!(m, k; x_sol=x_sol, u_sol=u_sol, K_sol=K_sol)
    N = length(u_sol)
    nu = input_dimension(m)
    @show Int(floor((k-1)/S))
    ind = (Int(floor((k-1)/S))) % N + 1

    x = get_minimal_state(m)
    x[1] -= 0.15
    u = u_sol[ind] - K_tv[ind] * (x - x_sol[ind])
    u = clamp.(u, -0.2, +0.2)
    u = [zeros(6); u[7:end]] / S
    set_input!(m, SVector{nu}(u))
end

x1_roll = deepcopy(xref[1])
x1_roll[3] += 0.3
x1_roll[4] += 0.1
set_minimal_state!(mech_sim, x1_roll)

Main.@elapsed storage = simulate!(mech_sim, 3.0, ctrl!,
    record=true,
    opts=SolverOptions(rtol=1e-5, btol=1e-4, undercut=5.0),
    )
Dojo.visualize(mech_sim, storage, vis=vis, build=false)

















maximum(abs.(vcat(u_sol...)))

A = [zeros(n,n) for i = 1:T-1]
B = [zeros(n,m) for i = 1:T-1]
for i = 1:T-1
    dynamics_jacobian_state(A[i], env, deepcopy(x_sol)[i], deepcopy(u_sol)[i], zeros(0))
    dynamics_jacobian_input(B[i], env, deepcopy(x_sol)[i], deepcopy(u_sol)[i], zeros(0))
end
qt = [0.1; 0.005; 0.005;
    1e-1 * ones(3);
    1e-9 * ones(3);
    1e-9 * ones(3);
    fill([2, 1e-9], 12)...]

rt = timestep * 100.0 * ones(m)

Q = [Diagonal(qt) for i = 1:T]
R = [Diagonal(rt) for i = 1:T-1]

K_tv, P_tv = tvlqr(A,B,Q,R)

################################################################################
# Simulate accurately the TVLQR policy
################################################################################
# divide timestep
S = 5
mech_sim = get_mechanism(:quadruped,
    contact_body=true,
    timestep=timestep/S,
    gravity=gravity,
    friction_coefficient=friction_coefficient,
    damper=damper,
    spring=spring)

# controller
function ctrl!(m, k; x_sol=x_sol, u_sol=u_sol, K_sol=K_sol)
    N = length(u_sol)
    nu = input_dimension(m)
    @show Int(floor((k-1)/S))
    ind = (Int(floor((k-1)/S))) % N + 1

    x = get_minimal_state(m)
    x[1] -= 0.15
    u = u_sol[ind] - K_tv[ind] * (x - deepcopy(x_sol)[ind])
    u = clamp.(u, -2, +2)
    u = [zeros(6); u[7:end]] / S
    set_input!(m, SVector{nu}(u))
end

x1_roll = deepcopy(xref[1])
x1_roll[3] += 0.0
x1_roll[4] += 0.0
set_minimal_state!(mech_sim, x1_roll)

Main.@elapsed storage = simulate!(mech_sim, 3.0, ctrl!,
    record=true,
    opts=SolverOptions(rtol=1e-5, btol=1e-4, undercut=5.0),
    )
Dojo.visualize(mech_sim, storage, vis=vis, build=false)

convert_frames_to_video_and_gif("quadruped_tvlqr_forward")
