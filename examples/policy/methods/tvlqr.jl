function tvlqr(A, B, Q, R)
    T = length(Q)

    P = [zero(A[1]) for t = 1:T]
    K = [zero(B[1]') for t = 1:T-1]
    P[T] = Q[T]

    for t = T-1:-1:1
        K[t] = (R[t] + B[t]' * P[t+1] *  B[t]) \ (B[t]' * P[t+1] * A[t])
        P[t] = (Q[t] + K[t]' * R[t] * K[t]
                + (A[t] - B[t] * K[t])' * P[t+1] * (A[t] - B[t] * K[t]))
    end

    return K, P
end

function tvlqr(x, u, env;
        q_tracking=[0.3; 0.05; 0.05;
            5e-1 * ones(3);
            1e-6 * ones(3);
            1e-6 * ones(3);
            fill([2, 1e-6], 12)...],
        r_tracking=env.mechanism.timestep * 100 * ones(length(u[1])))
    nx = length(x[1])
    nu = length(u[1])

    A = [zeros(nx,nx) for i = 1:T-1]
    B = [zeros(nx,nu) for i = 1:T-1]
    for i = 1:T-1
        dynamics_jacobian_state(A[i], env, x[i], u[i], zeros(0))
        dynamics_jacobian_input(B[i], env, x[i], u[i], zeros(0))
    end

    Q = [Diagonal(q_tracking) for i = 1:T]
    R = [Diagonal(r_tracking) for i = 1:T-1]
    K_tv, P_tv = tvlqr(A,B,Q,R)
    return K_tv, P_tv
end
