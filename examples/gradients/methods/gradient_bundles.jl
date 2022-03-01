function gradient_bundle_1st(mechanism::Mechanism, z0, u0, idx;
    N::Int=100, 
    Σ=1e-6I,
    opts=SolverOptions(rtol=1e-10, btol=1e-10, undercut=1.0, no_progress_undercut=1.0))

    nx = 6
    nu = 3
    ∇ = zeros(nx,nu)

    for i = 1:N
        # u = normal_sample(u0, Σ)
        u = zeros(3)
        u[idx:idx] = normal_sample(u0[idx:idx], Σ)
        ∇x, ∇u = get_minimal_gradients!(mechanism, z0, u, 
            opts=opts)
        # ∇u = (u[idx] / mechanism.timestep >= 10.0) * ones(nx,nu) * mechanism.timestep
        ∇ += ∇u
    end
    return ∇ ./ N
end

function gradient_bundle_0th(mechanism::Mechanism, z0, u0, idx; 
    N::Int=100, 
    Σ=1e-6*I,
    opts=SolverOptions(rtol=1e-10, btol=1e-10, undercut=1.0, no_progress_undercut=1.0))
    nx = 6
    nu = 1
    δx = [zeros(nx) for i=1:N]
    δu = [zeros(nu) for i=1:N]

    step!(mechanism, z0, u0, opts=opts)
    z10 = get_next_state(mechanism)
    x10 = maximal_to_minimal(mechanism, z10)

    for i = 1:N
        u = zeros(3)
        u[idx:idx] = normal_sample(u0[idx:idx], Σ)
        step!(mechanism, z0, u, opts=opts)
        z1 = get_next_state(mechanism)
        x1 = maximal_to_minimal(mechanism, z1)
        # x11 = copy(x10)
        # x11[idx] = max(0, u[idx]  - 10.0) / 10
        # x11[3+idx] = max(0, u[idx - 10.0)
        δu[i] = (u - u0)[idx:idx]
        δx[i] = x1 - x10
    end
    ∇ = zeros(nx,3)
    ∇[:,idx:idx] = leastsquares(δx, δu)
    return ∇
end
