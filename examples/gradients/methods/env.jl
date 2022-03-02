function box2D_dojo(mechanism::Mechanism, F; 
    rtol=1e-10, 
    btol=1e-10, 
    undercut=1.0, 
    no_progress_undercut=1.0, 
    mode::Symbol=:friction)

    (mode == :friction) && (idx = 1)
    (mode == :impact) && (idx = 2)
    timestep= mechanism.timestep
    opts_grad = SolverOptions(
        rtol=rtol, 
        btol=btol, 
        undercut=undercut, 
        no_progress_undercut=no_progress_undercut, 
        verbose=false)

    initialize!(mechanism, :box2D, 
        x=[0.0, 0.0], 
        v=[0.0, 0.0], 
        θ=0.0, 
        ω=0.0)
    z0 = get_maximal_state(mechanism)
    u0 = zeros(3)
    u0[idx] += F 

    ∇x, ∇u = get_minimal_gradients!(mechanism, z0, u0, 
        opts=opts_grad)

    step!(mechanism, z0, u0, opts=opts_grad)

    z1 = get_next_state(mechanism)
    x0 = maximal_to_minimal(mechanism, z0)
    x1 = maximal_to_minimal(mechanism, z1)

    ∂x∂F = ∇u[idx,idx]
    Δx = x1[idx] - x0[idx]
    return Δx, ∂x∂F
end

function box2D_gradientbundle(mechanism::Mechanism, F; 
    N::Int=100, 
    Σ=1e-6*I, 
    rtol=1e-10, 
    btol=1e-10, 
    undercut=1.5, 
    no_progress_undercut=10.0, 
    mode::Symbol=:friction)

    (mode == :friction) && (idx = 1)
    (mode == :impact) && (idx = 2)
    timestep= mechanism.timestep
    opts_grad = SolverOptions(
            rtol=rtol, 
            btol=btol, 
            undercut=undercut, 
            no_progress_undercut=no_progress_undercut)

    initialize!(mechanism, :box2D, 
        x=[0.0, 0.0], 
        v=[0.0, 0.0], 
        θ=0.0, 
        ω=0.0)
    z0 = get_maximal_state(mechanism)
    u0 = zeros(3)
    u0[idx] += F

    GB0 = gradient_bundle_0th(mechanism, z0, u0, idx, 
        opts=opts_grad, 
        N=N, 
        Σ=Σ)
    GB1 = gradient_bundle_1st(mechanism, z0, u0, idx, 
        opts=opts_grad, 
        N=N, 
        Σ=Σ)
    step!(mechanism, z0, u0, 
        opts=opts_grad)
    z1 = get_next_state(mechanism)
    x0 = maximal_to_minimal(mechanism, z0)
    x1 = maximal_to_minimal(mechanism, z1)

    GB0xF = GB0[idx,idx]
    GB1xF = GB1[idx,idx]
    Δx = x1[idx] - x0[idx]
    return Δx, GB0xF, GB1xF
end