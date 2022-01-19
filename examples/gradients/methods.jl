using Symbolics

nx = 6
nu = 1
@variables δx[1:nx] δu[1:nu] Bf[1:nx*nu]

function cost(δx, δu, Bf)
    B = reshape(Bf, nx, nu)
    return 1/2 * (δx - B*δu)' * (δx - B*δu)
end

c = cost(δx, δu, Bf)
dc = Symbolics.gradient(c, Bf)
ddc = Symbolics.hessian(c, Bf)

c_fct = eval(Symbolics.build_function([c], δx, δu, Bf)[1])
dc_fct = eval(Symbolics.build_function(dc, δx, δu, Bf)[1])
ddc_fct = eval(Symbolics.build_function(ddc, δx, δu, Bf)[1])

δx0 = rand(nx)
δu0 = rand(nu)
Bf0 = rand(nx*nu)

c_fct(δx0, δu0, Bf0)
dc_fct(δx0, δu0, Bf0)
ddc_fct(δx0, δu0, Bf0)

function leastsquares(δx, δu)
    N = length(δx)
    nx = length(δx[1])
    nu = length(δu[1])
    Bf = zeros(nx*nu)
    e = 0.0
    grad = zeros(nx*nu)
    hess = zeros(nx*nu, nx*nu)
    for i = 1:N
        e += c_fct(δx[i], δu[i], Bf)[1]
        grad += dc_fct(δx[i], δu[i], Bf)
        hess += ddc_fct(δx[i], δu[i], Bf)
    end
    Bf = Bf - hess \ grad
    B = reshape(Bf, (nx, nu))
    return B
end

function normal_sample(μ, Σ)
    n = length(μ)
    x = μ + sqrt(Σ) * randn(n)
    return x
end

function gradient_bundle_1st(mechanism::Mechanism, z0, u0, idx; N::Int=100, Σ=1e-6I,
        opts=InteriorPointOptions(rtol=1e-10, btol=1e-10, undercut=1.0, no_progress_undercut=1.0))
    nx = 6
    nu = 3
    ∇ = zeros(nx,nu)

    for i = 1:N
        # u = normal_sample(u0, Σ)
        u = zeros(3)
        u[idx:idx] = normal_sample(u0[idx:idx], Σ)
        ∇x, ∇u = getMinGradients!(mechanism, z0, u, opts=opts)
        # ∇u = (u[idx] / mechanism.Δt >= 10.0) * ones(nx,nu) * mechanism.Δt
        ∇ += ∇u
    end
    return ∇ ./ N
end

function gradient_bundle_0th(mechanism::Mechanism, z0, u0, idx; N::Int=100, Σ=1e-6*I,
        opts=InteriorPointOptions(rtol=1e-10, btol=1e-10, undercut=1.0, no_progress_undercut=1.0))
    nx = 6
    nu = 1
    δx = [zeros(nx) for i=1:N]
    δu = [zeros(nu) for i=1:N]

    step!(mechanism, z0, u0, opts=opts)
    z10 = getNextState(mechanism)
    x10 = max2min(mechanism, z10)

    for i = 1:N
        u = zeros(3)
        u[idx:idx] = normal_sample(u0[idx:idx], Σ)
        step!(mechanism, z0, u, opts=opts)
        z1 = getNextState(mechanism)
        x1 = max2min(mechanism, z1)
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

function box2d_dojo(mechanism::Mechanism, F; rtol=1e-10, btol=1e-10, undercut=1.0, no_progress_undercut=1.0, mode::Symbol=:friction)
    (mode == :friction) && (idx = 1)
    (mode == :impact) && (idx = 2)
    Δt = mechanism.Δt
    opts_grad = InteriorPointOptions(rtol=rtol, btol=btol, undercut=undercut, no_progress_undercut=no_progress_undercut, verbose=false)

    initialize!(mechanism, :box2d, x=[0,0.], v=[0,0.], θ=0.0, ω=0.0)
    z0 = getMaxState(mechanism)
    u0 = zeros(3)
    u0[idx] += F 

    ∇x, ∇u = getMinGradients!(mechanism, z0, u0, opts=opts_grad)

    step!(mechanism, z0, u0, opts=opts_grad)

    z1 = getNextState(mechanism)
    x0 = max2min(mechanism, z0)
    x1 = max2min(mechanism, z1)

    ∂x∂F = ∇u[idx,idx]
    Δx = x1[idx] - x0[idx]
    return Δx, ∂x∂F
end

function box2d_gradientbundle(mechanism::Mechanism, F; N::Int=100, Σ=1e-6*I, rtol=1e-10, btol=1e-10, undercut=1.5, no_progress_undercut=10.0, mode::Symbol=:friction)
    (mode == :friction) && (idx = 1)
    (mode == :impact) && (idx = 2)
    Δt = mechanism.Δt
    opts_grad = InteriorPointOptions(rtol=rtol, btol=btol, undercut=undercut, no_progress_undercut=no_progress_undercut)

    initialize!(mechanism, :box2d, x=[0,0.], v=[0,0.], θ=0.0, ω=0.0)
    z0 = getMaxState(mechanism)
    u0 = zeros(3)
    u0[idx] += F

    GB0 = gradient_bundle_0th(mechanism, z0, u0, idx, opts=opts_grad, N=N, Σ=Σ)
    GB1 = gradient_bundle_1st(mechanism, z0, u0, idx, opts=opts_grad, N=N, Σ=Σ)
    step!(mechanism, z0, u0, opts=opts_grad)
    z1 = getNextState(mechanism)
    x0 = max2min(mechanism, z0)
    x1 = max2min(mechanism, z1)

    GB0xF = GB0[idx,idx]
    GB1xF = GB1[idx,idx]
    Δx = x1[idx] - x0[idx]
    return Δx, GB0xF, GB1xF
end

function write_csv(names, data, filename)
    open(joinpath(@__DIR__, filename), "w") do f
        write(f, join(names, ",") * "\n")
        for tup in data
            write(f, join(tup, ",") * "\n")
        end
    end
    return nothing
end
