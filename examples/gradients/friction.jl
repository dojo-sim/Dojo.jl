# Utils
function module_dir()
    return joinpath(@__DIR__, "..", "..")
end

# Activate package
using Pkg
Pkg.activate(module_dir())

using MeshCat

# Open visualizer
vis = Visualizer()
open(vis)

# Include new files
include(joinpath(module_dir(), "examples", "loader.jl"))


mech = getmechanism(:box2d, Δt=0.01, g=-9.81, cf=2.0, contact=true, conetype=:soc);
initialize!(mech, :box2d, x=[0,0.], v=[0,0.], θ=0.0, ω=0.0)

# function controller!(mechanism, k)
#     setControl!(mech, [0.0,0,0])
#     return
# end
#
# storage = simulate!(mech, 3.00, controller!, record = true, verbose = true,
#     opts=InteriorPointOptions(rtol=3e-4, btol=3e-4))
# visualize(mech, storage, vis = vis)

mech = getmechanism(:box2d, Δt=0.01, g=-10.0, cf=0.5, contact=true, conetype=:soc);
initialize!(mech, :box2d, x=[0,0.], v=[0,0.], θ=0.0, ω=0.0)
z0 = getMaxState(mech)
nu = controldim(mech)
u0 = zeros(nu)
lift_force = total_mass(mech) * 10.0
u0 = [0.0, 0.1, 0.0]
opts_grad = InteriorPointOptions(rtol=1e-10, btol=1e-10)
∇x, ∇u = getMinGradients!(mech, z0, u0, opts=opts_grad)
∇u





function box2d_displacement(mechanism::Mechanism, F; btol=1e-10, mode::Symbol=:friction)
    (mode == :friction) && (idx = 1)
    (mode == :impact) && (idx = 2)
    Δt = mechanism.Δt
    opts_grad = InteriorPointOptions(rtol=1e-10, btol=btol, undercut=1.5)

    initialize!(mechanism, :box2d, x=[0,0.], v=[0,0.], θ=0.0, ω=0.0)
    z0 = getMaxState(mechanism)
    u0 = zeros(3)
    u0[idx] += F * Δt

    ∇x, ∇u = getMinGradients!(mechanism, z0, u0, opts=opts_grad)
    GB0 = gradient_bundle_0th(mechanism, z0, u0, opts=opts_grad)
    GB1 = gradient_bundle_1st(mechanism, z0, u0, opts=opts_grad)
    step!(mechanism, z0, u0, opts=opts_grad)
    z1 = getNextState(mechanism)
    x0 = max2min(mechanism, z0)
    x1 = max2min(mechanism, z1)

    ∂x∂F = ∇u[idx,idx]/Δt
    GB0xF = GB0[idx,idx]/Δt
    GB1xF = GB1[idx,idx]/Δt
    Δx = x1[idx] - x0[idx]
    return Δx, ∂x∂F, GB0xF, GB1xF
end

mech = getmechanism(:box2d, Δt=0.05, g=-10.0, cf=1.0, contact=true, conetype=:linear);
btol = 1e-4
Fs = 0:0.05:20
plt = plot(layout=(2,1))
mode = :friction
# mode = :impact

for btol = 1e-10
    x  = [box2d_displacement(mech, F, btol=btol, mode=mode)[1] for F = Fs]
    ∇x = [box2d_displacement(mech, F, btol=btol, mode=mode)[2] for F = Fs]
    plot!(plt[1,1], Fs,  x, legend=false, xlabel="F", ylabel="x", linewidth=3.0)
    plot!(plt[2,1], Fs, ∇x, legend=false, xlabel="F", ylabel="∇x", linewidth=3.0)
    display(plt)
end

# for btol in exp.(log(10)* (-10:1:-3))
for btol = 3e-4
    x  = [box2d_displacement(mech, F, btol=btol, mode=mode)[1] for F = Fs]
    ∇x = [box2d_displacement(mech, F, btol=btol, mode=mode)[2] for F = Fs]
    # plot!(plt[1,1], Fs,  x, legend=false, xlabel="F", ylabel="x", linewidth=3.0, linestyle=:dot)
    plot!(plt[2,1], Fs, ∇x, legend=false, xlabel="F", ylabel="∇x", linewidth=3.0, linestyle=:dot)
    display(plt)
end

display(plt)


function normal_sample(μ, Σ)
    n = length(μ)
    x = μ + sqrt(Σ) * randn(n)
    return x
end

function gradient_bundle_1st(mechanism::Mechanism, z0, u0; N::Int=100, Σ=1e-6I,
        opts=InteriorPointOptions(rtol=1e-10, btol=1e-10, undercut=1.5))
    nx = 6
    nu = 3
    ∇ = zeros(nx,nu)

    for i = 1:N
        u = normal_sample(u0, Σ)
        ∇x, ∇u = getMinGradients!(mechanism, z0, u, opts=opts)
        ∇ += ∇u
    end
    return ∇ ./ N
end

function gradient_bundle_0th(mechanism::Mechanism, z0, u0; N::Int=100, Σ=1e-6*I,
        opts=InteriorPointOptions(rtol=1e-10, btol=1e-10, undercut=1.5))
    nx = 6
    nu = 3
    δx = [zeros(nx) for i=1:N]
    δu = [zeros(nu) for i=1:N]

    step!(mechanism, z0, u0, opts=opts)
    z10 = getNextState(mechanism)
    x10 = max2min(mechanism, z10)

    for i = 1:N
        u = normal_sample(u0, Σ)
        step!(mechanism, z0, u, opts=opts)
        z1 = getNextState(mechanism)
        x1 = max2min(mechanism, z1)
        δu[i] = u - u0
        δx[i] = x1 - x10
    end
    ∇ = leastsquares(δx, δu)
    return ∇
end

∇1 = gradient_bundle_1st(mech, z0, u0, N=100, Σ=1e-6I)
∇0 = gradient_bundle_0th(mech, z0, u0, N=100, Σ=1e-6I)


using Symbolics

nx = 6
nu = 3
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

δu1 = [rand(3) for k=1:10]
B1 = rand(6,3)
Bf1 = vec(B1)
δx1 = [B1*δu for δu in δu1]
δx1 = [rand(6) for k=1:10]

leastsquares(δx1, δu1)




# function integrate_displacement(f0::T, ∇f::Vector{T}, x) where T
#     n = length(x)
#     f = [f0; zeros(n-1)]
#     for k = 2:n
#         f[k] = f[k-1] + (∇f[k-1] + ∇f[k])/2 * (x[k] - x[k-1])
#     end
#     return f
# end
