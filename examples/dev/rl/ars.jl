# Utils
function module_dir()
    return joinpath(@__DIR__, "..", "..")
end

# Activate package
using Pkg
Pkg.activate(module_dir())

# Load packages
using MeshCat

# Open visualizer
vis = Visualizer()
open(vis)

# Include new files
include(joinpath(module_dir(), "examples", "loader.jl"))

Δt0 = 0.05
g0 = -0.81
u0 = [0.0; 0.0; g0 * Δt0]
u_mask = [0 0 0 1 0 0 0;
          0 0 0 0 1 0 0;
          0 0 0 0 0 0 1]
nxm = 26
nx = 15
nu = 3
K0 = zeros(nu,nx)

mech = getmechanism(:hopper, Δt = Δt0, g = g0, spring = 0.0, damper = 10.0)
initialize!(mech, :hopper, leg_length_nominal = 0.5, altitude = 0.1)
function controller!(mechanism, k, K = zeros(3,26))
    xm = getMaxState(mechanism)
    u = K * xm
    um = u_mask' * u
    off = 0
    for (i,eqc) in enumerate(collect(mechanism.eqconstraints))
        nu = controldim(eqc)
        setForce!(mechanism, eqc, SVector{nu}(um[off .+ (1:nu)]))
        off += nu
    end
    return
end
@elapsed storage = simulate!(mech, 4.00, controller!, record = true, solver = :mehrotra!, verbose = false)
visualize(mech, storage, vis = vis)



function evaluate_policy(mechanism::Mechanism, K::AbstractMatrix, μ::AbstractVector, Σ; H = 1.05, btol::T = 1e-4,
        vis::Visualizer = vis, storage_ref::Storage = storage_ref, display::Bool = false) where {T}
    function controller!(mechanism, k)
        xm = getMaxState(mechanism)
        x = max2min(mechanism, xm)
        u = Δt0 * K * Diagonal(1 ./ sqrt.(diag(Σ))) * (x - μ)
        u = clamp.(u, -0.5Δt0, 0.5Δt0)
        um = u_mask' * u
        off = 0
        for (i,eqc) in enumerate(collect(mechanism.eqconstraints))
            nu = controldim(eqc)
            setForce!(mechanism, eqc, SVector{nu}(um[off .+ (1:nu)]))
            off += nu
        end
        return
    end

    initialize!(mechanism, :hopper, leg_length_nominal = 0.5, altitude = 0.0)
    storage = simulate!(mechanism, H, controller!, record = true, solver = :mehrotra!, btol = btol, verbose = false)
    display && visualize(mechanism, storage, vis = vis)

    cost, μt, Σt = cost_function(mechanism, storage, storage_ref)
    μ = μ + 0.01 * μt
    Σ = Σ + 0.01 * Σt
    return cost, μ, Σ
end

function cost_function(mechanism::Mechanism, storage::Storage{T,N}, storage_ref::Storage{T,Nr};
        Qx = Diagonal(ones(3)), Qv = Diagonal(1e-3ones(3)), Qq = Diagonal(1e-2ones(3)), Qϕ = Diagonal(1e-3ones(3))) where {T,N,Nr}
    Nb = length(storage.x)
    c = 0.0
    z = []
    for t = 1:N
        zt = zeros(0)
        for i = 1:Nb
            m = 1 + (i==1) * 10
            Δx = storage.x[i][t] - storage_ref.x[i][Nr]
            Δv = storage.v[i][t] - storage_ref.v[i][Nr]
            Δq = storage.q[i][t] * inv(storage_ref.q[i][Nr])
            Δq = [Δq.x, Δq.y, Δq.z]
            Δϕ = storage.ω[i][t] - storage_ref.ω[i][Nr]
            c += 0.5 * Δx' * Qx * Δx * m
            c += 0.5 * Δv' * Qv * Δv
            c += 0.5 * Δq' * Qq * Δq
            c += 0.5 * Δϕ' * Qϕ * Δϕ
            push!(zt, storage.x[i][t]...)
            push!(zt, storage.v[i][t]...)
            push!(zt, vector(storage.q[i][t])...)
            push!(zt, storage.ω[i][t]...)
        end
        push!(z, zt)
    end
    x = [max2min(mechanism, zt) for zt in z]
    μ = mean(x)
    Σ = cov(x)
    return c, μ, Σ
end

function ars(mechanism::Mechanism{T}; storage_ref::Storage = storage_ref, H::T = 2.05, btol::T = 1e-4,
        α::T = 0.01, ν::T = 0.001, Nsample::Int = 4, iter::Int = 30, vis::Visualizer = vis) where {T}
    K = 0.001 * (rand(nu,nx) .- 0.5)
    μ = zeros(nx)
    Σ = Matrix(Diagonal(ones(nx)))
    for i = 1:iter
        Δ = zeros(nu, nx)
        for s = 1:Nsample
            δK = randn(nu,nx)
            rp, μ, Σ = evaluate_policy(mechanism, K + ν*δK, μ, Σ; H = H, btol = btol, vis = vis, storage_ref = storage_ref, display = false)
            rm, μ, Σ = evaluate_policy(mechanism, K - ν*δK, μ, Σ; H = H, btol = btol, vis = vis, storage_ref = storage_ref, display = false)
            Δ += (rp - rm) * δK
            @show rp
            # @show norm(K)
        end
        K = K - α/Nsample * Δ
        (i%10 == 0) && evaluate_policy(mechanism, K, μ, Σ; H = H, vis = vis, storage_ref = storage_ref, display = true)
    end
    return K
end

r0 = evaluate_policy(mech, K0, display = true, btol = 1e-2)
Kopt = ars(mech, btol = 1e-3)
