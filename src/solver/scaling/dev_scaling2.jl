using LinearAlgebra
using Test
using Plots
using Random

include("util.jl")

################################################################################
# Quadratic Cone Program
################################################################################

nx = 6
ns = 3
nb = 2

Random.seed!(10)
P = rand(nx, nx)
P = P'*P
c = rand(nx)
D = rand(ns, nx)
h = rand(ns)
A = rand(nb, nx)
b = rand(nb)


################################################################################
# Objective and Constraints
################################################################################
Jmat(n) = Diagonal([1; -ones(n-1)])
val_soc(s) = s[1]^2 - s[2:end]'*s[2:end]
cone_neutral(n) = [1; zeros(n-1)]
cone_prod(s, z) = [s'*z; s[1]*z[2:end] + z[1]*s[2:end]]
function cone_mat(s)
    n = length(s)
    mat = Matrix(Diagonal(s[1]*ones(n)))
    mat[1,2:end] .= s[2:end]
    mat[2:end,1] .= s[2:end]
    return mat
end
function soc_step_length(λ::AbstractVector{T}, Δ::AbstractVector{T};
        τ::T = 0.99, ϵ::T = 1e-14) where {T}
    # check Section 8.2 CVXOPT
    # The CVXOPT linear and quadratic cone program solvers

    # Adding to slack ϵ to make sure that we never get out of the cone
    λ0 = λ[1]
    λ_λ = max(λ0^2 - λ[2:end]' * λ[2:end], 1e-25)
    if λ_λ < 0.0
        @show λ_λ
        @warn "should always be positive"
    end
    λ_λ += ϵ
    λ_Δ = λ0 * Δ[1] - λ[2:end]' * Δ[2:end] + ϵ

    ρs = λ_Δ / λ_λ
    ρv = Δ[2:end] / sqrt(λ_λ)
    ρv -= (λ_Δ / sqrt(λ_λ) + Δ[1]) / (λ0 / sqrt(λ_λ) + 1) * λ[2:end] / λ_λ
    # we make sre that the inverse always exists with ϵ,
    # if norm(ρv) - ρs) is negative (Δ is pushing towards a more positive cone)
        # the computation is ignored and we get the maximum value for α = 1.0
    # else we have α = τ / norm(ρv) - ρs)
    # we add ϵ to the denumerator to ensure strict positivity and avoid 1e-16 errors.
    α = 1.0
    if norm(ρv) - ρs > 0.0
        α = min(α, τ / (norm(ρv) - ρs))
    end
    return α
end

function soc_nt_scaling(s::AbstractVector{T}, z::AbstractVector{T}) where {T}
    # http://www.seas.ucla.edu/~vandenbe/publications/coneprog.pdf
    n = length(s)
    @assert length(z) == n
    J = Jmat(n)
    e = cone_neutral(n)

    z̄ = 1 / sqrt(z'*J*z) * z
    s̄ = 1 / sqrt(s'*J*s) * s
    γ = sqrt((1 + z̄'*s̄) / 2)
    ωb = 1 / (2γ) * (s̄ + J*z̄)

    v = 1 / sqrt(2 * (ωb[1] + 1)) * (ωb + e)
    W̄ = 2 * v * v' - J
    W̄i = 2 * J * v * v' * J - J

    r = ((s'*J*s)/(z'*J*z))^0.25
    W = r * W̄
    Wi = 1/r * W̄i

    λ = W * z
    # @assert norm(W̄i - inv(W̄), Inf) < 1e-8
    # @assert norm(Wi - inv(W), Inf) < 1e-8
    # @assert norm(λ - inv(W)*s) < 1e-8
    return W, Wi, λ
end

function ort_nt_scaling(s::AbstractVector{T}, z::AbstractVector{T}) where {T}
    # http://www.seas.ucla.edu/~vandenbe/publications/coneprog.pdf
    W = Diagonal(sqrt.(s) ./ sqrt.(z))
    Wi = Diagonal(sqrt.(z) ./ sqrt.(s))
    λ = W*z
    return W, Wi, λ
end

function unpack(Δ)
    Δx = Δ[1:nx]
    Δy = Δ[nx+1:nx+nb]
    Δz = Δ[nx+nb+1:nx+nb+ns]
    Δs = Δ[nx+nb+ns+1:nx+nb+2ns]
    return Δx, Δy, Δz, Δs
end

function jacobian(s, z; scaling = true)
    W, Wi, λ = soc_nt_scaling(s, z)
    mλ = cone_mat(λ)
    jac = [P A' D' zeros(nx,ns);
         A zeros(nb,nb) zeros(nb,ns) zeros(nb,ns);
         D zeros(ns,nb) zeros(ns,ns) Diagonal(ones(ns));
         zeros(ns,nx) zeros(ns,nb) mλ*W    mλ*Wi';
    ]
    ujac = [P A' D' zeros(nx,ns);
        A zeros(nb,nb) zeros(nb,ns) zeros(nb,ns);
        D zeros(ns,nb) zeros(ns,ns) Diagonal(ones(ns));
        zeros(ns,nx) zeros(ns,nb) cone_mat(s) cone_mat(z);
        ]
    scaling && return jac
    return ujac
end

function residual(x, y, z, s)
    jac = jacobian(s, z)
    res = jac[1:nx+nb+ns, 1:nx+nb+ns] * [x; y; z] + [c; -b; -h] + [zeros(nx+nb); s]
    return res
end

function affine(s, z, res; scaling = true)
    jac = jacobian(s, z, scaling = scaling)
    W, Wi, λ = soc_nt_scaling(s, z)
    Δa = - jac \ [res; cone_prod(λ, λ)]
    uΔa = - jac \ [res; cone_prod(s, z)]
    scaling && return Δa
    return uΔa
end

function correction(s, z, Δsa, Δza, σ, μ, res; scaling = true, αsa = 1.0, αza = 1.0)
    jac = jacobian(s, z, scaling = scaling)
    n = length(s)
    e = cone_neutral(n)
    W, Wi, λ = soc_nt_scaling(s, z)
    Δ = - jac \ [(1-σ)*res; cone_prod(λ, λ) + cone_prod(Wi'*αsa*Δsa, αza*W*Δza) - σ*μ*e]
    uΔ = - jac \ [(1-σ)*res; cone_prod(s, z) + cone_prod(αsa*Δsa, αza*Δza) - σ*μ*e]
    scaling && return Δ
    return uΔ
end

function qp_solve(x, y, z, s; scaling = true, iter::Int = 10,
        btol = 1e-8, rtol = 1e-8, μ_mult = 1.0)
    m = 1
    n = length(s)
    e = cone_neutral(n)

    α = 0.0
    μ = 0.0
    σ = 0.0
    Δ = zeros(nx+nb+2ns)
    println("#########################################")
    for i = 1:iter
        # Residual
        res = residual(x, y, z, s)
        μ = s'*z / m * μ_mult
        W, Wi, λ = soc_nt_scaling(s, z)

        # Termination
        cvio_s = min(val_soc(s), 0)
        cvio_z = min(val_soc(z), 0)
        cvio = max(cvio_s, cvio_z)
        bvio = norm(z'*s, Inf)
        rvio = norm(res, Inf)
        println(i, "   bvio", scn(bvio), "   rvio", scn(rvio),
            "   α", scn(α), "   μ", scn(μ), "   σ", scn(σ),
            "   Δ∞", scn(norm(Δ, Inf)), "   vs", scn(val_soc(s)),
            "   vz", scn(val_soc(z)))
        if bvio < btol && rvio < rtol
            return x, y, z, s
        end

        # Affine direction
        Δa = affine(s, z, res, scaling = scaling)
        Δxa, Δya, Δza, Δsa = unpack(Δa)

        # Step size and Centering
        if false && scaling
            αsa = soc_step_length(λ, Wi'*Δsa, τ = 1.0)
            αza = soc_step_length(λ, W*Δza, τ = 1.0)
        else
            αsa = soc_step_length(s, Δsa, τ = 1.0)
            αza = soc_step_length(z, Δza, τ = 1.0)
        end
        # αsa0 = soc_step_length(λ, Wi'*Δsa, τ = 1.00)
        # αza0 = soc_step_length(λ, W*Δza, τ = 1.00)
        # αsa1 = soc_step_length(s, Δsa, τ = 1.00)
        # αza1 = soc_step_length(z, Δza, τ = 1.00)
        # @show norm(αsa0 - αsa1)
        # @show norm(αza0 - αza1)
        αa = min(αsa, αza)
        # αa = max(αsa, αza)
        ρ = (s + αa*Δsa)' * (z + αa*Δza) / (s'*z)
        # ρ = (s + αsa*Δsa)' * (z + αza*Δza) / (s'*z)
        σ = clamp(ρ, 0, 1)^3

        # Correction
        Δ = correction(s, z, Δsa, Δza, σ, μ, res, scaling = scaling)
        Δx, Δy, Δz, Δs = unpack(Δ)

        # Update iterates
        if false && scaling
            αs = soc_step_length(λ, Wi'*Δs, τ = 0.99)
            αz = soc_step_length(λ, W*Δz, τ = 0.99)
        else
            αs = soc_step_length(s, Δs, τ = 0.99)
            αz = soc_step_length(z, Δz, τ = 0.99)
        end
        # αs0 = soc_step_length(λ, Wi'*Δs, τ = 0.99)
        # αz0 = soc_step_length(λ, W*Δz, τ = 0.99)
        # αs1 = soc_step_length(s, Δs, τ = 0.99)
        # αz1 = soc_step_length(z, Δz, τ = 0.99)
        # @show norm(αs0 - αs1)
        # @show norm(αz0 - αz1)

        α = min(αs, αz)
        # α = max(αs, αz)

        # Display cones
        plt = plot(layout = (2,2))
        plot_cone(s, αsa*Δsa, plt = plt[1,1], linewidth = 2.0, show = false, xlabel = "S")
        plot_cone(s, α*Δs,    plt = plt[1,1], linewidth = 4.0, show = false, xlabel = "S")
        plot_cone(z, αza*Δza, plt = plt[1,2], linewidth = 2.0, show = false, xlabel = "Z")
        plot_cone(z, α*Δz,    plt = plt[1,2], linewidth = 4.0, show = false, xlabel = "Z")
        plot_cone(s, αsa*Δsa, plt = plt[2,1], linewidth = 2.0, show = false, xlabel = "S", zoom = true)
        plot_cone(s, α*Δs,    plt = plt[2,1], linewidth = 4.0, show = false, xlabel = "S", zoom = true)
        plot_cone(z, αza*Δza, plt = plt[2,2], linewidth = 2.0, show = false, xlabel = "Z", zoom = true)
        plot_cone(z, α*Δz,    plt = plt[2,2], linewidth = 4.0, show = false, xlabel = "Z", zoom = true)
        # plot_cone(s, 1e0*Δsa, plt = plt[1,1], linewidth = 2.0, show = false, xlabel = "S")
        # plot_cone(s, 1e0*Δs,   plt = plt[1,1], linewidth = 4.0, show = false, xlabel = "S")
        # plot_cone(z, 1e0*Δza, plt = plt[1,2], linewidth = 2.0, show = false, xlabel = "Z")
        # plot_cone(z, 1e0*Δz,   plt = plt[1,2], linewidth = 4.0, show = false, xlabel = "Z")
        display(plt)
        sleep(0.15)
        x = x + α*Δx
        y = y + α*Δy
        z = z + α*Δz
        # z = z + αz*Δz
        s = s + α*Δs
        # s = s + αs*Δs
    end

    return x, y, z, s
end

Random.seed!(100)
x0 = rand(nx)
y0 = rand(nb)
z0 = [1.0, -(1.0-1e-10), 0.0]
s0 = [1.0, -(1.0-1e-10), 0.0]
z0 = [1e-0, +0.1, -0.1]
s0 = [1e-0, -0.1, +0.1]
θ = π/10
R = [cos(θ) sin(θ); -sin(θ) cos(θ)]
s0 = [ss[1]; R*ss[2:3]]
z0 = [zs[1]; R*zs[2:3]]
qp_solve(x0, y0, z0, s0, iter = 30, scaling = true, μ_mult = 1.0)
qp_solve(x0, y0, z0, s0, iter = 30, scaling = false, μ_mult = 1.0)
# xs, ys, zs, ss = qp_solve(x0, y0, z0, s0, iter = 30, scaling = true)
# h -= ss
# Ws = nd_scaling(ss, zs)
# Ws * zs - inv(Ws) * ss
# zs[1]
# ss[1]
val_soc(ss)
val_soc(zs)

################################################################################
# Test
################################################################################
@test norm(Jmat(3) - Diagonal([1, -1, -1]), Inf) < 1e-8
@test norm(cone_neutral(3) - [1, 0, 0], Inf) < 1e-8
s0 = [3, 0, 1.]
z0 = [3, 2, 1.]
@test norm(val_soc(s0) - 8, Inf) < 1e-8
@test norm(cone_prod(s0, z0) - [10, 6, 6], Inf) < 1e-8
@test norm(cone_prod(s0, z0) - cone_prod(z0, s0), Inf) < 1e-8
@test norm(cone_prod(s0, z0) - cone_mat(s0)*z0, Inf) < 1e-8
Δs0 = [0, 0, 4.]
@test norm(soc_step_length(s0, Δs0, τ = 1.0) - 0.5, Inf) < 1e-8
Δs0 = [0, 0, 1.]
@test norm(soc_step_length(s0, Δs0, τ = 1.0) - 1.0, Inf) < 1e-8
z1 = [3, 2, 1.]
Δz1 = [0, 0, 1.]
@test norm(soc_step_length(z1, Δz1, τ = 1.0) - 1.0, Inf) < 1e-8

nd_scaling(s0, z0)
nd_λ(s0, z0)

Δ = rand(nx+nb+2ns)
@test norm(Δ - vcat(unpack(Δ)...), Inf) < 1e-8
