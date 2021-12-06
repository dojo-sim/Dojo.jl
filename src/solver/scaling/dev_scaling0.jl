using LinearAlgebra

################################################################################
# Quadratic Cone Program
################################################################################

nx = 6
ns = 3
nb = 2

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
# f(x) = 0.5 * x'*P*x + c'*x
# ce(x) = A*x - b
# ci(x, s) = D*x + s - h
# val_soc(s) = s[1]^2 - s[2:end]'*s[2:end]
# ϕ(s) = - 0.5 * log(s[1]^2 - s[2:end]'*s[2:end])
jordan(n) = Diagonal([1; -ones(n-1)])

# function g(s)
#     n = length(s)
#     grad = - 1/(s'*J*s) * J * s
#     return grad
# end
# function H(s)
#     n = length(s)
#     J = - I(n)
#     J[1,1] = 1
#     hess = 1 / (s'*J*s)^2 * (2*J*s*s'*J - (s'*J*s)*J)
#     return hess
# end

cone_prod(s, z) = [s'*z; s[1]*z[2:end] + z[1]*s[2:end]]
cone_neutral(n) = [1; zeros(n-1)]
# function cone_inverse(s)
#     n = length(s)
#     J = jordan(n)
#     si = 1 / (s'*J*s) * J *s
#     return si
# end
# function cone_square(s)
#     ss = [s'*s; 2s[1]*s[2:end]]
#     return ss
# end
# function cone_sqrt(s)
#     n = length(s)
#     J = jordan(n)
#     t = sqrt(s'*J*s)
#     sq = 1 / sqrt(2 * (s[1] + t)) * [s[1] + t; s[2:end]]
#     return sq
# end

# function w_scaling(s, z)
#     z_ = cone_inverse(cone_sqrt(z))
#     w = H(cone_sqrt(z)) * sqrt(H(z_)*s)
#     return w
# end
#
# function W_scaling(s, z)
#     w = w_scaling(s, z)
#     W = H(cone_inverse(cone_sqrt(w)))
#     return W
# end

function nd_point(s, z)
    n = length(s)
    J = jordan(n)
    z̄ = 1 / sqrt(z'*J*z) * z
    s̄ = 1 / sqrt(s'*J*s) * s
    γ = sqrt((1 + z̄'*s̄) / 2)
    ωb = 1 / (2γ) * (s̄ + J*z̄)
    @assert norm(ωb'*J*ωb - 1) < 1e-8
    @assert norm(ωb'*z̄ - γ) < 1e-8
    @assert norm(ωb'*J*s̄ - γ) < 1e-8
    r = sqrt((s'*J*s)/(z'*J*z))
    ω = sqrt(r) * ωb
    return ω
end

function nd_scaling(s, z)
    n = length(s)
    J = jordan(n)
    e = cone_neutral(n)
    z̄ = 1 / sqrt(z'*J*z) * z
    s̄ = 1 / sqrt(s'*J*s) * s
    γ = sqrt((1 + z̄'*s̄) / 2)
    ωb = 1 / (2γ) * (s̄ + J*z̄)

    v = 1 / sqrt(2 * (ωb[1] + 1)) * (ωb + e)
    W̄ = 2 * v * v' - J
    W̄i = 2 * J * v * v' * J - J
    @assert norm(W̄ * (2J*ωb*ωb'*J - J) * W̄ - I(n)) < 1e-8
    r = sqrt((s'*J*s)/(z'*J*z))
    W = sqrt(r) * W̄
    return W
end

function nd_λ(s, z)
    n = length(s)
    J = jordan(n)
    e = cone_neutral(n)
    z̄ = 1 / sqrt(z'*J*z) * z
    s̄ = 1 / sqrt(s'*J*s) * s
    γ = sqrt((1 + z̄'*s̄) / 2)
    ωb = 1 / (2γ) * (s̄ + J*z̄)

    v = 1 / sqrt(2 * (ωb[1] + 1)) * (ωb + e)
    W̄ = 2 * v * v' - J
    W̄i = 2 * J * v * v' * J - J

    λb = W̄ * z̄
    @assert norm(λb - W̄i * s̄) < 1e-8
    @assert norm(λb - J*W̄*J*s̄) < 1e-8
    @assert norm(λb'*J*λb - 1) < 1e-8
    λ = W * z
    @assert norm(λ - inv(nd_scaling(s,z))*s) < 1e-8


    return λ
end


s = [1, 0.5, 0]
z = [1, 0.5, 0]
ω = nd_point(s, z)
W = nd_scaling(s, z)
λ = nd_λ(s, z)



function solve(x, y, z, s; iter::Int = 10, btol = 1e-3, rtol = 1e-3)
    m = 1
    n = length(s)
    e = cone_neutral(n)

    for i = 1:iter
        # Residual
        res = residual(x, y, z, s)
        μ = s'*z / m

        # Termination
        cvio_s = min(val_soc(s), 0)
        cvio_v = min(val_soc(z), 0)
        cvio = max(cvio_s, cvio_z)
        bvio = norm(z'*s, Inf)
        rvio = norm(res, Inf)
        if bvio < bol && rvio < rtol
            return x, y, z, s
        end

        # Affine direction
        jac = jacobian(s, z)
        λ = nd_λ(s, z)
        Δa = - jac \ [res; cone_prod(λ, λ)]
        Δxa, Δya, Δza, Δsa = unpack(Δa)

        # Step size and Centering
        αs = soc_step_length(s, Δsa, τ = 1.0)
        αz = soc_step_length(z, Δza, τ = 1.0)
        α = min(αs, αz)
        ρ = (s + α*Δsa)' * (z + α*Δza) / (s'*z)
        σ = clamp(ρ, 0, 1)^3

        # Correction
        λ = nd_λ(s, z)
        Δ = - jac \ [(1-σ)*res; cone_prod(λ, λ) + γ*cone_prod(inv(W)'*Δsa, W*Δza) - σ*μ*e]
        Δx, Δy, Δz, Δs = unpack(Δ)

        # Update iterates
        αs = soc_step_length(λ, inv(W)'*Δs, τ = 0.99)
        αz = soc_step_length(λ, W*Δz, τ = 0.99)
        α = min(αs, αz)
        x = x + α*Δx
        y = y + α*Δy
        z = z + α*Δz
        s = s + α*Δs
    end

    return x, y, z, s
end

x0 = rand(nx)
y0 = rand(nb)
z0 = rand(ns)
s0 = rand(ns)
solve(x0, y0, z0, s0, iter = 1)

function unpack(Δ)
    Δx = Δ[1:nx]
    Δy = Δ[nx+1:nx+nb]
    Δz = Δ[nx+nb+1:nx+nb+ns]
    Δs = Δ[nx+nb+ns+1:nx+nb+2ns]
    return Δx, Δy, Δz, Δs
end

function residual(x, y, z, s)
    jac = jacobian(s, z)
    res = jac[1:nx+nb+ns, 1:nx+nb+ns] * [x; y; z] + [c; -b; -h] + [zeros(nx+nb); s]
    return res
end

function jacobian(s, z)
    W = nd_scaling(s, z)
    λ = nd_λ(s, z)
    mλ = mat_prod(λ)
    jac = [P A' G';
         A zeros(nb,nb) zeros(nb,ns) zeros(nb,ns);
         D zeros(ns,nb) zeros(ns,ns) I(ns);
         zeros(ns,nx) zeros(ns,ns) mλ*W mλ*inv(W)']
    return jac
end

function mat_prod(s)
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
    ρv = Δ[2:end] / sqrt(λ_λ)dev_scaling
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


s = [2; rand(2)]
z = [3; rand(2)]
val_soc(s) > 0
val_soc(z) > 0
c = cone_prod(s, z)
c - mat_prod(s) * z
