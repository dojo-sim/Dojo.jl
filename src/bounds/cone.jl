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
    @assert norm(W̄i - inv(W̄), Inf) < 1e-8
    @assert norm(Wi - inv(W), Inf) < 1e-8
    @assert norm(λ - inv(W)*s) < 1e-8
    return W, Wi, λ
end

function ort_nt_scaling(s::AbstractVector{T}, z::AbstractVector{T}) where {T}
    # http://www.seas.ucla.edu/~vandenbe/publications/coneprog.pdf
    W = Diagonal(sqrt.(s) ./ sqrt.(z))
    Wi = Diagonal(sqrt.(z) ./ sqrt.(s))
    λ = W*z
    return W, Wi, λ
end



function ∇cone_product(u::AbstractVector{T}) where {T}
    n = length(u)
    U = zeros(n,n)
    U += u[1] * I(n)
    U[1,2:end] += u[2:end]
    U[2:end, 1] += u[2:end]
    return U
end

@inline function ∇cone_product(u::SVector{3,T}) where {T}
    SMatrix{3,3,T,9}(u[1], u[2], u[3], u[2], u[1], 0, u[3], 0, u[1])
end

function cone_product(u::AbstractVector{T}, v::AbstractVector{T}) where {T}
    [u'*v; u[1] * v[2:end] + v[1] * u[2:end]]
end

function cone_product(u::SVector{N,T}, v::SVector{N,T}) where {N,T}
    vcat(u'*v, u[1] * v[SVector{N-1}(2:end)] + v[1] * u[SVector{N-1}(2:end)])
end
