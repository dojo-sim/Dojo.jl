function cone_product(u::AbstractVector{T}, v::AbstractVector{T}) where {T}
    [u'*v; u[1] * v[2:end] + v[1] * u[2:end]]
end

function cone_product(u::SVector{N,T}, v::SVector{N,T}) where {N,T}
    vcat(u'*v, u[1] * v[SVector{N-1}(2:end)] + v[1] * u[SVector{N-1}(2:end)])
end

@inline function ∇cone_product(u::SVector{3,T}) where {T}
    SMatrix{3,3,T,9}(u[1], u[2], u[3], u[2], u[1], 0, u[3], 0, u[1])
end
