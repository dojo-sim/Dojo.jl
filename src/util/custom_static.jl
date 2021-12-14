function Base.convert(T::Type{Vector{SVector{N2,T1}}}, M::Matrix) where {N2,T1}
    N1 = size(M)[1]
    @assert size(M)[2] == N2
    Mout = [szeros(T1, N2) for i = 1:N1]
    for i = 1:N1
        Mout[i] = convert(SVector{N2,T1}, M[i,:])
    end
    return Mout
end


Base.:*(u::LinearAlgebra.AdjointAbsVec, v::SVector{1,T}) where T = u * v[1]
Base.:*(u::AbstractVector, v::SVector{1,T}) where T = u * v[1]


Base.hcat(E::UniformScaling, v::SVector{N,T}) where {T, N} = hcat(SMatrix{N,N,T,N*N}(E), v)
Base.hcat(v::SVector{N,T}, E::UniformScaling) where {T, N} = hcat(v, SMatrix{N,N,T,N*N}(E))
Base.hcat(E::UniformScaling, v::Adjoint{T,SVector{N,T}}) where {T, N} = @error "not implemented" # shcat(E.λ, v)
Base.hcat(v::Adjoint{T,SVector{N,T}}, E::UniformScaling) where {T, N} = @error "not implemented" # shcat(v, E.λ)
Base.hcat(E::UniformScaling, A::SMatrix{N1,N2,T,N1N2}) where {T, N1, N2, N1N2} = hcat(SMatrix{N1,N1,T,N1*N1}(E), A)
Base.hcat(A::SMatrix{N1,N2,T,N1N2}, E::UniformScaling) where {T, N1, N2, N1N2} = hcat(A, SMatrix{N1,N1,T,N1*N1}(E))

Base.vcat(E::UniformScaling, v::SVector{N,T}) where {T, N} = svcat(E.λ, v)
Base.vcat(v::SVector{N,T}, E::UniformScaling) where {T, N} = svcat(v, E.λ)
Base.vcat(E::UniformScaling, v::Adjoint{T,SVector{N,T}}) where {T, N} = vcat(SMatrix{N,N,T,N*N}(E), v)
Base.vcat(v::Adjoint{T,SVector{N,T}}, E::UniformScaling) where {T, N} = vcat(v, SMatrix{N,N,T,N*N}(E))
Base.vcat(E::UniformScaling, A::SMatrix{N1,N2,T,N1N2}) where {T, N1, N2, N1N2} = vcat(SMatrix{N2,N2,T,N2*N2}(E), A)
Base.vcat(A::SMatrix{N1,N2,T,N1N2}, E::UniformScaling) where {T, N1, N2, N1N2} = vcat(A, SMatrix{N2,N2,T,N2*N2}(E))

Base.convert(::Type{SMatrix{N,N,T,N2}}, E::UniformScaling) where {T, N, N2} = SMatrix{N,N,T,N2}(E)

@inline svcat(a::T) where T = SA[a]
@inline svcat(a::StaticArray) = a
@inline svcat(a::T, b::T) where T = SA[a; b]
@inline svcat(a::StaticArray, b::StaticArray) = vcat(a, b)
@inline svcat(a::StaticArray{Tuple{N},T,1}, b::T) where {T,N} = vcat(a, SA[b])
@inline svcat(a::T, b::StaticArray{Tuple{N},T,1}) where {T,N} = vcat(SA[a], b)
@inline svcat(a, b, c...) = svcat(svcat(a, b), svcat(c...))

@inline szeros(::Type{T}, N) where T = @SVector zeros(T, N)
@inline szeros(N)= @SVector zeros(N)
@inline szeros(::Type{T}, N1, N2) where T = @SMatrix zeros(T, N1, N2)
@inline szeros(N1, N2)= @SMatrix zeros(N1, N2)

@inline sones(::Type{T}, N) where T = @SVector ones(T, N)
@inline sones(N)= @SVector ones(N)
@inline sones(::Type{T}, N1, N2) where T = @SMatrix ones(T, N1, N2)
@inline sones(N1, N2)= @SMatrix ones(N1, N2)

@inline srand(::Type{T}, N) where T = @SVector rand(T, N)
@inline srand(N)= @SVector rand(N)
@inline srand(::Type{T}, N1, N2) where T = @SMatrix rand(T, N1, N2)
@inline srand(N1, N2)= @SMatrix rand(N1, N2)


sisnan(a::StaticArray) = any(isnan.(a))


# To fix StaticArray bug
zerodimstaticadjoint(A) = A'
zerodimstaticadjoint(::SMatrix{0,N,T,0}) where {T,N} = SMatrix{N,0,T,0}()
zerodimstaticadjoint(::SMatrix{N,0,T,0}) where {T,N} = SMatrix{0,N,T,0}()
Base.:*(x::LinearAlgebra.AdjointAbsVec, A::SMatrix{0,N,T,0}) where {T,N} = (zerodimstaticadjoint(A)*x')'
