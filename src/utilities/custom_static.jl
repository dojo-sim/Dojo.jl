Base.:*(u::AbstractVector, v::SVector{1,T}) where T = u * v[1]

Base.vcat(v::Adjoint{T,SVector{N,T}}, E::UniformScaling) where {T, N} = vcat(v, SMatrix{N,N,T,N*N}(E))
Base.vcat(A::SMatrix{N1,N2,T,N1N2}, E::UniformScaling) where {T, N1, N2, N1N2} = vcat(A, SMatrix{N2,N2,T,N2*N2}(E))

svcat(a::T, b::T) where T = SA[a; b]
svcat(a::StaticArray, b::StaticArray) = vcat(a, b)
svcat(a::StaticArray{Tuple{N},T,1}, b::T) where {T,N} = vcat(a, SA[b])
svcat(a::T, b::StaticArray{Tuple{N},T,1}) where {T,N} = vcat(SA[a], b)

szeros(::Type{T}, N) where T = @SVector zeros(T, N)
szeros(N)= @SVector zeros(N)
szeros(::Type{T}, N1, N2) where T = @SMatrix zeros(T, N1, N2)
szeros(N1, N2)= @SMatrix zeros(N1, N2)

sones(::Type{T}, N) where T = @SVector ones(T, N)
sones(N)= @SVector ones(N)

srand(N)= @SVector rand(N)

# TODO: check StaticArray bug fix, then remove
zerodimstaticadjoint(A) = A'
zerodimstaticadjoint(::SMatrix{0,N,T,0}) where {T,N} = SMatrix{N,0,T,0}()
