Rotations.UnitQuaternion(w::T, v::StaticVector{3,T}, normalize::Bool = true) where T = UnitQuaternion{T}(w, v[1], v[2], v[3], normalize)
Rotations.UnitQuaternion(w::T, v::Vector{T}, normalize::Bool = true) where T = (@assert length(v)==3; UnitQuaternion{T}(w, v[1], v[2], v[3], normalize))
Rotations.UnitQuaternion(v::StaticVector{3,T}) where T = pure_quaternion(v)
Rotations.UnitQuaternion(v::Vector) = (@assert length(v)==3; pure_quaternion(v))

@inline imag(q::UnitQuaternion) = Rotations.vector(q)

qrotate(q1::UnitQuaternion,q2::UnitQuaternion) = q2 * q1 / q2
vrotate(v::Vector,q::UnitQuaternion) = imag(qrotate(pure_quaternion(v), q))
vrotate(v::StaticVector,q::UnitQuaternion) = q*v

@inline ∂vrotate∂p(p::AbstractVector, q::UnitQuaternion) = error("not implemented")
@inline ∂vrotate∂q(p::AbstractVector, q::UnitQuaternion) = VLmat(q) * Lmat(UnitQuaternion(p)) * Tmat() + VRᵀmat(q) * Rmat(UnitQuaternion(p))

rotation_vector(q::UnitQuaternion) = rotation_angle(q) * rotation_axis(q)

Lmat(q) = lmult(q)
Lᵀmat(q) = lmult(q)'
Rmat(q) = rmult(q)
Rᵀmat(q) = rmult(q)'

# Remove once added to Rotations.jl
function Base.:/(q::UnitQuaternion, w::Real)
    return UnitQuaternion(q.w/w, q.x/w, q.y/w, q.z/w, false)
end

Tmat(::Type{T}=Float64) where T = tmat(T)
Tᵀmat(::Type{T}=Float64) where T = tmat(T)
Vmat(::Type{T}=Float64) where T = vmat(T)
Vᵀmat(::Type{T}=Float64) where T = hmat(T)
Vmat(q::UnitQuaternion) = imag(q)

function VLmat(q::UnitQuaternion)
    SA[
        q.x  q.w -q.z  q.y;
        q.y  q.z  q.w -q.x;
        q.z -q.y  q.x  q.w;
    ]
end
function VLᵀmat(q::UnitQuaternion)
    SA[
        -q.x  q.w  q.z -q.y;
        -q.y -q.z  q.w  q.x;
        -q.z  q.y -q.x  q.w;
    ]
end
function VRmat(q::UnitQuaternion)
    SA[
        q.x  q.w  q.z -q.y;
        q.y -q.z  q.w  q.x;
        q.z  q.y -q.x  q.w;
    ]
end
function VRᵀmat(q::UnitQuaternion)
    SA[
        -q.x  q.w -q.z  q.y;
        -q.y  q.z  q.w -q.x;
        -q.z -q.y  q.x  q.w;
    ]
end

function LVᵀmat(q::UnitQuaternion)
    SA[
        -q.x -q.y -q.z;
         q.w -q.z  q.y;
         q.z  q.w -q.x;
        -q.y  q.x  q.w;
    ]
end
function LᵀVᵀmat(q::UnitQuaternion)
    SA[
         q.x  q.y  q.z;
         q.w  q.z -q.y;
        -q.z  q.w  q.x;
         q.y -q.x  q.w;
    ]
end
function RVᵀmat(q::UnitQuaternion)
    SA[
        -q.x -q.y -q.z;
         q.w  q.z -q.y;
        -q.z  q.w  q.x;
         q.y -q.x  q.w;
    ]
end
function RᵀVᵀmat(q::UnitQuaternion)
    SA[
         q.x  q.y  q.z;
         q.w -q.z  q.y;
         q.z  q.w -q.x;
        -q.y  q.x  q.w;
    ]
end

function ∂L∂qsplit(::Type{T}) where T
    return SA{T}[
        1 0 0 0
        0 1 0 0
        0 0 1 0
        0 0 0 1

        0 -1 0 0
        1 0 0 0
        0 0 0 1
        0 0 -1 0

        0 0 -1 0
        0 0 0 -1
        1 0 0 0
        0 1 0 0

        0 0 0 -1
        0 0 1 0
        0 -1 0 0
        1 0 0 0
    ]
end
function ∂Lᵀ∂qsplit(::Type{T}) where T
    return SA{T}[
        1 0 0 0
        0 -1 0 0
        0 0 -1 0
        0 0 0 -1

        0 1 0 0
        1 0 0 0
        0 0 0 -1
        0 0 1 0

        0 0 1 0
        0 0 0 1
        1 0 0 0
        0 -1 0 0

        0 0 0 1
        0 0 -1 0
        0 1 0 0
        1 0 0 0
    ]
end
function ∂R∂qsplit(::Type{T}) where T
    return SA{T}[
        1 0 0 0
        0 1 0 0
        0 0 1 0
        0 0 0 1

        0 -1 0 0
        1 0 0 0
        0 0 0 -1
        0 0 1 0

        0 0 -1 0
        0 0 0 1
        1 0 0 0
        0 -1 0 0

        0 0 0 -1
        0 0 -1 0
        0 1 0 0
        1 0 0 0
    ]
end
function ∂Rᵀ∂qsplit(::Type{T}) where T
    return SA{T}[
        1 0 0 0
        0 -1 0 0
        0 0 -1 0
        0 0 0 -1

        0 1 0 0
        1 0 0 0
        0 0 0 1
        0 0 -1 0

        0 0 1 0
        0 0 0 -1
        1 0 0 0
        0 1 0 0

        0 0 0 1
        0 0 1 0
        0 -1 0 0
        1 0 0 0
    ]
end

function slerp(q1,q2,h)
    s = params(q1)'*params(q2)
    if s < 0
        s = -s
        q2 = -q2
    end

    qdiff = q1\q2
    φdiff = rotation_angle(qdiff)
    udiff = rotation_axis(qdiff)
    φint = φdiff*h
    qint = UnitQuaternion(cos(φint/2),udiff*sin(φint/2),false)

    return q1*qint
end
