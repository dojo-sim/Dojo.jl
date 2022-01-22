Rotations.UnitQuaternion(w::T, v::StaticVector{3,T}, normalize::Bool = true) where T = UnitQuaternion{T}(w, v[1], v[2], v[3], normalize)
Rotations.UnitQuaternion(w::T, v::Vector{T}, normalize::Bool = true) where T = (@assert length(v)==3; UnitQuaternion{T}(w, v[1], v[2], v[3], normalize))
Rotations.UnitQuaternion(v::StaticVector{3,T}) where T = pure_quaternion(v)
Rotations.UnitQuaternion(v::Vector) = (@assert length(v)==3; pure_quaternion(v))

@inline imag(q::UnitQuaternion) = Rotations.vector(q)

qrotate(q1::UnitQuaternion,q2::UnitQuaternion) = q2 * q1 / q2
vrotate(v::Vector,q::UnitQuaternion) = imag(qrotate(pure_quaternion(v), q))
vrotate(v::StaticVector,q::UnitQuaternion) = q*v

@inline rotation_matrix(q::UnitQuaternion) = VRᵀmat(q) * LVᵀmat(q)
@inline ∂vrotate∂p(p::AbstractVector, q::UnitQuaternion) = VRᵀmat(q) * LVᵀmat(q)
@inline ∂vrotate∂q(p::AbstractVector, q::UnitQuaternion) = VLmat(q) * Lmat(UnitQuaternion(p)) * Tmat() + VRᵀmat(q) * Rmat(UnitQuaternion(p))

rotation_vector(q::UnitQuaternion) = rotation_angle(q) * rotation_axis(q)
@inline function axisangle2quaternion(v)
    θ = norm(v)
    q = UnitQuaternion(cos(θ/2), 1/2 * sinc(θ/(2π)) * v)
end

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

vector(q::UnitQuaternion) = SA[q.w, q.x, q.y, q.z]

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



################################################################################
# Matrix-Vector Product Jacobian
################################################################################
function ∂qVLmat(p) # 𝞉(VLmat(q)*p)/∂q
	SA[
    	0     p[1]  p[2]  p[3];
    	p[1]  0     p[3] -p[2];
    	p[2] -p[3]  0     p[1];
    	p[3]  p[2] -p[1]  0;
    ]
end

function ∂qLᵀVᵀmat(p) # 𝞉(VLmat(q)*p)/∂q
	SA[
    	0     p[1]  p[2]  p[3];
    	p[1]  0    -p[3]  p[2];
    	p[2]  p[3]  0    -p[1];
    	p[3] -p[2]  p[1]  0;
    ]
end

function ∂qVRmat(p) # 𝞉(VRmat(q)*p)/∂q
	SA[
		p[2]  p[1] -p[4]  p[3];
		p[3]  p[4]  p[1] -p[2];
		p[4] -p[3]  p[2]  p[1];
    ]
end

function ∂qRᵀVᵀmat(p) # 𝞉(RᵀVᵀmat(q)*p)/∂q
	SA[
    	p[2]  p[1]  p[4] -p[3];
    	p[3] -p[4]  p[1]  p[2];
    	p[4]  p[3] -p[2]  p[1];
    ]
end

function ∂pskew(λ) # 𝞉(skew(p)*λ)/∂p
	skew(-λ)
end
