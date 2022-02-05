Rotations.UnitQuaternion(w::T, v::StaticVector{3,T}, normalize::Bool = true) where T = UnitQuaternion{T}(w, v[1], v[2], v[3], normalize)
Rotations.UnitQuaternion(w::T, v::Vector{T}, normalize::Bool = true) where T = (@assert length(v)==3; UnitQuaternion{T}(w, v[1], v[2], v[3], normalize))
Rotations.UnitQuaternion(v::StaticVector{3,T}) where T = pure_quaternion(v)
Rotations.UnitQuaternion(v::Vector) = (@assert length(v)==3; pure_quaternion(v))

@inline imag(q::UnitQuaternion) = Rotations.vector(q)

qrotate(q1::UnitQuaternion,q2::UnitQuaternion) = q2 * q1 / q2
vrotate(v::Vector,q::UnitQuaternion) = imag(qrotate(pure_quaternion(v), q))
vrotate(v::StaticVector,q::UnitQuaternion) = q*v

@inline rotation_matrix(q::UnitQuaternion) = VRᵀmat(q) * LVᵀmat(q)
# ∂(rotation_matrix(q)*p)/∂q
@inline ∂qrotation_matrix(q::UnitQuaternion, p::AbstractVector) =
 	∂qVRᵀmat(LVᵀmat(q) * p) + VRᵀmat(q) * ∂qLVᵀmat(p)
# ∂(rotation_matrix(inv(q))*p)/∂q
@inline ∂qrotation_matrix_inv(q::UnitQuaternion, p::AbstractVector) =
 	∂qrotation_matrix(inv(q), p) * Tmat()

@inline ∂vrotate∂p(p::AbstractVector, q::UnitQuaternion) = VRᵀmat(q) * LVᵀmat(q)
@inline ∂vrotate∂q(p::AbstractVector, q::UnitQuaternion) = VLmat(q) * Lmat(UnitQuaternion(p)) * Tmat() + VRᵀmat(q) * Rmat(UnitQuaternion(p))

rotation_vector(q::UnitQuaternion) = rotation_angle(q) * rotation_axis(q)
# function rotation_vector1(q::UnitQuaternion{T}) where T
#     angle = 2 * acos(q.w)
# 	ϵ = max(1-q.w*q.w, 1e-10)
#     x = q.x / sqrt(ϵ)
#     y = q.y / sqrt(ϵ)
#     z = q.z / sqrt(ϵ)
#     return angle * SVector{3,T}(x, y, z)
# end
function ∂qrotation_vector(q::UnitQuaternion{T}) where T
    # v = amp * Vmat(q)
	ϵ = max(1-q.w*q.w, 1e-10)
    amp = 2 * acos(q.w) / sqrt(ϵ)
    ∇amp = SVector{4,T}(2q.w*acos(q.w)*(sqrt(ϵ)^-3) - 2(sqrt(ϵ)^-2), 0, 0, 0)
    ∇ = Vmat(q) * ∇amp' + amp * Vmat()
    return ∇
end

@inline function axis_angle_to_quaternion(v)
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
function ∂qVLmat(p::AbstractVector) # 𝞉(VLmat(q)*p)/∂q
	SA[
    	0     p[1]  p[2]  p[3];
    	p[1]  0     p[3] -p[2];
    	p[2] -p[3]  0     p[1];
    	p[3]  p[2] -p[1]  0;
    ]
end

function ∂qLVᵀmat(p::AbstractVector) # 𝞉(∂qLVᵀmat(q)*p)/∂q
	SA[
    	0    -p[1] -p[2] -p[3];
    	p[1]  0     p[3] -p[2];
    	p[2] -p[3]  0     p[1];
    	p[3]  p[2] -p[1]  0;
    ]
end

function ∂qVLᵀmat(p::AbstractVector) # 𝞉(VLᵀmat(q)*p)/∂q
	SA[
		p[2] -p[1] -p[4]  p[3];
		p[3]  p[4] -p[1] -p[2];
		p[4] -p[3]  p[2] -p[1];
    ]
end

function ∂qLᵀVᵀmat(p::AbstractVector) # 𝞉(LᵀVᵀmat(q)*p)/∂q
	SA[
    	0     p[1]  p[2]  p[3];
    	p[1]  0    -p[3]  p[2];
    	p[2]  p[3]  0    -p[1];
    	p[3] -p[2]  p[1]  0;
    ]
end

function ∂qVRmat(p::AbstractVector) # 𝞉(VRmat(q)*p)/∂q
	SA[
		p[2]  p[1] -p[4]  p[3];
		p[3]  p[4]  p[1] -p[2];
		p[4] -p[3]  p[2]  p[1];
    ]
end

function ∂qRᵀVᵀmat(p::AbstractVector) # 𝞉(RᵀVᵀmat(q)*p)/∂q
	SA[
    	p[2]  p[1]  p[4] -p[3];
    	p[3] -p[4]  p[1]  p[2];
    	p[4]  p[3] -p[2]  p[1];
    ]
end

function ∂qVRᵀmat(p::AbstractVector) # 𝞉(RᵀVᵀmat(q)*p)/∂q
	SA[
    	p[2] -p[1]  p[4] -p[3];
    	p[3] -p[4] -p[1]  p[2];
    	p[4]  p[3] -p[2] -p[1];
    ]
end

function ∂qRᵀmat(p::AbstractVector) # 𝞉(Rᵀmat(q)*p)/∂q
	SA[
    	p[1]  p[2]  p[3]  p[4];
    	p[2] -p[1]  p[4] -p[3];
    	p[3] -p[4] -p[1]  p[2];
    	p[4]  p[3] -p[2] -p[1];
    ]
end

function ∂qLmat(p::AbstractVector) # 𝞉(Lmat(q)*p)/∂q
	SA[
    	p[1] -p[2] -p[3] -p[4];
    	p[2]  p[1]  p[4] -p[3];
    	p[3] -p[4]  p[1]  p[2];
    	p[4]  p[3] -p[2]  p[1];
    ]
end

function ∂pskew(λ) # 𝞉(skew(p)*λ)/∂p
	skew(-λ)
end


# using Symbolics
# @variables q_[1:4], p3_[1:3], p4_[1:4]
# qq_ = UnitQuaternion(q_, false)
# Symbolics.jacobian(Rᵀmat(qq_) * p4_, q_)
# Symbolics.jacobian(Lmat(qq_) * p4_, q_)
