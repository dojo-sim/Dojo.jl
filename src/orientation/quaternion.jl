Rotations.UnitQuaternion(w::T, v::StaticVector{3,T}, normalize::Bool = true) where T = UnitQuaternion{T}(w, v[1], v[2], v[3], normalize)
Rotations.UnitQuaternion(w::T, v::Vector{T}, normalize::Bool = true) where T = (@assert length(v)==3; UnitQuaternion{T}(w, v[1], v[2], v[3], normalize))
Rotations.UnitQuaternion(v::StaticVector{3,T}) where T = pure_quaternion(v)
# Rotations.UnitQuaternion(v::Vector) = (@assert length(v)==3; pure_quaternion(v))

imag(q::UnitQuaternion) = Rotations.vector(q)

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
vector(q::AbstractVector) = q 

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

################################################################################
# Matrix-Vector Product Jacobian
################################################################################
function ∂VLmat∂q(p::AbstractVector) # 𝞉(VLmat(q)*p)/∂q
	SA[
    	0     p[1]  p[2]  p[3];
    	p[1]  0     p[3] -p[2];
    	p[2] -p[3]  0     p[1];
    	p[3]  p[2] -p[1]  0;
    ]
end

function ∂LVᵀmat∂q(p::AbstractVector) # 𝞉(∂LVᵀmat∂q(q)*p)/∂q
	SA[
    	0    -p[1] -p[2] -p[3];
    	p[1]  0     p[3] -p[2];
    	p[2] -p[3]  0     p[1];
    	p[3]  p[2] -p[1]  0;
    ]
end

function ∂VLᵀmat∂q(p::AbstractVector) # 𝞉(VLᵀmat(q)*p)/∂q
	SA[
		p[2] -p[1] -p[4]  p[3];
		p[3]  p[4] -p[1] -p[2];
		p[4] -p[3]  p[2] -p[1];
    ]
end

function ∂LᵀVᵀmat∂q(p::AbstractVector) # 𝞉(LᵀVᵀmat(q)*p)/∂q
	SA[
    	0     p[1]  p[2]  p[3];
    	p[1]  0    -p[3]  p[2];
    	p[2]  p[3]  0    -p[1];
    	p[3] -p[2]  p[1]  0;
    ]
end

function ∂VRmat∂q(p::AbstractVector) # 𝞉(VRmat(q)*p)/∂q
	SA[
		p[2]  p[1] -p[4]  p[3];
		p[3]  p[4]  p[1] -p[2];
		p[4] -p[3]  p[2]  p[1];
    ]
end

function ∂RᵀVᵀmat∂q(p::AbstractVector) # 𝞉(RᵀVᵀmat(q)*p)/∂q
	SA[
    	p[2]  p[1]  p[4] -p[3];
    	p[3] -p[4]  p[1]  p[2];
    	p[4]  p[3] -p[2]  p[1];
    ]
end

function ∂VRᵀmat∂q(p::AbstractVector) # 𝞉(RᵀVᵀmat(q)*p)/∂q
	SA[
    	p[2] -p[1]  p[4] -p[3];
    	p[3] -p[4] -p[1]  p[2];
    	p[4]  p[3] -p[2] -p[1];
    ]
end

function ∂Rᵀmat∂q(p::AbstractVector) # 𝞉(Rᵀmat(q)*p)/∂q
	SA[
    	p[1]  p[2]  p[3]  p[4];
    	p[2] -p[1]  p[4] -p[3];
    	p[3] -p[4] -p[1]  p[2];
    	p[4]  p[3] -p[2] -p[1];
    ]
end

function ∂Lmat∂q(p::AbstractVector) # 𝞉(Lmat(q)*p)/∂q
	SA[
    	p[1] -p[2] -p[3] -p[4];
    	p[2]  p[1]  p[4] -p[3];
    	p[3] -p[4]  p[1]  p[2];
    	p[4]  p[3] -p[2]  p[1];
    ]
end

function ∂skew∂p(λ) # 𝞉(skew(p)*λ)/∂p
	skew(-λ)
end