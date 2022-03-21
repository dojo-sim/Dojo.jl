MeshCat.LinearMap(q::Quaternion) = MeshCat.LinearMap(rotation_matrix(q))
MeshCat.js_quaternion(q::Quaternion) = [q.v1, q.v2, q.v3, q.s]

RotX(θ) = Quaternion(cos(θ/2), sin(θ/2), 0, 0)
RotY(θ) = Quaternion(cos(θ/2), 0, sin(θ/2), 0)
RotZ(θ) = Quaternion(cos(θ/2), 0, 0, sin(θ/2))

quateltype(x) = eltype(x) # TODO not super elegant
quateltype(::Quaternion{T}) where T = T

vector(q::Quaternion) = SA[q.s, q.v1, q.v2, q.v3]
vector(q::AbstractVector) = q

function Lmat(q::Quaternion)
    SA[
        q.s  -q.v1 -q.v2 -q.v3;
        q.v1  q.s  -q.v3  q.v2;
        q.v2  q.v3  q.s  -q.v1;
        q.v3 -q.v2  q.v1  q.s;
    ]
end

function Rmat(q::Quaternion)
    SA[
        q.s  -q.v1 -q.v2 -q.v3;
        q.v1  q.s   q.v3 -q.v2;
        q.v2 -q.v3  q.s   q.v1;
        q.v3  q.v2 -q.v1  q.s;
    ]
end

Lᵀmat(q) = Lmat(q)'
Rᵀmat(q) = Rmat(q)'

function Tmat(::Type{T}=Float64) where T
    SA{T}[
        1  0  0  0;
        0 -1  0  0;
        0  0 -1  0;
        0  0  0 -1;
    ]
end

function Vmat(::Type{T}=Float64) where T
    SA{T}[
        0 1 0 0;
        0 0 1 0;
        0 0 0 1;
    ]
end

function Vᵀmat(::Type{T}=Float64) where T
    SA{T}[
        0 0 0;
        1 0 0;
        0 1 0;
        0 0 1;
    ]
end

Vmat(q::Quaternion) = SA[q.v1, q.v2, q.v3]

function VLmat(q::Quaternion)
    SA[
        q.v1  q.s -q.v3  q.v2;
        q.v2  q.v3  q.s -q.v1;
        q.v3 -q.v2  q.v1  q.s;
    ]
end

function VLᵀmat(q::Quaternion)
    SA[
        -q.v1  q.s   q.v3 -q.v2;
        -q.v2 -q.v3  q.s   q.v1;
        -q.v3  q.v2 -q.v1  q.s;
    ]
end

function VRmat(q::Quaternion)
    SA[
        q.v1  q.s  q.v3 -q.v2;
        q.v2 -q.v3  q.s  q.v1;
        q.v3  q.v2 -q.v1  q.s;
    ]
end

function VRᵀmat(q::Quaternion)
    SA[
        -q.v1  q.s -q.v3  q.v2;
        -q.v2  q.v3  q.s -q.v1;
        -q.v3 -q.v2  q.v1  q.s;
    ]
end

function LVᵀmat(q::Quaternion)
    SA[
        -q.v1 -q.v2 -q.v3;
         q.s -q.v3  q.v2;
         q.v3  q.s -q.v1;
        -q.v2  q.v1  q.s;
    ]
end

function LᵀVᵀmat(q::Quaternion)
    SA[
         q.v1  q.v2  q.v3;
         q.s  q.v3 -q.v2;
        -q.v3  q.s  q.v1;
         q.v2 -q.v1  q.s;
    ]
end

function RVᵀmat(q::Quaternion)
    SA[
        -q.v1 -q.v2 -q.v3;
         q.s  q.v3 -q.v2;
        -q.v3  q.s  q.v1;
         q.v2 -q.v1  q.s;
    ]
end

function RᵀVᵀmat(q::Quaternion)
    SA[
         q.v1  q.v2  q.v3;
         q.s -q.v3  q.v2;
         q.v3  q.s -q.v1;
        -q.v2  q.v1  q.s;
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

function skew(p)
    SA[
    	 0    -p[3]  p[2];
    	 p[3]  0    -p[1];
    	-p[2]  p[1]  0;
    ]
end

function ∂skew∂p(λ) # 𝞉(skew(p)*λ)/∂p
	skew(-λ)
end
