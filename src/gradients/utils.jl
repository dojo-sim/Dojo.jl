################################################################################
# Inertia
################################################################################
function lift_inertia(j::SVector{6,T}) where T
    J = SMatrix{3,3,T,9}(
        [j[1] j[2] j[3];
         j[2] j[4] j[5];
         j[3] j[5] j[6]])
end
function flatten_inertia(J::SMatrix{3,3,T,9}) where T
    j = SVector{6,T}([J[1,1], J[1,2], J[1,3], J[2,2], J[2,3], J[3,3]])
end
function ∂inertia(p) #∂(J*p)/∂flatten(J)
    SA[
        p[1]  p[2]  p[3]  0     0     0;
        0     p[1]  0     p[2]  p[3]  0;
        0     0     p[1]  0     p[2]  p[3];
    ]
end


function getλJoint(joint::JointConstraint{T,N,Nc}, i::Int) where {T,N,Nc}
    n1 = 1
    for j = 1:i-1
        n1 += impulses_length([joint.translational, joint.rotational][j])
    end
    n2 = n1 - 1 + impulses_length([joint.translational, joint.rotational][i])

    λi = SVector{n2-n1+1,T}(joint.impulses[2][n1:n2])
    return λi
end


function attitude_jacobian(data::AbstractVector, Nb::Int)
    G = zeros(0,0)
    for i = 1:Nb
        x2, v15, q2, ϕ15 = unpack_data(data[13*(i-1) .+ (1:13)])
        q2 = UnitQuaternion(q2..., false)
        G = cat(G, I(6), LVᵀmat(q2), I(3), dims = (1,2))
    end
    ndata = length(data)
    nu = ndata - size(G)[1]
    G = cat(G, I(nu), dims = (1,2))
    return G
end

function unpack_data(data::AbstractVector)
    x2 = data[SVector{3,Int}(1:3)]
    v15 = data[SVector{3,Int}(4:6)]
    q2 = data[SVector{4,Int}(7:10)]
    ϕ15 = data[SVector{3,Int}(11:13)]
    return x2, v15, q2, ϕ15
end



################################################################################
# Index and Dimensions
################################################################################
function control_dimension(joint::JointConstraint{T,N,Nc}; ignore_floating_base::Bool = false) where {T,N,Nc}
    ignore_floating_base && (N == 0) && return 0
    N̄ = 0
    for (i, element) in enumerate([joint.translational, joint.rotational])
        N̄ += control_dimension(element)
    end
    return N̄
end

function control_dimension(joint::Joint{T,N}) where {T,N}
    return 3 - N
end

function control_dimension(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}; ignore_floating_base::Bool = false) where {T,Nn,Ne,Nb,Ni}
    nu = 0
    for joint in mechanism.joints
        nu += control_dimension(joint, ignore_floating_base = ignore_floating_base)
    end
    return nu
end

function minimal_dimension(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}) where {T,Nn,Ne,Nb,Ni}
    nx = 0
    free_rot_base = false # we are going to check if the link attached to the base has free orientation
    nx = 2 * control_dimension(mechanism, ignore_floating_base = false)
    free_rot_base && (nx += 1)
    return nx
end

maximal_dimension(mechanism::Mechanism{T,Nn,Ne,Nb}; attjac::Bool=false) where {T,Nn,Ne,Nb} = attjac ? 12Nb : 13Nb
