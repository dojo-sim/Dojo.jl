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

function ∂Jp∂J(p) #∂(J*p)/∂flatten(J)
    SA[
        p[1]  p[2]  p[3]  0     0     0;
        0     p[1]  0     p[2]  p[3]  0;
        0     0     p[1]  0     p[2]  p[3];
    ]
end

function attitude_jacobian(data::AbstractVector, Nb::Int)
    G = zeros(0,0)
    for i = 1:Nb
        x2, v15, q2, ω15 = unpack_data(data[13 * (i-1) .+ (1:13)])
        q2 = Quaternion(q2...)
        G = cat(G, sI(6), LVᵀmat(q2), sI(3), dims = (1,2))
    end
    ndata = length(data)
    nu = ndata - size(G)[1]
    G = cat(G, I(nu), dims = (1,2))
    return G
end

function unpack_data(data::AbstractVector)
    x2 = data[SA[1;2;3]]
    v15 = data[SA[4;5;6]]
    q2 = data[SA[7;8;9;10]]
    ω15 = data[SA[11;12;13]]
    return x2, v15, q2, ω15
end






