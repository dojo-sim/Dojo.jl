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
