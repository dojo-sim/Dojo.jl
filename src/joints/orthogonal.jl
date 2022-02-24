function orthogonal_rows(axis::AbstractVector)
    if norm(axis) > 0
        axis = normalize(axis)
    end
    A = svd(skew(axis)).Vt
    inds = SA[1; 2; 3]
    V1 = A[1,inds]'
    V2 = A[2,inds]'
    V3 = axis' # instead of A[3,:] for correct sign: abs(axis) = abs(A[3,:])

    return V1, V2, V3
end

function orthogonal_columns(axis::AbstractVector)
    V1, V2, V3 = orthogonal_rows(axis)
    return V1', V2', V3'
end