function mrp(q::Quaternion)
    q̄ = vector(q)
    return q̄[SUnitRange(2,4)] ./ (q̄[1] + 1.0)
end

function mrp(q::AbstractVector)
    return q[SUnitRange(2,4)] ./ (q[1] + 1.0)
end

function dmrpdq(q::AbstractVector)
    s = q[1]
    v = q[SA[2;3;4]]

    d1 = 1 / (s + 1)^2
    di = 1 / (s + 1)

    [-v[1] * d1  di 0.0 0.0;
     -v[2] * d1  0.0 di 0.0;
     -v[3] * d1  0.0 0.0 di]
end

function axis(q)
    m = mrp(q)
    mag = norm(m)
    if mag > 0
        # θ = 4.0 * atan(mag)
        n = m ./ mag
    else
        # θ = 0.0
        n = SVector{3}(1.0, 0.0, 0.0)
    end
    return n
end

function daxisdq(q::AbstractVector)
    m = mrp(q)
    dmrpdq(q) ./ norm(m) - m ./ norm(m)^2 * transpose(m ./ norm(m)) * dmrpdq(q)
end


function angle(q)
    m = mrp(q)
    mag = norm(m)
    if mag > 0
        θ = 4.0 * atan(mag)
    else
        θ = 0.0
    end
    return θ
end

function dangledq(q::AbstractVector)
    m = mrp(q)
    4.0 * 1.0 / (1.0 + norm(m)^2.0) * transpose(m ./ norm(m)) * dmrpdq(q)
end

function axis_angle(q)
    return axis(q), angle(q)
end

function rotation_vector(q)
    n, θ = axis_angle(q)
    return θ .* n
end

function drotation_vectordq(q::AbstractVector)
    θ = angle(q)
    if θ != 0.0
        a = axis(q)
        return a * dangledq(q) + θ * daxisdq(q)
    else
        [
            0.0  2.0  0.0  0.0;
            0.0  0.0  2.0  0.0;
            0.0  0.0  0.0  2.0;
        ]
    end
end

drotation_vectordq(q::Quaternion) = drotation_vectordq(vector(q))
