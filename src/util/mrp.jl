function mrp(q::UnitQuaternion)
    q̄ = vector(q)
    return q̄[2:4] ./ (q̄[1] + 1.0) 
end

function mrp(q::AbstractVector)
    return q[2:4] ./ (q[1] + 1.0)
end

function ∂mrp∂q(q::AbstractVector)
    s = q[1] 
    v = q[2:4]

    d1 = 1 / (s + 1)^2
    di = 1 / (s + 1)

    [-v[1] * d1  di 0.0 0.0; 
     -v[2] * d1 0.0 di 0.0;
     -v[3] * d1 0.0 0.0 di]
end

function axis(q)
    m = mrp(q)
    mag = norm(m)
    if mag > 0
        # θ = 4.0 * atan(mag) 
        n = m ./ mag
    else 
        # θ = 0.0
        n = [1.0; 0.0; 0.0]
    end
    return n
end

function ∂axis∂q(q::AbstractVector)
    m = mrp(q)
    ∂mrp∂q(q) ./ norm(m) - m ./ norm(m)^2 * transpose(m ./ norm(m)) * ∂mrp∂q(q)
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

function ∂angle∂q(q::AbstractVector) 
    m = mrp(q)
    4.0 * 1.0 / (1.0 + norm(m)^2.0) * transpose(m ./ norm(m)) * ∂mrp∂q(q)
end

function axis_angle(q)
    return axis(q), angle(q)
end

function rotation_vector(q)
    n, θ = axis_angle(q)
    return θ .* n 
end

function ∂rotation_vector∂q(q::AbstractVector)
    θ = angle(q)
    if θ != 0.0
        a = axis(q)
        return a * ∂angle∂q(q) + θ * ∂axis∂q(q) 
    else
        FiniteDiff.finite_difference_jacobian(rotation_vector, vector(q))
    end
end

∂rotation_vector∂q(q::UnitQuaternion) = ∂rotation_vector∂q(vector(q))
