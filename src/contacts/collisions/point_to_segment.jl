# capsule contact point origin sphere-capsule
function contact_point_segment(ps, xc, qc, ka, kb)
    # capsule origin a 
    ca = contact_point_origin(xc, qc, ka) 

    # capsule origin b 
    cb = contact_point_origin(xc, qc, kb) 

    # point differences
    dab = cb - ca
    das = ps - ca

    t = dot(das, dab) / dot(dab, dab)
    # t_clamp = min(max(t, 0.0), 1.0)

    # closest point 
    if t <= 0.0 
        p = ca 
    elseif t >= 1.0 
        p = cb 
    else
        p = ca + t * dab
    end

    return p
end

# capsule contact point origin sphere-capsule
function ∂contact_point_segment∂p(ps, xc, qc, ka, kb)
    # capsule origin a 
    ca = contact_point_origin(xc, qc, ka) 

    # capsule origin b 
    cb = contact_point_origin(xc, qc, kb) 

    # point differences
    dab = cb - ca
    das = ps - ca

    t = dot(das, dab) ./ dot(dab, dab)
    # t_clamp = min(max(t, 0.0), 1.0)

    # closest point 
    if t <= 0.0 
        # p = ca 
        return szeros(3, 3)
    elseif t >= 1.0 
        # p = cb 
        return szeros(3, 3)
    else
        # p = ca + t * dab
        return dab * dab' ./ dot(dab, dab)
    end

    return p
end

# capsule contact point origin sphere-capsule
function ∂contact_point_segment∂x(ps, xc, qc, ka, kb)
    # capsule origin a 
    ca = contact_point_origin(xc, qc, ka) 

    # capsule origin b 
    cb = contact_point_origin(xc, qc, kb) 

    # point differences
    dab = cb - ca
    das = ps - ca

    t = dot(das, dab) ./ dot(dab, dab)
    # t_clamp = min(max(t, 0.0), 1.0)

    # closest point 
    if t <= 0.0 
        # p = ca 
        return ∂contact_point_origin∂x(xc, qc, ka) 
    elseif t >= 1.0 
        # p = cb 
        return ∂contact_point_origin∂x(xc, qc, kb) 
    else
        # p = ca + t * dab
        X = ∂contact_point_origin∂x(xc, qc, ka) 
        X += t * (∂contact_point_origin∂x(xc, qc, kb) - ∂contact_point_origin∂x(xc, qc, ka)) 

        t1 = 1.0 / (dab' * dab) * dab' * (-∂contact_point_origin∂x(xc, qc, ka))
        t2 = 1.0 / (dab' * dab) * das' * (∂contact_point_origin∂x(xc, qc, kb) - ∂contact_point_origin∂x(xc, qc, ka))
        t3 = -2.0 * das' * dab / (dab' * dab)^2 * dab' * (∂contact_point_origin∂x(xc, qc, kb) - ∂contact_point_origin∂x(xc, qc, ka)) 

        X += dab * (t1 + t2 + t3)
        return X
    end
end

# capsule contact point origin sphere-capsule
function ∂contact_point_segment∂q(ps, xc, qc, ka, kb)
    # capsule origin a 
    ca = contact_point_origin(xc, qc, ka) 

    # capsule origin b 
    cb = contact_point_origin(xc, qc, kb) 

    # point differences
    dab = cb - ca
    das = ps - ca

    t = dot(das, dab) ./ dot(dab, dab)
    # t_clamp = min(max(t, 0.0), 1.0)

    # closest point 
    if t <= 0.0 
        # p = ca 
        return ∂contact_point_origin∂q(xc, qc, ka) 
    elseif t >= 1.0 
        # p = cb 
        return ∂contact_point_origin∂q(xc, qc, kb) 
    else
        # p = ca + t * dab
        Q = ∂contact_point_origin∂q(xc, qc, ka) 
        Q += t * (∂contact_point_origin∂q(xc, qc, kb) - ∂contact_point_origin∂q(xc, qc, ka)) 

        t1 = 1.0 / (dab' * dab) * dab' * (-∂contact_point_origin∂q(xc, qc, ka))
        t2 = 1.0 / (dab' * dab) * das' * (∂contact_point_origin∂q(xc, qc, kb) - ∂contact_point_origin∂q(xc, qc, ka))
        t3 = -2.0 * das' * dab ./ (dab' * dab)^2 * dab' * (∂contact_point_origin∂q(xc, qc, kb) - ∂contact_point_origin∂q(xc, qc, ka)) 
        Q += dab * (t1 + t2 + t3)

        return Q
    end
end