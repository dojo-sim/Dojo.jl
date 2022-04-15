# capsule contact point origin sphere-capsule
function contact_point_segment(ps, xc, qc, ka, kb)
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
        p = ca 
    elseif t >= 1.0 
        p = cb 
    else
        p = ca + t * dab
    end

    return p
end