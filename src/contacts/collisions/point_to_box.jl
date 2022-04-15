function contact_point_box(ps, xc, qc, kx, ky, kz)
    # corner kinematics
    cx = contact_point_origin(szeros(3), qc, kx)
    cy = contact_point_origin(szeros(3), qc, ky)
    cz = contact_point_origin(szeros(3), qc, kz)

    v1 = 2.0 * cx 
    v2 = 2.0 * cy 
    v3 = 2.0 * cz

    # origin is opposite corner
    origin = xc - cx - cy - cz

    tx = dot(ps - origin, v1) / dot(v1, v1) 
    ty = dot(ps - origin, v2) / dot(v2, v2) 
    tz = dot(ps - origin, v3) / dot(v3, v3)

    coc = origin
    
    if tx > 1.0 
        coc += v1 
    elseif tx < 0.0 
        nothing 
    else 
        coc += tx * v1 
    end

    if ty > 1.0 
        coc += v2
    elseif ty < 0.0 
        nothing 
    else 
        coc += ty * v2
    end

    if tz > 1.0 
        coc += v3 
    elseif tz < 0.0 
        nothing 
    else 
        coc += tz * v3 
    end

    return coc
end


function ∂contact_point_box∂p(ps, xc, qc, kx, ky, kz)
    cx = contact_point_origin(szeros(3), qc, kx)
    cy = contact_point_origin(szeros(3), qc, ky)
    cz = contact_point_origin(szeros(3), qc, kz)

    v1 = 2.0 * cx 
    v2 = 2.0 * cy 
    v3 = 2.0 * cz

    origin = xc - cx - cy - cz

    tx = dot(ps - origin, v1) / dot(v1, v1) 
    ty = dot(ps - origin, v2) / dot(v2, v2) 
    tz = dot(ps - origin, v3) / dot(v3, v3)

    P = szeros(3, 3)
    
    if tx > 1.0 
        # coc += v1 
        # X += 
    elseif tx < 0.0 
        # nothing 
        # X += 
    else 
        # coc += tx * v1 
        P += v1 * 1.0 / dot(v1, v1) * v1' 
    end

    if ty > 1.0 
        # coc += v2
    elseif ty < 0.0 
        # nothing 
    else 
        # coc += ty * v2
        P += v2 * 1.0 / dot(v2, v2) * v2'
    end

    if tz > 1.0 
        # coc += v3 
    elseif tz < 0.0 
        # nothing 
    else 
        # coc += tz * v3
        P += v3 * 1.0 / dot(v3, v3) * v3' 
    end

    return P
    # FiniteDiff.finite_difference_jacobian(p -> contact_point_box(p, xc, qc, kx, ky, kz), ps)
end

function ∂contact_point_box∂x(ps, xc, qc, kx, ky, kz)
    cx = contact_point_origin(szeros(3), qc, kx)
    cy = contact_point_origin(szeros(3), qc, ky)
    cz = contact_point_origin(szeros(3), qc, kz)

    v1 = 2.0 * cx 
    v2 = 2.0 * cy 
    v3 = 2.0 * cz

    origin = xc - cx - cy - cz

    tx = dot(ps - origin, v1) / dot(v1, v1) 
    ty = dot(ps - origin, v2) / dot(v2, v2) 
    tz = dot(ps - origin, v3) / dot(v3, v3)

    X = 1.0 * I(3)
    
    if tx > 1.0 
        # coc += v1 
    elseif tx < 0.0 
        # nothing 
    else 
        # coc += tx * v1 
        X += v1 * 1.0 / dot(v1, v1) * v1' * -1.0
    end

    if ty > 1.0 
        # coc += v2
    elseif ty < 0.0 
        # nothing 
    else 
        # coc += ty * v2
        X += v2 * 1.0 / dot(v2, v2) * v2' * -1.0
    end

    if tz > 1.0 
        # coc += v3 
    elseif tz < 0.0 
        # nothing 
    else 
        # coc += tz * v3
        X += v3 * 1.0 / dot(v3, v3) * v3' * -1.0 
    end

    return X
    # FiniteDiff.finite_difference_jacobian(x -> contact_point_box(ps, x, qc, kx, ky, kz), xc)
end

function ∂contact_point_box∂q(ps, xc, qc, kx, ky, kz)
    cx = contact_point_origin(szeros(3), qc, kx)
    cy = contact_point_origin(szeros(3), qc, ky)
    cz = contact_point_origin(szeros(3), qc, kz)

    v1 = 2.0 * cx 
    v2 = 2.0 * cy 
    v3 = 2.0 * cz

    ∂v1∂q = 2.0 * ∂contact_point_origin∂q(szeros(3), qc, kx)
    ∂v2∂q = 2.0 * ∂contact_point_origin∂q(szeros(3), qc, ky)
    ∂v3∂q = 2.0 * ∂contact_point_origin∂q(szeros(3), qc, kz)

    origin = xc - cx - cy - cz

    tx = dot(ps - origin, v1) / dot(v1, v1) 
    ty = dot(ps - origin, v2) / dot(v2, v2) 
    tz = dot(ps - origin, v3) / dot(v3, v3)

    Q = -0.5 * (∂v1∂q + ∂v2∂q + ∂v3∂q)

    if tx > 1.0 
        # coc += v1 
        Q += ∂v1∂q
    elseif tx < 0.0 
        # nothing 
    else 
        # coc += tx * v1 
        Q += tx * ∂v1∂q
        Q += v1 * (ps - xc)' ./ dot(v1, v1) * ∂v1∂q
        Q += v1 * -1.0 * dot(ps - xc, v1) ./ dot(v1, v1)^2 * 2.0 * v1' * ∂v1∂q
        Q += v1 * 1.0 / dot(v1, v1) * v1' * -0.5 * ∂v1∂q
    end

    if ty > 1.0 
        # coc += v2
        Q += ∂v2∂q
    elseif ty < 0.0 
        # nothing 
    else 
        # coc += ty * v2
        Q += tx * ∂v2∂q
        Q += v2 * (ps - xc)' ./ dot(v2, v2) * ∂v2∂q
        Q += v2 * -1.0 * dot(ps - xc, v2) ./ dot(v2, v2)^2 * 2.0 * v2' * ∂v2∂q
        Q += v2 * 1.0 / dot(v2, v2) * v2' * -0.5 * ∂v2∂q
    end

    if tz > 1.0 
        # coc += v3 
        Q += ∂v3∂q
    elseif tz < 0.0 
        # nothing 
    else 
        # coc += tz * v3 
        Q += tx * ∂v3∂q
        Q += v3 * (ps - xc)' ./ dot(v3, v3) * ∂v3∂q
        Q += v3 * -1.0 * dot(ps - xc, v3) ./ dot(v3, v3)^2 * 2.0 * v3' * ∂v3∂q
        Q += v3 * 1.0 / dot(v3, v3) * v3' * -0.5 * ∂v3∂q
    end

    return Q
    # FiniteDiff.finite_difference_jacobian(q -> contact_point_box(ps, xc, Quaternion(q..., false), kx, ky, kz), vector(qc))
end