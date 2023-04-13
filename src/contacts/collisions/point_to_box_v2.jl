function klamp(t, a, b) 
    if t <= a 
        return a 
    elseif b <= t 
        return b 
    else
        da = abs(t - a) 
        db = abs(t - b)

        if da < db 
            return a 
        else
            return b
        end
        # return t
    end
end

function contact_point_box(ps, xc, qc, kx, ky, kz)
    # rotate point into box frame
    ps_box_frame = vector_rotate(ps - xc, inv(qc)) 

    # find box point in box frame
    px = clamp(ps_box_frame[1], -0.5 * kx, 0.5 * kx) 
    py = clamp(ps_box_frame[2], -0.5 * ky, 0.5 * ky)
    pz = clamp(ps_box_frame[3], -0.5 * kz, 0.5 * kz)

    # box point in world frame
    return vector_rotate([px; py; pz], qc) + xc
end

function ∂contact_point_box∂p(ps, xc, qc, kx, ky, kz)
    
    FiniteDiff.finite_difference_jacobian(p -> contact_point_box(p, xc, qc, kx, ky, kz), ps)
end

function ∂contact_point_box∂x(ps, xc, qc, kx, ky, kz)
    
    FiniteDiff.finite_difference_jacobian(x -> contact_point_box(ps, x, qc, kx, ky, kz), xc)
end

function ∂contact_point_box∂q(ps, xc, qc, kx, ky, kz)
    
    FiniteDiff.finite_difference_jacobian(q -> contact_point_box(ps, xc, Quaternion(q...), kx, ky, kz), vector(qc))
end
