function initialize!(collision::CapsuleCapsuleCollision, xa, qa, xb, qb)
    # parent points
    ha = collision.height_parent 
    pa1 = xa + rotation_matrix(qa) * [0.0; 0.0;  0.5 * ha]
    pa2 = xa + rotation_matrix(qa) * [0.0; 0.0; -0.5 * ha]

    # child points
    hb = collision.height_child
    pb1 = xb + rotation_matrix(qb) * [0.0; 0.0;  0.5 * hb]
    pb2 = xb + rotation_matrix(qb) * [0.0; 0.0; -0.5 * hb]

    # initialize variables
    collision.ip.z .= [0.5; 0.5; 1.0; 1.0; 0.0; 0.0; 1.0; 1.0; 1.0; 1.0] #TODO: allocation free
    
    # initialize data
    collision.ip.θ .= [pa1; pa2; pb1; pb2] #TODO: allocation free

    return 
end 

function line_segment_points(collision::CapsuleCapsuleCollision, x, q, h)
    p1 = x + rotation_matrix(q) * [0.0; 0.0;  0.5 * h]
    p2 = x + rotation_matrix(q) * [0.0; 0.0; -0.5 * h]
    return p1, p2 
end

function capsule_contact_origin(relative::Symbol, collision::CapsuleCapsuleCollision, xa, qa, xb, qb) 
    # initialize solver
    initialize!(collision, xa, qa, xb, qb) 

    # optimize to find closest points on line segments
    interior_point_solve!(collision.ip)

    if relative == :parent 
        ta = collision.ip.z[1]
        ha = collision.height_parent
        pa1, pa2 = line_segment_points(collision, xa, qa, ha)
        pa = pa1 + ta * (pa2 - pa1)
        return pa
    elseif relative == :child
        tb = collision.ip.z[2]
        hb = collision.height_child
        pb1, pb2 = line_segment_points(collision, xb, qb, hb)
        pb = pb1 + tb * (pb2 - pb1)
        return pb
    elseif relative == :difference
        ta = collision.ip.z[1]
        ha = collision.height_parent
        pa1, pa2 = line_segment_points(collision, xa, qa, ha)
        pa = pa1 + ta * (pa2 - pa1)
        pa = pa1 + ta * (pa2 - pa1)
        tb = collision.ip.z[2]
        hb = collision.height_child
        pb1, pb2 = line_segment_points(collision, xb, qb, hb)
        pb = pb1 + tb * (pb2 - pb1)
        return pa - pb
    end
end

function ∂capsule_contact_origin∂x(relative::Symbol, jacobian::Symbol, collision::CapsuleCapsuleCollision, xa, qa, xb, qb)

end
   
function ∂capsule_contact_origin∂q(relative::Symbol, jacobian::Symbol, collision::CapsuleCapsuleCollision, xa, qa, xb, qb)

end

# distance
function distance(collision::CapsuleCapsuleCollision, xp, qp, xc, qc)
    # difference
    d = capsule_contact_origin(:difference, collision, xp, qp, xc, qc)
     
    return norm(d, 2) - (collision.radius_parent + collision.radius_child)
end

function ∂distance∂x(gradient::Symbol, collision::CapsuleCapsuleCollision, xp, qp, xc, qc)
    
end

function ∂distance∂q(gradient::Symbol, collision::CapsuleCapsuleCollision, xp, qp, xc, qc)
   
end

# contact point in world frame
function contact_point(relative::Symbol, collision::CapsuleCapsuleCollision, xp, qp, xc, qc) 
    
end

function ∂contact_point∂x(relative::Symbol, jacobian::Symbol, collision::CapsuleCapsuleCollision, xp, qp, xc, qc)
   
end

function ∂contact_point∂q(relative::Symbol, jacobian::Symbol, collision::CapsuleCapsuleCollision, xp, qp, xc, qc)

end

function contact_normal(collision::CapsuleCapsuleCollision, xp, qp, xc, qc)
    # unnormalized direction 
    dir = capsule_contact_origin(:difference, collision, xp, qp, xc, qc)

    # distance 
    dis = norm(dir, 2) - (collision.radius_parent + collision.radius_child)

    # normalized direction
    if dis >= 0.0
        return normalize(dir)'
    else 
        return -1.0 * normalize(dir)'
    end
end

function ∂contact_normal_transpose∂x(jacobian::Symbol, collision::CapsuleCapsuleCollision, xp, qp, xc, qc)
    # contact origin points
    cop = 
    coc = 

    # unnormalized direction 
    dir = cop - coc

    # Jacobians
    X = ∂normalize∂x(dir) * (∂∂x() - ∂∂x())
    
    # distance 
    dis = distance(collision, xp, qp, xc, qc)

    # normalized direction
    if dis >= 0.0
        return X
    else 
        return -1.0 * X
    end
end

function ∂contact_normal_transpose∂q(jacobian::Symbol, collision::CapsuleCapsuleCollision, xp, qp, xc, qc)
    # contact origin points
    cop = 
    coc = 

    # unnormalized direction 
    dir = cop - coc 

    Q = ∂normalize∂x(dir) * (∂∂q() - ∂∂q())

    # distance 
    dis = distance(collision, xp, qp, xc, qc)

    # normalized direction
    if dis >= 0.0
        return Q
    else 
        return -1.0 * Q
    end
end
