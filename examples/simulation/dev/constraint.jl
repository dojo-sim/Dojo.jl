


function sdf_capsule_capsule_constraint(xa, qa, ha, ra, xb, qb, hb, rb)
    sdf_capsule_capsule(xa, qa, ha, ra, xb, qb, hb, rb)
end

sdf_capsule_capsule_constraint(p1, o1, h1, r1, p2, o2, h2, r2)

function sdf_capsule_capsule_constraint_jacobian(relative, xa, qa, ha, ra, xb, qb, hb, rb; 
    attjac=false) 

    # get closest points 
    pa, pb = closest_points_capsule_capsule(xa, qa, ha, ra, xb, qb, hb, rb)
    d = pb - pa

    # Jacobians
    if relative == :parent 
        dx = FiniteDiff.finite_difference_jacobian(x -> closest_points_capsule_capsule(x, qa, ha, ra, xb, qb, hb, rb)[2], xa)
        dx -= FiniteDiff.finite_difference_jacobian(x -> closest_points_capsule_capsule(x, qa, ha, ra, xb, qb, hb, rb)[1], xa)
        dq = FiniteDiff.finite_difference_jacobian(q -> closest_points_capsule_capsule(xa, Quaternion(q...), ha, ra, xb, qb, hb, rb)[2], vector(qa))
        dq -= FiniteDiff.finite_difference_jacobian(q -> closest_points_capsule_capsule(xa, Quaternion(q...), ha, ra, xb, qb, hb, rb)[1], vector(qa))
        attjac && (Q *= LVᵀmat(qa))
    elseif relative == :child 
        dx = FiniteDiff.finite_difference_jacobian(x -> closest_points_capsule_capsule(xa, qa, ha, ra, x, qb, hb, rb)[2], xb)
        dx -= FiniteDiff.finite_difference_jacobian(x -> closest_points_capsule_capsule(xa, qa, ha, ra, x, qb, hb, rb)[1], xb)
        dq = FiniteDiff.finite_difference_jacobian(q -> closest_points_capsule_capsule(xa, qa, ha, ra, xb, Quaternion(q...), hb, rb)[2], vector(qb))
        dq -= FiniteDiff.finite_difference_jacobian(q -> closest_points_capsule_capsule(xa, qa, ha, ra, xb, Quaternion(q...), hb, rb)[1], vector(qb))
        attjac && (Q *= LVᵀmat(qb))
    end

    return d' ./ norm(d) * [dx dq]
end

Jp = constraint_jacobian(:parent, p1, o1, h1, r1, p2, o2, h2, r2, attjac=false)
Jc = constraint_jacobian(:child, p1, o1, h1, r1, p2, o2, h2, r2, attjac=false)

using FiniteDiff
Jpx = FiniteDiff.finite_difference_jacobian(x -> constraint(x, o1, h1, r1, p2, o2, h2, r2), p1)
Jpq = FiniteDiff.finite_difference_jacobian(q -> constraint(p1, Quaternion(q...), h1, r1, p2, o2, h2, r2), vector(o1))

Jcx = FiniteDiff.finite_difference_jacobian(x -> constraint(p1, o1, h1, r1, x, o2, h2, r2), p2)
Jcq = FiniteDiff.finite_difference_jacobian(q -> constraint(p1, o1, h1, r1, p2, Quaternion(q...), h2, r2), vector(o2))

norm(Jp - [Jpx Jpq])
norm(Jc - [Jcx Jcq])

function force_mapping(relative, xa, qa, ha, ra, xb, qb, hb, rb)
    pa, pb = closest_points_capsule_capsule(xa, qa, ha, ra, xb, qb, hb, rb)

    na = (pb - pa) ./ norm(pb - pa)
    if norm(na - [1.0; 0.0; 0.0]) > 1.0e-8
        ta1 = RotX(0.5 * π) * na # add edge case
    else
        ta1 = RotY(0.5 * π) * na 
    end
    ta2 = cross(na, ta1)

    if relative == :parent 
        return [ta1 ta2 na]
    elseif relative == :child 
        return -1.0 * [ta1 ta2 na]
    end
end

function force_mapping_jacobian(relative, jacobian, xa, qa, ha, ra, xb, qb, hb, rb, λ;
    attjac=false)
    if jacobian == :parent 
        dx = FiniteDiff.finite_difference_jacobian(x -> force_mapping(relative, x, qa, ha, ra, xb, qb, hb, rb) * λ, xa)
        dq = FiniteDiff.finite_difference_jacobian(q -> force_mapping(relative, xa, Quaternion(q...), ha, ra, xb, qb, hb, rb) * λ, vector(qa))
        attjac && (dq *= LVᵀmat(qa))
    elseif jacobian == :child 
        dx = FiniteDiff.finite_difference_jacobian(x -> force_mapping(relative, xa, qa, ha, ra, x, qb, hb, rb) * λ, xb)
        dq = FiniteDiff.finite_difference_jacobian(q -> force_mapping(relative, xa, qa, ha, ra, xb, Quaternion(q...), hb, rb) * λ, vector(qb))
        attjac && (dq *= LVᵀmat(qa))
    end
    return [dx dq]
end

force_mapping(:parent, p1, o1, h1, r1, p2, o2, h2, r2)
force_mapping(:child, p1, o1, h1, r1, p2, o2, h2, r2)

force_mapping_jacobian(:parent, :parent, p1, o1, h1, r1, p2, o2, h2, r2, ones(3), attjac=false)
force_mapping_jacobian(:parent, :child, p1, o1, h1, r1, p2, o2, h2, r2, ones(3), attjac=false)
force_mapping_jacobian(:child, :parent, p1, o1, h1, r1, p2, o2, h2, r2, ones(3), attjac=false)
force_mapping_jacobian(:child, :child, p1, o1, h1, r1, p2, o2, h2, r2, ones(3), attjac=false)

function torque_mapping(relative, xa, qa, ha, ra, xb, qb, hb, rb)
    pa, pb = closest_points_capsule_capsule(xa, qa, ha, ra, xb, qb, hb, rb)
    F = force_mapping(relative, xa, qa, ha, ra, xb, qb, hb, rb)

    if relative == :parent 
        ra = pa - xa 
        τa = rotation_matrix(inv(qa)) * skew(ra) * F
        return τa
    elseif relative == :child 
        rb = pb - xb 
        τb = rotation_matrix(inv(qb)) * skew(rb) * F
        return τb
    end
end

function torque_mapping_jacobian(relative, jacobian, xa, qa, ha, ra, xb, qb, hb, rb, λ;
    attjac=false)
    if jacobian == :parent 
        dx = FiniteDiff.finite_difference_jacobian(x -> torque_mapping(relative, x, qa, ha, ra, xb, qb, hb, rb) * λ, xa)
        dq = FiniteDiff.finite_difference_jacobian(q -> torque_mapping(relative, xa, Quaternion(q...), ha, ra, xb, qb, hb, rb) * λ, vector(qa))
        attjac && (dq *= LVᵀmat(qa))
    elseif jacobian == :child 
        dx = FiniteDiff.finite_difference_jacobian(x -> torque_mapping(relative, xa, qa, ha, ra, x, qb, hb, rb) * λ, xb)
        dq = FiniteDiff.finite_difference_jacobian(q -> torque_mapping(relative, xa, qa, ha, ra, xb, Quaternion(q...), hb, rb) * λ, vector(qb))
        attjac && (dq *= LVᵀmat(qa))
    end
    return [dx dq]
end

torque_mapping(:parent, p1, o1, h1, r1, p2, o2, h2, r2)
torque_mapping(:child, p1, o1, h1, r1, p2, o2, h2, r2)

torque_mapping_jacobian(:parent, :parent, p1, o1, h1, r1, p2, o2, h2, r2, ones(3), attjac=false)
torque_mapping_jacobian(:parent, :child, p1, o1, h1, r1, p2, o2, h2, r2, ones(3), attjac=false)
torque_mapping_jacobian(:child, :parent, p1, o1, h1, r1, p2, o2, h2, r2, ones(3), attjac=false)
torque_mapping_jacobian(:child, :child, p1, o1, h1, r1, p2, o2, h2, r2, ones(3), attjac=false)

function impulse_map(relative, xa, qa, ha, ra, xb, qb, hb, rb)
    [
        force_mapping(relative, xa, qa, ha, ra, xb, qb, hb, rb);
        torque_mapping(relative, xa, qa, ha, ra, xb, qb, hb, rb);
    ]
end

function impulse_map_jacobian_configuration(relative, jacobian, xa, qa, ha, ra, xb, qb, hb, rb, λ)
    [
        force_mapping_jacobian(relative, jacobian, xa, qa, ha, ra, xb, qb, hb, rb, λ);
        torque_mapping_jacobian(relative, jacobian, xa, qa, ha, ra, xb, qb, hb, rb, λ);
    ]
end

function tangential_velocity(x3a, q3a, x2a, q2a, ha, ra, x3b, q3b, x2b, q2b, hb, rb, h)
    # contact point kinematics
    ka, kb = kinematics_capsule_capsule(x3a, q3a, ha, ra, x3b, q3b, hb, rb)
    # capsule a points
    p3a = x3a + vector_rotate(ka, q3a)
    p3b = x3b + vector_rotate(kb, q3b)
    # capsule b points
    p2a = x2a + vector_rotate(ka, q2a) 
    p2b = x2b + vector_rotate(kb, q2b)
    # velocities
    va = (p3a - p2a) ./ h
    vb = (p3b - p2b) ./ h 
    # relative velocity
    Δv = vb - va 
    # projection matrix 
    na = (p3b - p3a) ./ norm(p3b - p3a)
    if norm(na - [1.0; 0.0; 0.0]) > 1.0e-8
        ta1 = RotX(0.5 * π) * na # add edge case
    else
        ta1 = RotY(0.5 * π) * na 
    end
    ta2 = cross(na, ta1)
    P = [ta1'; ta2'] 
    # tangential velocity
    return P * Δv
end

h = 0.1
x_shift = [0.0; 0.0; 1.0]
tangential_velocity(p1, o1, p1 + x_shift * h, o1, h1, r1, p2, o2, p2, o2, h2, r2, h)

## constraints 
function constraint(model, s::AbstractVector{T}, γ::AbstractVector{T},
    x3a::AbstractVector{T}, q3a::Quaternion{T}, x2a::AbstractVector{T}, q2a::Quaternion{T},
    x3b::AbstractVector{T}, q3b::Quaternion{T}, x2b::AbstractVector{T}, q2b::Quaternion{T},
    h::T) where T

    # transforms the velocities of the origin of the link into velocities
    vp = tangential_velocity(x3a, q3a, x2a, q2a, h1, r1, x3b, q3b, x2b, q2b, h2, r2, h)

    SVector{4,T}(
        sdf_capsule_capsule(x3a, q3a, ha, ra, x3b, q3b, hb, rb),
        model.friction_coefficient * γ[1] - γ[2],
        (vp - s[@SVector [3,4]])...)
end

