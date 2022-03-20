"""
    Collision 

    abstract type defining interaction between two bodies
"""
abstract type Collision{T,O,I,OI} end 

# contact point origin 
function contact_point_origin(x, q, k) 
    x + vector_rotate(k, q)
end

∂contact_point_origin∂x(x, q, k) = 1.0 * I(3)
∂contact_point_origin∂q(x, q, k) = ∂vector_rotate∂q(k, q)
∂contact_point_origin∂k(x, q, k) = ∂vector_rotate∂p(k, q)

# normal projection (from child to parent)
function contact_normal(collision::Collision, xp, qp, xc, qc)
    # contact points
    cop = contact_point(:parent, collision, xp, qp, xc, qc) 
    coc = contact_point(:child,  collision, xp, qp, xc, qc)
 
    # unnormalized direction 
    dir = cop - coc

    # distance 
    dis = distance(collision, xp, qp, xc, qc)

    # normalized direction
    if dis >= 0.0
        return normalize(dir)'
    else 
        return -1.0 * normalize(dir)'
    end
end

function ∂contact_normal_transpose∂x(jacobian::Symbol, collision::Collision, xp, qp, xc, qc)
    # contact origin points
    cop = contact_point(:parent, collision, xp, qp, xc, qc) 
    coc = contact_point(:child,  collision, xp, qp, xc, qc)

    # unnormalized direction 
    dir = cop - coc

    # Jacobians
    X = ∂normalize∂x(dir) * (∂contact_point∂x(:parent, jacobian, collision, xp, qp, xc, qc) - ∂contact_point∂x(:child, jacobian, collision, xp, qp, xc, qc))
    
    if jacobian == :parent 
        FD = FiniteDiff.finite_difference_jacobian(x -> contact_normal(collision, x, qp, xc, qc)', xp)
    elseif jacobian == :child 
        FD = FiniteDiff.finite_difference_jacobian(x -> contact_normal(collision, xp, qp, x, qc)', xc)
    end

    # distance 
    dis = distance(collision, xp, qp, xc, qc)

    @assert norm((dis >= 0.0 ? 1.0 : -1.0) * X - FD, Inf) < 1.0e-4

    # normalized direction
    if dis >= 0.0
        return X
    else 
        return -1.0 * X
    end
end

function ∂contact_normal_transpose∂q(jacobian::Symbol, collision::Collision, xp, qp, xc, qc)
    # contact origin points
    cop = contact_point(:parent, collision, xp, qp, xc, qc) 
    coc = contact_point(:child,  collision, xp, qp, xc, qc)

    # unnormalized direction 
    dir = cop - coc 

    Q = ∂normalize∂x(dir) * (∂contact_point∂q(:parent, jacobian, collision, xp, qp, xc, qc) - ∂contact_point∂q(:child, jacobian, collision, xp, qp, xc, qc))

    # Jacobians
    if jacobian == :parent 
        FD = FiniteDiff.finite_difference_jacobian(q -> contact_normal(collision, xp, UnitQuaternion(q..., false), xc, qc)', vector(qp))
    elseif jacobian == :child 
        FD = FiniteDiff.finite_difference_jacobian(q -> contact_normal(collision, xp, qp, xc, UnitQuaternion(q..., false))', vector(qc))
    end

    # distance 
    dis = distance(collision, xp, qp, xc, qc)

    @assert norm((dis >= 0.0 ? 1.0 : -1.0) * Q - FD, Inf) < 1.0e-4

    # normalized direction
    if dis >= 0.0
        return Q
    else 
        return -1.0 * Q
    end
end

# tangent projection (in child frame)
function contact_tangent(collision::Collision, xp, qp, xc, qc)
    # normal
    n = vec(contact_normal(collision, xp, qp, xc, qc))

    # tangents
    v1 = [1.0; 0.0; 0.0] # candidate
    v2 = skew(v1) * n # tangent
    if !(norm(v2) > 1.0e-6)
        @warn "edge case!"
        v1 = [0.0; 1.0; 0.0] # candidate
        v2 = skew(v1) * n # tangent
    end
    v3 = skew(v2) * n # tangent

    return [v2'; v3']
    # return szeros(0, 3)
end

# contact_tangent * λ
function ∂contact_tangent_jvp∂x(jacobian::Symbol, collision::Collision, xp, qp, xc, qc, λ)
    @assert length(λ) == 3
    
    if jacobian == :parent 
        T1 = λ' * FiniteDiff.finite_difference_jacobian(x -> contact_tangent(collision, x, qp, xp, qc)[1, :]', xp)
        T2 = λ' * FiniteDiff.finite_difference_jacobian(x -> contact_tangent(collision, x, qp, xp, qc)[2, :]', xp)
        return [T1; T2]
    elseif jacobian == :child 
        T1 = λ' * FiniteDiff.finite_difference_jacobian(x -> contact_tangent(collision, xp, qp, x, qc)[1, :]', xc)
        T2 = λ' * FiniteDiff.finite_difference_jacobian(x -> contact_tangent(collision, xp, qp, x, qc)[2, :]', xc)
        return [T1; T2]
    end
end

# λ' * contact_tangent
function ∂contact_tangent_vjp∂x(jacobian::Symbol, collision::Collision, xp, qp, xc, qc, λ)
    @assert length(λ) == 2 || length(λ) == 4 || length(λ) == 0
    if jacobian == :parent 
        return FiniteDiff.finite_difference_jacobian(x -> (λ' * contact_tangent(collision, x, qp, xc, qc))', xp)'
    elseif jacobian == :child 
        return FiniteDiff.finite_difference_jacobian(x -> (λ' * contact_tangent(collision, xp, qp, x, qc))', xc)'
    end
end

# contact_tangent * λ
function ∂contact_tangent_jvp∂q(jacobian::Symbol, collision::Collision, xp, qp, xc, qc, λ)
    @assert length(λ) == 3
    if jacobian == :parent 
        T1 = λ' * FiniteDiff.finite_difference_jacobian(q -> contact_tangent(collision, xp, UnitQuaternion(q..., false), xc, qc)[1, :]', vector(qp))
        T2 = λ' * FiniteDiff.finite_difference_jacobian(q -> contact_tangent(collision, xp, UnitQuaternion(q..., false), xc, qc)[2, :]', vector(qp))
        return [T1; T2]
    elseif jacobian == :child 
        T1 = λ' * FiniteDiff.finite_difference_jacobian(q -> contact_tangent(collision, xp, qp, xc, UnitQuaternion(q..., false))[1, :]', vector(qc))
        T2 = λ' * FiniteDiff.finite_difference_jacobian(q -> contact_tangent(collision, xp, qp, xc, UnitQuaternion(q..., false))[2, :]', vector(qc))
        return [T1; T2]
    end
end

# λ' * contact_tangent
function ∂contact_tangent_vjp∂q(jacobian::Symbol, collision::Collision, xp, qp, xc, qc, λ)
    @assert length(λ) == 2 || length(λ) == 4 || length(λ) == 0
    if jacobian == :parent 
        return FiniteDiff.finite_difference_jacobian(q -> (λ' * contact_tangent(collision, xp, UnitQuaternion(q..., false), xc, qc))', vector(qp))'
    elseif jacobian == :child 
        return FiniteDiff.finite_difference_jacobian(q -> (λ' * contact_tangent(collision, xp, qp, xc, UnitQuaternion(q..., false)))', vector(qc))'
    end
end


