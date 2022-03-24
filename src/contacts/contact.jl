"""
    Contact{T,N} 

    Abstract type containing contact information associated with Body objects.
"""
abstract type Contact{T,N} end

# constraint Jacobians
function constraint_jacobian_configuration(relative::Symbol, model::Contact, 
    xp::AbstractVector, vp::AbstractVector, qp::Quaternion, ϕp::AbstractVector, 
    xc::AbstractVector, vc::AbstractVector, qc::Quaternion, ϕc::AbstractVector, 
    timestep)

    # distance Jacobian 
    ∂d∂x = ∂distance∂x(relative, model.collision, xp, qp, xc, qc)
    ∂d∂q = ∂distance∂q(relative, model.collision, xp, qp, xc, qc)

    # contact point velocity Jacobian 
    ∂vt∂x = ∂relative_tangential_velocity∂x(relative, model, xp, qp, vp, ϕp, xc, qc, vc, ϕc)
    ∂vt∂q = ∂relative_tangential_velocity∂q(relative, model, xp, qp, vp, ϕp, xc, qc, vc, ϕc)

    X = [
            ∂d∂x;
            szeros(1, 3);
            model.friction_parameterization * ∂vt∂x;
        ]

    Q = [
            ∂d∂q;
            szeros(1, 4);
            model.friction_parameterization * ∂vt∂q
        ]

    return [X Q]
end

function constraint_jacobian_velocity(relative::Symbol, model::Contact, 
    xp::AbstractVector, vp::AbstractVector, qp::Quaternion, ϕp::AbstractVector, 
    xc::AbstractVector, vc::AbstractVector, qc::Quaternion, ϕc::AbstractVector, 
    timestep)

    # distance Jacobian
    ∂d∂x = ∂distance∂x(relative, model.collision, xp, qp, xc, qc)
    ∂d∂q = ∂distance∂q(relative, model.collision, xp, qp, xc, qc)

    # relative tangential velocity Jacobian 
    ∂vt∂q = ∂relative_tangential_velocity∂q(relative, model, xp, qp, vp, ϕp, xc, qc, vc, ϕc)
    ∂vt∂v = ∂relative_tangential_velocity∂v(relative, model, xp, qp, vp, ϕp, xc, qc, vc, ϕc)
    ∂vt∂ϕ = ∂relative_tangential_velocity∂ϕ(relative, model, xp, qp, vp, ϕp, xc, qc, vc, ϕc)

    # recover current orientation 
    if relative == :parent 
        x = next_position(xp, -vp, timestep)
        ∂x∂v = linear_integrator_jacobian_velocity(x, vp, timestep)
        q = next_orientation(qp, -ϕp, timestep)
        ∂q∂ϕ = rotational_integrator_jacobian_velocity(q, ϕp, timestep);
    elseif relative == :child 
        x = next_position(xc, -vc, timestep)
        ∂x∂v = linear_integrator_jacobian_velocity(x, vc, timestep)
        q = next_orientation(qc, -ϕc, timestep)
        ∂q∂ϕ = rotational_integrator_jacobian_velocity(q, ϕc, timestep);
    end
    
    V = [
            ∂d∂x * ∂x∂v;
            szeros(1, 3);
            model.friction_parameterization * ∂vt∂v
        ]
   
    Ω = [
            ∂d∂q * ∂q∂ϕ;
            szeros(1, 3);
            model.friction_parameterization * (∂vt∂ϕ + ∂vt∂q * ∂q∂ϕ);
        ]

    return [V Ω]
end

function impulse_map(relative::Symbol, model::Contact, pbody::Node, cbody::Node, timestep)
    # configurations
    xp, qp = next_configuration(pbody.state, timestep)
    xc, qc = next_configuration(cbody.state, timestep)

    # mapping
    X = force_mapping(relative, model, xp, qp, xc, qc)
    c = contact_point(relative, model.collision, xp, qp, xc, qc)

    if relative == :parent 
        r = c - xp
        Q = rotation_matrix(inv(qp)) * skew(r) * X
    elseif relative == :child 
        r = c - xc
        Q = rotation_matrix(inv(qc)) * skew(r) * X
    end

    return [
                X; 
                Q;
           ]
end

function impulse_map_jacobian(relative::Symbol, jacobian::Symbol, model::Contact, pbody::Node, cbody::Node, λ, timestep)

    # configurations
    xp, qp = next_configuration(pbody.state, timestep)
    xc, qc = next_configuration(cbody.state, timestep)

    # mapping
    X = force_mapping(relative, model, xp, qp, xc, qc)

    # force Jacobian 
    Xx = ∂force_mapping_jvp∂x(relative, jacobian, model, xp, qp, xc, qc, λ)
    Xq = ∂force_mapping_jvp∂q(relative, jacobian, model, xp, qp, xc, qc, λ)

    # contact point
    c = contact_point(relative, model.collision, xp, qp, xc, qc)

    # torque Jacobian
    if relative == :parent 
        r = c - xp
        q = qp
    elseif relative == :child 
        r = c - xc
        q = qc
    end

    Qx = rotation_matrix(inv(q)) * skew(r) * Xx 
    Qx -= rotation_matrix(inv(q)) * skew(X * λ) * (∂contact_point∂x(relative, jacobian, model.collision, xp, qp, xc, qc) - (relative == jacobian ? 1.0 : 0.0) * I(3)) 
   
    Qq = rotation_matrix(inv(q)) * skew(r) * Xq 
    Qq -= rotation_matrix(inv(q)) * skew(X * λ) * ∂contact_point∂q(relative, jacobian, model.collision, xp, qp, xc, qc) 
    Qq += ∂rotation_matrix∂q(inv(q), skew(r) * X * λ) * Tmat()

    return [
                Xx Xq; 
                Qx Qq;
           ]
end

# force mapping 
function force_mapping(relative::Symbol, model::Contact, 
    xp::AbstractVector, qp::Quaternion, 
    xc::AbstractVector, qc::Quaternion)

    X = [
        contact_normal(model.collision, xp, qp, xc, qc)' szeros(3, 1) contact_tangent(model.collision, xp, qp, xc, qc)' * model.friction_parameterization'
    ]

    if relative == :parent 
        return X
    elseif relative == :child 
        return -1.0 * X
    end
end

# mapping * λ
function ∂force_mapping_jvp∂x(relative::Symbol, jacobian::Symbol,
    model::Contact, 
    xp::AbstractVector, qp::Quaternion, 
    xc::AbstractVector, qc::Quaternion,
    λ::AbstractVector)

    # normal
    X = ∂contact_normal_transpose∂x(jacobian, model.collision, xp, qp, xc, qc) * λ[1]

    # tangent
    λ_tangent = model.friction_parameterization' * λ[SUnitRange(3, length(λ))]
    X += ∂contact_tangent_one_transpose∂x(jacobian, model.collision, xp, qp, xc, qc) * λ_tangent[1]
    X += ∂contact_tangent_two_transpose∂x(jacobian, model.collision, xp, qp, xc, qc) * λ_tangent[2]

    if relative == :parent 
        return X 
    elseif relative == :child 
        return -1.0 * X 
    end
end

# mapping * λ
function ∂force_mapping_jvp∂q(relative::Symbol, jacobian::Symbol,
    model::Contact, 
    xp::AbstractVector, qp::Quaternion, 
    xc::AbstractVector, qc::Quaternion,
    λ::AbstractVector)

    # normal
    X = ∂contact_normal_transpose∂q(jacobian, model.collision, xp, qp, xc, qc) * λ[1]

    # tangent
    λ_tangent = model.friction_parameterization' * λ[SUnitRange(3, length(λ))]
    X += ∂contact_tangent_one_transpose∂q(jacobian, model.collision, xp, qp, xc, qc) * λ_tangent[1]
    X += ∂contact_tangent_two_transpose∂q(jacobian, model.collision, xp, qp, xc, qc) * λ_tangent[2]

    if relative == :parent 
        return X 
    elseif relative == :child 
        return -1.0 * X 
    end
end

# utilities
Base.length(model::Contact{T,N}) where {T,N} = N
neutral_vector(model::Contact{T,N}) where {T,N} = sones(T, Int(N / 2))
cone_degree(model::Contact{T,N}) where {T,N} = Int(N / 2)




