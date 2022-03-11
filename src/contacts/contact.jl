"""
    Contact{T,N} 

    Abstract type containing contact information associated with Body objects.
"""
abstract type Contact{T,N} end

# constraint Jacobians
function constraint_jacobian_configuration(relative::Symbol, model::Contact, 
    xp::AbstractVector, vp::AbstractVector, qp::UnitQuaternion, ϕp::AbstractVector, 
    xc::AbstractVector, vc::AbstractVector, qc::UnitQuaternion, ϕc::AbstractVector, 
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
            ∂vt∂x;
        ]

    Q = [
            ∂d∂q;
            szeros(1, 4);
            ∂vt∂q
        ]

    return [X Q]
end

function constraint_jacobian_velocity(relative::Symbol, model::Contact, 
    xp::AbstractVector, vp::AbstractVector, qp::UnitQuaternion, ϕp::AbstractVector, 
    xc::AbstractVector, vc::AbstractVector, qc::UnitQuaternion, ϕc::AbstractVector, 
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
        q = next_orientation(qp, -ϕp, timestep)
        ∂q∂ϕ = rotational_integrator_jacobian_velocity(q, ϕp, timestep);
    elseif relative == :child 
        q = next_orientation(qc, -ϕc, timestep)
        ∂q∂ϕ = rotational_integrator_jacobian_velocity(q, ϕc, timestep);
    end
    
    V = [
            ∂d∂x * timestep;
            szeros(1, 3);
            ∂vt∂v
        ]
   
    Ω = [
            ∂d∂q * ∂q∂ϕ;
            szeros(1, 3);
            ∂vt∂ϕ + ∂vt∂q * ∂q∂ϕ;
        ]

    return [V Ω]
end

function impulse_map(relative::Symbol, model::Contact, pbody::Node, cbody::Node, timestep)
    # configurations
    xp, qp = next_configuration(pbody.state, timestep)
    xc, qc = next_configuration(cbody.state, timestep)

    # mapping
    X = force_mapping(relative, model, xp, qp, xc, qc)
    offset = model.collision.contact_normal' * model.collision.contact_radius
    Q = - X * qp * skew(model.collision.contact_origin - vector_rotate(offset, inv(qp))) 

    return [X'; Q']
end

# force mapping 
function force_mapping(relative::Symbol, model::Contact, 
    xp::AbstractVector, qp::UnitQuaternion, 
    xc::AbstractVector, qc::UnitQuaternion)

    X = [
            model.collision.contact_normal;
            szeros(1, 3);
            model.collision.contact_tangent
        ]

    return X
end

# utilities
Base.length(model::Contact{T,N}) where {T,N} = N
neutral_vector(model::Contact{T,N}) where {T,N} = sones(T, Int(N / 2))
cone_degree(model::Contact{T,N}) where {T,N} = Int(N / 2)




