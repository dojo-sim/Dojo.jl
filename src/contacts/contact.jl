"""
    Contact{T,N} 

    Abstract type containing contact information associated with Body objects.
"""
abstract type Contact{T,N} end

# constraint Jacobians
function constraint_jacobian_configuration(model::Contact, 
    xp::AbstractVector, vp::AbstractVector, qp::UnitQuaternion, ϕp::AbstractVector, 
    xc::AbstractVector, vc::AbstractVector, qc::UnitQuaternion, ϕc::AbstractVector, 
    timestep)

    relative = :parent 

    # distance Jacobian 
    ∂d∂x = ∂distance∂x(relative, model.collision, xp, qp, xc, qc)
    ∂d∂q = ∂distance∂q(relative, model.collision, xp, qp, xc, qc)

    # contact point velocity Jacobian 
    ∂vt∂x = ∂relative_tangential_velocity∂x(relative, model, xp, qp, vp, ϕp, xc, qc, vc, ϕc)#∂contact_point_velocity∂x(model, xp, qp, vp, ϕp)
    ∂vt∂q = ∂relative_tangential_velocity∂q(relative, model, xp, qp, vp, ϕp, xc, qc, vc, ϕc)#∂contact_point_velocity∂q(model, xp, qp, vp, ϕp)

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

function constraint_jacobian_velocity(model::Contact, 
    xp::AbstractVector, vp::AbstractVector, qp::UnitQuaternion, ϕp::AbstractVector, 
    xc::AbstractVector, vc::AbstractVector, qc::UnitQuaternion, ϕc::AbstractVector, 
    timestep)

    relative = :parent

    # distance Jacobian
    ∂d∂x = ∂distance∂x(relative, model.collision, xp, qp, xc, qc)
    ∂d∂q = ∂distance∂q(relative, model.collision, xp, qp, xc, qc)

    # contact point velocity Jacobian 
    ∂vt∂q = ∂relative_tangential_velocity∂q(relative, model, xp, qp, vp, ϕp, xc, qc, vc, ϕc)
    ∂vt∂v = ∂relative_tangential_velocity∂v(relative, model, xp, qp, vp, ϕp, xc, qc, vc, ϕc)
    ∂vt∂ϕ = ∂relative_tangential_velocity∂ϕ(relative, model, xp, qp, vp, ϕp, xc, qc, vc, ϕc)

    # recover current orientation 
    if relative == :parent 
        qp_prev = next_orientation(qp, -ϕp, timestep)
        ∂q_prev∂ϕ = rotational_integrator_jacobian_velocity(qp_prev, ϕp, timestep);
    elseif relative == :child 
        qp_prev = next_orientation(qc, -ϕc, timestep)
        ∂q_prev∂ϕ = rotational_integrator_jacobian_velocity(qp_prev, ϕc, timestep);
    end
    
    V = [
            ∂d∂x * timestep;
            szeros(1, 3);
            ∂vt∂v
        ]
   
    Ω = [
            ∂d∂q * ∂q_prev∂ϕ;
            szeros(1, 3);
            ∂vt∂ϕ + ∂vt∂q * ∂q_prev∂ϕ;
        ]

    return [V Ω]
end

# impulses
function impulse_map(model::Contact, pbody::Node, cbody::Node, timestep)
    impulse_map(model, 
        next_configuration(pbody.state, timestep)..., 
        next_configuration(cbody.state, timestep)...)
end

function impulse_map(model::Contact, xp::AbstractVector, qp::UnitQuaternion, xc::AbstractVector, qc::UnitQuaternion,)
    X = force_mapping(model, xp, qp, xc, qc)
    offset = model.collision.contact_normal' * model.collision.contact_radius
    Q = - X * qp * skew(model.collision.contact_origin - vector_rotate(offset, inv(qp))) 
    return [X'; Q']
end

# force mapping 
function force_mapping(model::Contact, xp::AbstractVector, qp::UnitQuaternion, xc::AbstractVector, qc::UnitQuaternion)
    X = [model.collision.contact_normal;
         szeros(1,3);
         model.collision.contact_tangent]
    return X
end

# utilities
Base.length(model::Contact{T,N}) where {T,N} = N
neutral_vector(model::Contact{T,N}) where {T,N} = sones(T, Int(N / 2))
cone_degree(model::Contact{T,N}) where {T,N} = Int(N / 2)




