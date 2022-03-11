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

    # distance Jacobian 
    ∂d∂x = ∂distance∂xp(model.collision, xp, qp, xc, qc)
    ∂d∂q = ∂distance∂qp(model.collision, xp, qp, xc, qc)

    # contact point velocity Jacobian 
    ∂v∂x = ∂contact_point_velocity∂x(model, xp, qp, vp, ϕp)
    ∂v∂q = ∂contact_point_velocity∂q(model, xp, qp, vp, ϕp)

    X = [
            ∂d∂x;
            szeros(1, 3);
            model.collision.contact_tangent * ∂v∂x;
        ]

    Q = [
            ∂d∂q;
            szeros(1, 4);
            model.collision.contact_tangent * ∂v∂q
        ]

    return [X Q]
end

function constraint_jacobian_velocity(model::Contact, 
    xp::AbstractVector, vp::AbstractVector, qp::UnitQuaternion, ϕp::AbstractVector, 
    xc::AbstractVector, vc::AbstractVector, qc::UnitQuaternion, ϕc::AbstractVector, 
    timestep)

    # distance Jacobian
    ∂d∂x = ∂distance∂xp(model.collision, xp, qp, xc, qc)
    ∂d∂q = ∂distance∂qp(model.collision, xp, qp, xc, qc)

    # contact point velocity Jacobian 
    ∂v∂q = ∂contact_point_velocity∂q(model, xp, qp, vp, ϕp)
    ∂v∂v = ∂contact_point_velocity∂v(model, xp, qp, vp, ϕp)
    ∂v∂ϕ = ∂contact_point_velocity∂ϕ(model, xp, qp, vp, ϕp)

    # recover current orientation 
    qp_prev = next_orientation(qp, -ϕp, timestep)
    
    V = [
            ∂d∂x * timestep;
            szeros(1, 3);
            model.collision.contact_tangent * ∂v∂v
        ]
   
    Ω = [
            ∂d∂q * rotational_integrator_jacobian_velocity(qp_prev, ϕp, timestep);
            szeros(1, 3);
            model.collision.contact_tangent * (∂v∂ϕ + ∂v∂q * rotational_integrator_jacobian_velocity(qp_prev, ϕp, timestep))
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




