"""
    Contact{T,N} 

    Abstract type containing contact information associated with Body objects.
"""
abstract type Contact{T,N} end

# constraint Jacobians
function constraint_jacobian_configuration(model::Contact, 
    xp3::AbstractVector, vp25::AbstractVector, qp3::UnitQuaternion, ϕp25::AbstractVector, 
    xc3::AbstractVector, vc25::AbstractVector, qc3::UnitQuaternion, ϕc25::AbstractVector, 
    timestep)

    # distance Jacobian 
    ∂d∂x = ∂distance∂xp(model.collision, xp3, qp3, xc3, qc3)
    ∂d∂q = ∂distance∂qp(model.collision, xp3, qp3, xc3, qc3)

    # contact point velocity Jacobian 
    ∂v∂x3 = ∂contact_point_velocity∂x(model, xp3, qp3, vp25, ϕp25)
    ∂v∂q3 = ∂contact_point_velocity∂q(model, xp3, qp3, vp25, ϕp25)

    X = [
            ∂d∂x;
            szeros(1, 3);
            model.collision.surface_projector * ∂v∂x3;
        ]

    Q = [
            ∂d∂q;
            szeros(1, 4);
            model.collision.surface_projector * ∂v∂q3
        ]

    return [X Q]
end

function constraint_jacobian_velocity(model::Contact, 
    xp3::AbstractVector, vp25::AbstractVector, qp3::UnitQuaternion, ϕp25::AbstractVector, 
    xc3::AbstractVector, vc25::AbstractVector, qc3::UnitQuaternion, ϕc25::AbstractVector, 
    timestep)

    # distance Jacobian
    ∂d∂x = ∂distance∂xp(model.collision, xp3, qp3, xc3, qc3)
    ∂d∂q = ∂distance∂qp(model.collision, xp3, qp3, xc3, qc3)

    # contact point velocity Jacobian 
    ∂v∂q3 = ∂contact_point_velocity∂q(model, xp3, qp3, vp25, ϕp25)
    ∂v∂v25 = ∂contact_point_velocity∂v(model, xp3, qp3, vp25, ϕp25)
    ∂v∂ϕ25 = ∂contact_point_velocity∂ϕ(model, xp3, qp3, vp25, ϕp25)

    # recover current orientation 
    qp2 = next_orientation(qp3, -ϕp25, timestep)
    
    V = [
            ∂d∂x * timestep;
            szeros(1, 3);
            model.collision.surface_projector * ∂v∂v25
        ]
   
    Ω = [
            ∂d∂q * rotational_integrator_jacobian_velocity(qp2, ϕp25, timestep);
            szeros(1, 3);
            model.collision.surface_projector * (∂v∂ϕ25 + ∂v∂q3 * rotational_integrator_jacobian_velocity(qp2, ϕp25, timestep))
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
    offset = model.collision.surface_normal_projector' * model.collision.contact_radius
    Q = - X * qp * skew(model.collision.contact_point - vector_rotate(offset, inv(qp))) 
    return [X'; Q']
end

# force mapping 
function force_mapping(model::Contact, xp::AbstractVector, qp::UnitQuaternion, xc::AbstractVector, qc::UnitQuaternion)
    X = [model.collision.surface_normal_projector;
         szeros(1,3);
         model.collision.surface_projector]
    return X
end

# utilities
Base.length(model::Contact{T,N}) where {T,N} = N
neutral_vector(model::Contact{T,N}) where {T,N} = sones(T, Int(N / 2))
cone_degree(model::Contact{T,N}) where {T,N} = Int(N / 2)




