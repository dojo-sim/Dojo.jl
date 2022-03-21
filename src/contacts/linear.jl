"""
    LinearContact{T,N} <: Contact{T,N}

    contact object for impact and friction with a linearized friction cone

    friction_coefficient: value of friction coefficient
    collision: Collision
"""
mutable struct LinearContact{T,N} <: Contact{T,N}
    friction_coefficient::T
    friction_parameterization::SMatrix{4,2,T,8}
    collision::Collision{T,4,3,12}
end

function LinearContact(body::Body{T}, normal::AbstractVector, friction_coefficient; 
    contact_origin=szeros(T, 3), 
    contact_radius=0.0) where T

    # contact directions
    V1, V2, V3 = orthogonal_columns(normal)
    A = [V1 V2 V3]
    Ainv = inv(A)
    contact_normal = Ainv[3, SA[1; 2; 3]]'
    contact_tangent = Ainv[SA[1; 2], SA[1; 2; 3]]

    # friction parameterization
    parameterization = SA{T}[
         0.0  1.0
         0.0 -1.0
         1.0  0.0
        -1.0  0.0
    ]

    # collision 
    collision = SphereFlatCollision(parameterization * contact_tangent, contact_normal, SVector{3}(contact_origin), contact_radius)
    
    LinearContact{Float64,12}(friction_coefficient, parameterization, collision)
end

function constraint_jacobian(contact::ContactConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs<:LinearContact{T,N},N½}
    friction_coefficient = contact.model.friction_coefficient
    γ = contact.impulses[2] + REG * neutral_vector(contact.model)
    s = contact.impulses_dual[2] + REG * neutral_vector(contact.model)

    ∇s1 = Diagonal(γ) # 6x6
    ∇s2 = Diagonal(-sones(T, 6))
    ∇s = vcat(∇s1, ∇s2) # 12x6

    ∇γ1 = Diagonal(s) # 6x6
    ∇γ2 = @SMatrix[0.0                  0.0  0.0  0.0  0.0  0.0;
                   friction_coefficient 0.0 -1.0 -1.0 -1.0 -1.0;
                   0.0                  1.0  0.0  0.0  0.0  0.0;
                   0.0                  1.0  0.0  0.0  0.0  0.0;
                   0.0                  1.0  0.0  0.0  0.0  0.0;
                   0.0                  1.0  0.0  0.0  0.0  0.0;]

    ∇γ = vcat(∇γ1, ∇γ2) # 12x6

    return [∇s ∇γ]
end

function constraint(mechanism, contact::ContactConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs<:LinearContact{T,N},N½}
    # contact model 
    model = contact.model

    # parent
    pbody = get_body(mechanism, contact.parent_id)
    xp, vp, qp, ϕp = next_configuration_velocity(pbody.state, mechanism.timestep)

    # child 
    cbody = get_body(mechanism, contact.child_id)
    xc, vc, qc, ϕc = next_configuration_velocity(cbody.state, mechanism.timestep)

    # distance 
    d = distance(model.collision, xp, qp, xc, qc)

    # relative tangential velocity
    vt = relative_tangential_velocity(model, xp, qp, vp, ϕp, xc, qc, vc, ϕc)

    # unpack contact variables
    γ = contact.impulses[2][1]
    sγ = contact.impulses_dual[2][1]
    ψ = contact.impulses[2][2]
    sψ = contact.impulses_dual[2][2]
    β = contact.impulses[2][@SVector [3,4,5,6]]
    sβ = contact.impulses_dual[2][@SVector [3,4,5,6]]

    SVector{N½,T}(
        d - sγ,
        model.friction_coefficient * γ - sum(β) - sψ,
        (model.friction_parameterization * vt + ψ * sones(4) - sβ)...)
end

