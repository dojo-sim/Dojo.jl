"""
    LinearContact{T,N} <: Contact{T,N}

    contact object for impact and friction with a linearized friction cone

    friction_coefficient: value of friction coefficient
    surface_projector: mapping from world frame to surface tangent frame 
    surface_normal_projector: inverse/complement of surface_projector
    contact_point: position of contact on Body relative to center of mass 
    offset: position of contact relative to contact_point
"""
mutable struct LinearContact{T,N} <: Contact{T,N}
    friction_coefficient::T
    surface_projector::SMatrix{4,3,T,12}
    surface_normal_projector::Adjoint{T,SVector{3,T}} # inverse matrix
    contact_point::SVector{3,T}
    offset::SVector{3,T}

    function LinearContact(body::Body{T}, normal::AbstractVector, friction_coefficient; 
        contact_point=szeros(T, 3), 
        offset::AbstractVector=szeros(T, 3)) where T
        V1, V2, V3 = orthogonal_columns(normal)
        A = [V1 V2 V3]
        Ainv = inv(A)
        surface_normal_projector = Ainv[3, SA[1; 2; 3]]'
        surface_projector = SA{T}[
             1.0  0.0  0.0
            -1.0  0.0  0.0
             0.0  1.0  0.0
             0.0 -1.0  0.0
        ]
        new{Float64,12}(friction_coefficient, surface_projector, surface_normal_projector, contact_point, offset)
    end
end

function constraint(mechanism, contact::ContactConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs<:LinearContact{T,N},N½}
    # contact model 
    model = contact.model

    # parent
    pbody = get_body(mechanism, contact.parent_id)
    xp3, vp25, qp3, ϕp25 = next_configuration_velocity(pbody.state, mechanism.timestep)

    # child 
    cbody = get_body(mechanism, contact.child_id)
    xc3, vc25, qc3, ϕc25 = next_configuration_velocity(cbody.state, mechanism.timestep)

    # distance 
    d = distance(model, xp3, qp3, xc3, qc3)

    # relative tangential velocity
    vt = relative_tangential_velocity(model, xp3, qp3, vp25, ϕp25, xc3, qc3, vc25, ϕc25)

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
        (vt + ψ * sones(4) - sβ)...)
end

function complementarity_jacobian(contact::ContactConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs<:LinearContact{T,N},N½}
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

