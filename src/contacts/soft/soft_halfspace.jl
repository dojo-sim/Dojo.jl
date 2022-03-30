"""
    SoftHalfSpaceCollision

    collision between a soft contact and a flat surface

    contact_tangent: mapping from world frame to surface tangent frame
    contact_normal: inverse/complement of contact_tangent
    contact_origin: position of contact on Body relative to center of mass
    contact_radius: radius of contact
"""
mutable struct SoftHalfSpaceCollision{T,O,I,OI,N} <: SoftCollision{T,O,I,OI}
    collider::SoftCollider{T,N}
    contact_tangent::SMatrix{O,I,T,OI}
    contact_normal::Adjoint{T,SVector{I,T}}
    contact_origin::SVector{I,T}
end

function inside(collision::SoftHalfSpaceCollision{T}, p) where T
    origin = szeros(T,3)
    normal = collision.contact_normal
    c = normal * (p - origin)
    return c <= 0.0
end

function overlap(collision::SoftHalfSpaceCollision{T,O,I,OI,N}, xp, qp, xc, qc) where {T,O,I,OI,N}
    collider = collision.collider
    center_of_mass = collider.center_of_mass
    Ψ = Vector{T}()
    active_particles = []
    barycenter = szeros(T,3)

    for i = 1:N
        particle = collider.particles[i]
        p = xp + Dojo.vector_rotate(particle - center_of_mass, qp)
        if inside(collision, p)
            push!(Ψ, collider.weights[i])
            push!(active_particles, particle)
        end
    end

    num_active = length(Ψ)
    ψ = sum(Ψ)
    (ψ > 0) && (barycenter = sum(Vector{SVector{3,T}}([Ψ[i] * active_particles[i] for i=1:num_active])) / ψ)
    normal = contact_normal(collision, barycenter)
    return ψ, barycenter, normal
end

# normal projection (from child to parent)
function contact_normal(collision::SoftHalfSpaceCollision, p)
    return collision.contact_normal[1,:]
end

function constraint(model::SoftContact{T},
        xp::AbstractVector, vp::AbstractVector, qp::Quaternion, ϕp::AbstractVector,
        xc::AbstractVector, vc::AbstractVector, qc::Quaternion, ϕc::AbstractVector,
        timestep; intermediate::Bool=false) where T

    collision = model.collision
    collider = collision.collider
    opts = collider.options
    smoothing = opts.coulomb_smoothing
    regularizer = opts.coulomb_regularizer
    x2p = next_position(xp, -vp, timestep)
    q2p = next_orientation(qp, -ϕp, timestep)
    x2c = next_position(xc, -vc, timestep)
    q2c = next_orientation(qc, -ϕc, timestep)

    # collision
    ψ, barycenter, normal = overlap(collision, x2p, q2p, x2c, q2c)
    # velocities
    vc = vp + vector_rotate(skew(collider.center_of_mass - barycenter) * ϕp, q2p) # TODO could be q3p
    vc_normal = normal * normal' * vc
    vc_tangential = vc - vc_normal
    # impact
    F_impact = ψ * opts.impact_spring * normal
    F_impact -= ψ * opts.impact_damper * vc_normal
    # friction
    F_friction = -opts.sliding_drag * norm(F_impact) * vc_tangential
    F_friction += opts.sliding_friction * norm(F_impact) * coulomb_direction(vc_tangential, smoothing, regularizer)
    F_contact = F_impact + F_friction

    τ_contact = -opts.rolling_drag * norm(F_impact) * ϕp
    τ_contact += opts.rolling_friction * norm(F_impact) * coulomb_direction(ϕp, smoothing, regularizer)

    # constraint
    impulse = timestep * [F_contact; τ_contact]
    !intermediate && return impulse
    return ψ, barycenter, normal, vc_normal, vc_tangential, F_impact, F_friction, impulse
end

function constraint_jacobian_velocity(relative::Symbol, model::SoftContact{T},
    xp::AbstractVector, vp::AbstractVector, qp::Quaternion, ϕp::AbstractVector,
    xc::AbstractVector, vc::AbstractVector, qc::Quaternion, ϕc::AbstractVector,
    timestep) where T

    collision = model.collision
    collider = collision.collider
    opts = collider.options
    smoothing = opts.coulomb_smoothing
    regularizer = opts.coulomb_regularizer
    ψ, barycenter, contact_normal, vc_normal, vc_tangential, F_impact, F_friction, impulse = constraint(model,
        xp, vp, qp, ϕp,
        xc, vc, qc, ϕc,
        timestep, intermediate=true)
    # recover current orientation
    q2p = next_orientation(qp, -ϕp, timestep)
    ∂vc∂v = Diagonal(sones(3))
    ∂vc∂ϕ = rotation_matrix(q2p) * skew(collider.center_of_mass - barycenter)
    if relative == :parent
        # impact
        ∂F_impact∂vc = -ψ * opts.impact_damper * contact_normal * contact_normal'
        FV_impact = ∂F_impact∂vc * ∂vc∂v
        FΩ_impact = ∂F_impact∂vc * ∂vc∂ϕ
        # friction
        ∂F_friction∂vc =  -opts.sliding_drag * vc_tangential * (FV_impact * F_impact/(norm(F_impact) + 1e-20))'
        ∂F_friction∂vc += -opts.sliding_drag * norm(F_impact) * (Diagonal(sones(3)) - contact_normal * contact_normal')
        ∂F_friction∂vc +=  opts.sliding_friction *
            coulomb_direction(vc_tangential, smoothing, regularizer) *
            (FV_impact * F_impact/(norm(F_impact) + 1e-20))'
        ∂F_friction∂vc +=  opts.sliding_friction * norm(F_impact) *
            ∂coulomb_direction∂v(vc_tangential, smoothing, regularizer) *
            (Diagonal(sones(3)) - contact_normal * contact_normal')
        FV_friction = ∂F_friction∂vc * ∂vc∂v
        FΩ_friction = ∂F_friction∂vc * ∂vc∂ϕ

        τV_contact = (-opts.rolling_drag * ϕp + opts.rolling_friction * coulomb_direction(ϕp, smoothing, regularizer)) *
            (FV_impact * F_impact/(norm(F_impact) + 1e-20))'
        τΩ_contact = τV_contact * ∂vc∂ϕ
        τΩ_contact += norm(F_impact) * (-opts.rolling_drag * Diagonal(sones(3)) +
            opts.rolling_friction * ∂coulomb_direction∂v(ϕp, smoothing, regularizer))
        # contact
        V = [FV_impact + FV_friction; τV_contact]
        Ω = [FΩ_impact + FΩ_friction; τΩ_contact]
    elseif relative == :child
        V = szeros(T,N,3)
        Ω = szeros(T,N,3)
    end
    return timestep * [V Ω]
end


#
#
# # distance
# function distance(collision::SoftHalfSpaceCollision, xp, qp, xc, qc)
#     collision.contact_normal * (xp + vector_rotate(collision.contact_origin, qp)) - collision.contact_radius
# end
#
# function ∂distance∂x(gradient::Symbol, collision::SoftHalfSpaceCollision, xp, qp, xc, qc)
#     if gradient == :parent
#         return collision.contact_normal
#     elseif gradient == :child
#         return szeros(eltype(xc), 1, 3)
#     end
# end
#
# function ∂distance∂q(gradient::Symbol, collision::SoftHalfSpaceCollision, xp, qp, xc, qc)
#     if gradient == :parent
#         return collision.contact_normal * ∂vector_rotate∂q(collision.contact_origin, qp)
#     elseif gradient == :child
#         return szeros(eltype(qc), 1, 4)
#     end
# end
#
# # contact point in world frame
# function contact_point(relative::Symbol, collision::SoftHalfSpaceCollision, xp, qp, xc, qc)
#     if relative == :parent
#         return xp + vector_rotate(collision.contact_origin, qp) - collision.contact_normal' * collision.contact_radius
#     elseif relative == :child
#         projector = collision.contact_tangent' * collision.contact_tangent
#         return projector * (xp + vector_rotate(collision.contact_origin, qp))
#     end
# end
#
# function ∂contact_point∂x(relative::Symbol, jacobian::Symbol, collision::SoftHalfSpaceCollision, xp, qp, xc, qc)
#     if relative == :parent
#         if jacobian == :parent
#             return 1.0 * I(3)
#         elseif jacobian == :child
#             return szeros(eltype(xp), 3, 3)
#         end
#     elseif relative == :child
#         if jacobian == :parent
#             projector = collision.contact_tangent' * collision.contact_tangent
#             return projector
#         elseif jacobian == :child
#             return szeros(eltype(xp), 3, 3)
#         end
#     end
# end
#
# function ∂contact_point∂q(relative::Symbol, jacobian::Symbol, collision::SoftHalfSpaceCollision, xp, qp, xc, qc)
#     if relative == :parent
#         if jacobian == :parent
#             return ∂vector_rotate∂q(collision.contact_origin, qp)
#         end
#     elseif relative == :child
#         if jacobian == :parent
#             projector = collision.contact_tangent' * collision.contact_tangent
#             return projector * ∂vector_rotate∂q(collision.contact_origin, qp)
#         elseif jacobian == :child
#             return szeros(eltype(qp), 3, 4)
#         end
#     end
# end
#
# # normal projection (from child to parent)
# function contact_normal(collision::SoftHalfSpaceCollision, xp, qp, xc, qc)
#     return collision.contact_normal
# end
#
# function ∂contact_normal_transpose∂x(jacobian::Symbol, collision::SoftHalfSpaceCollision, xp, qp, xc, qc)
#     return szeros(eltype(collision.contact_normal), 3, 3)
# end
#
# function ∂contact_normal_transpose∂q(jacobian::Symbol, collision::SoftHalfSpaceCollision, xp, qp, xc, qc)
#     return szeros(eltype(collision.contact_normal), 3, 4)
# end
#
# # tangent projection
# function contact_tangent(collision::SoftHalfSpaceCollision, xp, qp, xc, qc)
#     return collision.contact_tangent # {2,4} x 3
# end
#
# function ∂contact_tangent_one_transpose∂x(jacobian::Symbol, collision::SoftHalfSpaceCollision, xp, qp, xc, qc)
#     return szeros(eltype(collision.contact_tangent), 3, 3)
# end
#
# function ∂contact_tangent_two_transpose∂x(jacobian::Symbol, collision::SoftHalfSpaceCollision, xp, qp, xc, qc)
#     return szeros(eltype(collision.contact_tangent), 3, 3)
# end
#
# function ∂contact_tangent_one_transpose∂q(jacobian::Symbol, collision::SoftHalfSpaceCollision, xp, qp, xc, qc)
#     return szeros(eltype(collision.contact_tangent), 3, 4)
# end
#
# function ∂contact_tangent_two_transpose∂q(jacobian::Symbol, collision::SoftHalfSpaceCollision, xp, qp, xc, qc)
#     return szeros(eltype(collision.contact_tangent), 3, 4)
# end
