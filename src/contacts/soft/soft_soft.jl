"""
    SoftSoftCollision

    collision between a soft contact and another soft contact

    collider: soft parent object
    child_collider: soft child object
    collider_origin: position of contact on parent collider relative to its center of mass
    child_collider_origin: position of contact on child collider relative to its center of mass
    contact_tangent: mapping from world frame to surface tangent frame
"""
mutable struct SoftSoftCollision{T,O,I,OI,N,Nc} <: SoftCollision{T,O,I,OI,N}
    collider::SoftCollider{T,N}
    child_collider::SoftCollider{T,Nc}
    collider_origin::SVector{I,T}
    child_collider_origin::SVector{I,T}
    contact_tangent::SMatrix{O,I,T,OI}
end

function inside(collision::SoftSoftCollision{T}, p, xc, qc) where T
    # p = barycenter of the overlap expressed in the world frame
    # pc = center of the sphere expressed in the world frame
    pc = xc + Dojo.vector_rotate(collision.sphere_origin, qc)
    c = norm(p - pc)
    return c <= collision.contact_radius
end

# normal projection (from child to parent)
function contact_normal(collision::SoftSoftCollision, p, xc, qc)
    # p = barycenter of the overlap expressed in the world frame
    # pc = center of the sphere expressed in the world frame
    pc = xc + Dojo.vector_rotate(collision.sphere_origin, qc)
    normal = p - pc
    return normal ./ (1e-20 + norm(normal))
end

parent_origin(collision::SoftSoftCollision) = collision.collider_origin
child_origin(collision::SoftSoftCollision) = collision.child_collider_origin

"""
    ψ: overlap amount in mass unit.
    barycenter: barycenter of the overlapping particles expressed in the soft body's frame (i.e. parent frame)
    normal: normal of the contact point pointing towards the parent body, expressed in the world frame
"""
function overlap(collision::SoftSoftCollision{T,O,I,OI,Np,Nc}, xp, qp, xc, qc) where {T,O,I,OI,Np,Nc}
    colliderp = collision.collider
    colliderc = collision.child_collider
    center_of_massp = colliderp.center_of_mass
    center_of_massc = colliderc.center_of_mass

    ψp, barycenterpw, normal_pw = cross_overlap(colliderp, colliderc, xp, qp, xc, qc) # weigths of collider2 at particle1
    ψc, barycentercw, normal_cw = cross_overlap(colliderc, colliderp, xp, qp, xc, qc) # weigths of collider1 at particle2
    ψ = ψp + ψc
    normal_w = (ψp * normal_pw + ψc * -normal_cw)
    normal_w /= (1e-20 + norm(normal_w))
    barycenter_w = (ψp * barycenterpw + ψc * barycentercw) / (1e-20 + ψ)
    barycenter_p = Dojo.vector_rotate(barycenter_w - xp, inv(qp))
    return ψ, barycenter_p, normal_w
end

# function cross_overlap(collider1::SoftCollider{T,N1}, collider2::SoftCollider{T,N2}, x1, q1, x2, q2) where {T,N1,N2}
#     # returns
#     # the particles of collider 1 in world frame
#     # sum of cross weights of collider1 and collider2
#     # barycenter of contact of collider1 into collider2
#     # contact normal from collider1 to collider2
#     center_of_mass1 = collider1.center_of_mass
#     center_of_mass2 = collider2.center_of_mass
#
#     particles_w = zeros(Float32,N1,3)
#     particles_2 = zeros(Float32,N1,3)
#     for i = 1:N1
#         particle1 = collider1.particles[i]
#         pw = x1 + vector_rotate(particle1 - center_of_mass1, q1)
#         p2 = center_of_mass2 + vector_rotate(pw - x2, inv(q2))
#         particles_w[i,:] = pw
#         particles_2[i,:] = p2
#     end
#     densities = py"density_query"(collider2.nerf_object, convert.(Float32, particles_2))
#     weights = densities ./ sum(collider2.densities) * collider2.mass
#
#     ψ = 0.0
#     barycenter_w = zeros(3)
#     normal_w = zeros(3)
#     for i = 1:N1
#         weight_product = collider1.weights[i] * weights[i]
#         ψ += weight_product
#         barycenter_w += weight_product * particles_w[i,:]
#         normal_w += weight_product * Dojo.vector_rotate(-collider1.weight_gradients[i], q1)
#     end
#     barycenter_w /= (1e-20 + ψ)
#     normal_w /= (1e-20 + norm(normal_w))
#     return ψ, barycenter_w, normal_w
# end
