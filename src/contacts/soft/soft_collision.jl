"""
    Soft Collision

    abstract type defining interaction between a soft body and another body
"""
abstract type SoftCollision{T,O,I,OI,N} end


"""
    ψ: overlap amount in mass unit.
    barycenter: barycenter of the overlapping particles expressed in the soft body's frame (i.e. parent frame)
    normal: normal of the contact point pointing towards the parent body, expressed in the world frame
"""
function overlap(collision::SoftCollision{T,O,I,OI,N}, xp, qp, xc, qc) where {T,O,I,OI,N}
    collider = collision.collider
	num_active = 0
	ψ = 0.0
	barycenter = szeros(T,3)

    for i = 1:N
        particle = collider.particles[i]
        p = xp + vector_rotate(particle + collision.collider_origin, qp)
        if inside(collision, p, xc, qc)
			num_active += 1
			wi = collider.weights[i]
			ψ += wi
			barycenter += wi * particle
        end
    end
	(num_active > 0) && (barycenter /= ψ)
    p = xp + vector_rotate(barycenter + collision.collider_origin, qp)
    normal = contact_normal(collision, p, xc, qc)
    return ψ, barycenter, normal
end


# xp = 0.05*sones(3)
# qp = Quaternion(1,1,0,0.0)
# xc = 0.03*sones(3)
# qc = Quaternion(1,0,2,0.0)
# collision = mech.contacts[1].model.collision
# ψ, barycenter, normal = overlap(collision, xp, qp, xc, qc)
# ψ2, barycenter2, normal2 = overlap2(collision, xp, qp, xc, qc)
# # @benchmark $overlap($collision, $xp, $qp, $xc, $qc)
# ψ - ψ2
# norm(barycenter - barycenter2)
# norm(normal - normal2)
