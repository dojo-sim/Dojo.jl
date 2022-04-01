"""
    Soft Collision

    abstract type defining interaction between a soft body and another body
"""
abstract type SoftCollision{T,O,I,OI,N} end

function overlap(collision::SoftCollision{T,O,I,OI,N}, xp, qp, xc, qc) where {T,O,I,OI,N}
    collider = collision.collider
    center_of_mass = collider.center_of_mass
    Ψ = Vector{T}()
    active_particles = []
    barycenter = szeros(T,3)

    for i = 1:N
        particle = collider.particles[i]
        p = xp + Dojo.vector_rotate(particle - center_of_mass, qp)
        if inside(collision, p, xc, qc)
            push!(Ψ, collider.weights[i])
            push!(active_particles, particle)
        end
    end

    num_active = length(Ψ)
    ψ = sum(Ψ)
    (ψ > 0) && (barycenter = sum(Vector{SVector{3,T}}([Ψ[i] * active_particles[i] for i=1:num_active])) / ψ)

    p = xp + Dojo.vector_rotate(barycenter - center_of_mass, qp)
    normal = contact_normal(collision, p, xc, qc)
    return ψ, barycenter, normal
end
