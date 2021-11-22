"""
    setPosition!(body; x, q)

Set the position and orientation of a body.
"""
function setPosition!(body::Body; x::AbstractVector = SA[0;0;0], q::UnitQuaternion = one(UnitQuaternion))
    body.state.x2[1] = x
    body.state.q2[1] = q
    return
end

"""
    setPosition!(body1, body2; p1, p2, Δx, Δq)

Set the position and orientation of body2 relative to body1 at the connection points p1 and p2.
"""
function setPosition!(body1::Body, body2::Body;
        p1::AbstractVector = SA[0;0;0], p2::AbstractVector = SA[0;0;0],
        Δx::AbstractVector = SA[0;0;0], Δq::UnitQuaternion = one(UnitQuaternion)
    )

    q1 = body1.state.q2[1]
    q2 = body1.state.q2[1] * Δq
    x2 = body1.state.x2[1] + vrotate(p1 + Δx, q1) - vrotate(p2, q2)
    setPosition!(body2;x = x2,q = q2)
    return
end

"""
    setPosition!(origin, body2; p1, p2, Δx, Δq)

Set the position and orientation of body2 relative to the origin at the connection points p1 and p2.
"""
function setPosition!(body1::Origin, body2::Body;
        p1::AbstractVector = SA[0;0;0], p2::AbstractVector = SA[0;0;0],
        Δx::AbstractVector = SA[0;0;0], Δq::UnitQuaternion = one(UnitQuaternion)
    )

    q2 = Δq
    x2 = p1 + Δx - vrotate(p2, q2)
    setPosition!(body2;x = x2,q = q2)
    return
end


"""
    setVelocity!(body; v, ω)

Set the translational and angular velocity of a body.
"""
function setVelocity!(body::Body; v::AbstractVector = SA[0;0;0], ω::AbstractVector = SA[0;0;0])
    body.state.v15 = v
    body.state.ϕ15 = ω
    return
end

"""
    setVelocity!(body1, body2; p1, p2 Δv, Δω)

Set the translational and angular velocity of body2 relative to body1 at the connection points p1 and p2.
"""
function setVelocity!(body1::Body, body2::Body;
        p1::AbstractVector = SA[0;0;0], p2::AbstractVector = SA[0;0;0],
        Δv::AbstractVector = SA[0;0;0], Δω::AbstractVector = SA[0;0;0]
    )

    q1 = body1.state.q2[1]
    q2 = body2.state.q2[1]

    v1 = body1.state.v15
    ω1 = body1.state.ϕ15 # in local coordinates

    vp1 = v1 + vrotate(cross(ω1,p1),q1)
    ωp1 = vrotate(ω1,q1) # in world coordinates

    vp2 = vp1 + vrotate(Δv,q1)
    ωp2 = ωp1 + vrotate(Δω,q2) # in world coordinates

    v2 = vp2 + vrotate(cross(vrotate(ωp2,inv(q2)),-p2),q2)
    ω2 = vrotate(ωp2,inv(q2)) # in local coordinates

    setVelocity!(body2;v = v2,ω = ω2)
    return
end

"""
    setVelocity!(origin, body2; p1, p2 Δv, Δω)

Set the translational and angular velocity of body2 relative to the origin at the connection points p1 and p2.
"""
function setVelocity!(body1::Origin, body2::Body;
        p1::AbstractVector = SA[0;0;0], p2::AbstractVector = SA[0;0;0],
        Δv::AbstractVector = SA[0;0;0], Δω::AbstractVector = SA[0;0;0]
    )

    q2 = body2.state.q2[1]

    vp2 = Δv
    ωp2 = vrotate(Δω,q2) # in world coordinates

    v2 = vp2 + cross(ωp2,-p2)
    ω2 = vrotate(ωp2,inv(q2)) # in local coordinates

    setVelocity!(body2;v = v2,ω = ω2)
    return
end


function setForce!(body::Body;
        F::AbstractVector = SA[0;0;0], τ::AbstractVector = SA[0;0;0], p::AbstractVector = SA[0;0;0]
    )
    # F and p in local coordinates
    τ += torqueFromForce(F, p) # in local coordinates
    setForce!(body.state, vrotate(F,body.state.q2[1]), τ)
    return
end

@inline torqueFromForce(F::AbstractVector, p::AbstractVector) = cross(p, F)
