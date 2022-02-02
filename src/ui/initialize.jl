function set_position!(body::Body; x::AbstractVector = SA[0;0;0], q::UnitQuaternion = one(UnitQuaternion))
    body.state.x2[1] = x
    body.state.q2[1] = q
    return
end

function set_position!(body1::Body, body2::Body;
        p1::AbstractVector = SA[0;0;0], p2::AbstractVector = SA[0;0;0],
        Δx::AbstractVector = SA[0;0;0], Δq::UnitQuaternion = one(UnitQuaternion)
        # in body1's frame
    )

    q1 = body1.state.q2[1]
    q2 = body1.state.q2[1] * Δq
    x2 = body1.state.x2[1] + vrotate(p1 + Δx, q1) - vrotate(p2, q2)
    set_position!(body2; x = x2, q = q2)
    return
end

# function set_position!(body1::Origin, body2::Body;
#         p1::AbstractVector = SA[0;0;0], p2::AbstractVector = SA[0;0;0],
#         Δx::AbstractVector = SA[0;0;0], Δq::UnitQuaternion = one(UnitQuaternion)
#     )

#     q2 = Δq
#     x2 = p1 + Δx - vrotate(p2, q2)
#     set_position!(body2;x = x2,q = q2)
#     return
# end

function set_velocity!(body::Body; v::AbstractVector = SA[0;0;0], ω::AbstractVector = SA[0;0;0])
    body.state.v15 = v
    body.state.ϕ15 = ω
    return
end

function set_velocity!(body1::Body, body2::Body;
        p1::AbstractVector = SA[0;0;0], p2::AbstractVector = SA[0;0;0],
        Δv::AbstractVector = SA[0;0;0], Δω::AbstractVector = SA[0;0;0]
        # in body1's frame              in body1's frame
    )

    x1 = body1.state.x2[1]
    v1 = body1.state.v15
    q1 = body1.state.q2[1]
    ω1 = body1.state.ϕ15 # in local coordinates

    x2 = body2.state.x2[1]
    v2 = body2.state.v15
    q2 = body2.state.q2[1]
    ω2 = body2.state.ϕ15 # in local coordinates

    # Ω(B/W)b = Ra->b * [Ω(B/A)a + Ω(A/W)a]
    ω2 = vrotate(Δω + ω1, inv(q2) * q1)
    # V(cb,B/W)w =
    ω1w = vrotate(ω1, q1)
    ω2w = vrotate(ω2, q2)
    Δvw = vrotate(Δv, q1)
    cApB_w = (x2 + vrotate(p2, q2)) - x1
    pBcB_w = - vrotate(p2, q2)
    v2 = v1
    v2 += skew(ω1w) * cApB_w
    v2 += skew(ω2w) * pBcB_w
    v2 += Δvw
    set_velocity!(body2;v = v2,ω = ω2)
    return
end

# function set_velocity!(body1::Origin, body2::Body;
#         p1::AbstractVector = SA[0;0;0], p2::AbstractVector = SA[0;0;0],
#         Δv::AbstractVector = SA[0;0;0], Δω::AbstractVector = SA[0;0;0]
#     )

#     x2 = body2.state.x2[1]
#     v2 = body2.state.v15
#     q2 = body2.state.q2[1]
#     ω2 = body2.state.ϕ15 # in local coordinates

#     # Ω(B/W)b = Ra->b * [Ω(B/A)a + Ω(A/W)a]
#     ω2 = vrotate(Δω, inv(q2))
#     # V(cb,B/W)w =
#     ω2w = vrotate(ω2, q2)
#     Δvw = Δv
#     pBcB_w = - vrotate(p2, q2)
#     v2 = Δvw + skew(ω2w) * pBcB_w
#     set_velocity!(body2; v = v2, ω = ω2)
#     return
# end

function set_input!(body::Body;
        F::AbstractVector = SA[0;0;0], τ::AbstractVector = SA[0;0;0], p::AbstractVector = SA[0;0;0]
    )
    # F and p in local coordinates
    τ += torque_from_force(F, p) # in local coordinates
    set_input!(body.state, vrotate(F,body.state.q2[1]), τ)
    return
end

@inline torque_from_force(F::AbstractVector, p::AbstractVector) = cross(p, F)
