# t3r3
"""
    Fixed(body1::AbstractBody, body2; p1, p2, qoffset)

A fixed connection between two bodies.
"""
Fixed(body1::AbstractBody{T}, body2; p1 = szeros(T, 3), p2 = szeros(T, 3), qoffset = one(UnitQuaternion{T})) where T =
    Translational3{T}(body1, body2; p1, p2), Rotational3{T}(body1, body2; qoffset)

# t2r3
"""
    Prismatic(body1, body2, axis; p1, p2, qoffset, spring, damper)

A prismatic joint between two bodies.
"""
Prismatic(body1::AbstractBody{T}, body2, axis; p1 = szeros(T, 3), p2 = szeros(T, 3), qoffset = one(UnitQuaternion{T}), spring = zero(T), damper = zero(T)) where T =
    Translational2{T}(body1, body2; p1, p2, axis, spring, damper), Rotational3{T}(body1, body2; qoffset, spring, damper)

# t1r3
"""
    Planar(body1, body2, axis; p1, p2, qoffset, spring, damper)

A planar joint between two bodies.
"""
Planar(body1::AbstractBody{T}, body2, axis; p1 = szeros(T, 3), p2 = szeros(T, 3), qoffset = one(UnitQuaternion{T}), spring = zero(T), damper = zero(T)) where T =
    Translational1{T}(body1, body2; p1, p2, axis, spring, damper), Rotational3{T}(body1, body2; qoffset, spring, damper)

# t0r3
"""
    FixedOrientation(body1, body2; qoffset, spring, damper)

Fixed orientation between two bodies (chicken's head).
"""
FixedOrientation(body1::AbstractBody{T}, body2; qoffset = one(UnitQuaternion{T}), spring = zero(T), damper = zero(T)) where T =
    Translational0{T}(body1, body2; spring, damper), Rotational3{T}(body1, body2; qoffset, spring, damper)

# t3r2
"""
    Revolute(body1, body2, axis; p1, p2, qoffset, spring, damper)

A revolute joint between two bodies (pin, continuous, hinge joint).
"""
Revolute(body1::AbstractBody{T}, body2, axis; p1 = szeros(T, 3), p2 = szeros(T, 3), qoffset = one(UnitQuaternion{T}), spring = zero(T), damper = zero(T)) where T =
    Translational3{T}(body1, body2; p1, p2, spring, damper), Rotational2{T}(body1, body2; axis, qoffset, spring, damper)

# t2r2
"""
    Cylindrical(body1, body2, axis; p1, p2, qoffset, spring, damper)

A cylindrical joint between two bodies.
"""
Cylindrical(body1::AbstractBody{T}, body2, axis; p1 = szeros(T, 3), p2 = szeros(T, 3), qoffset = one(UnitQuaternion{T}), spring = zero(T), damper = zero(T)) where T =
    Translational2{T}(body1, body2; p1, p2, axis, spring, damper), Rotational2{T}(body1, body2; axis, qoffset, spring, damper)
    # Translational2{T}(body1, body2; p1 = p1, p2 = p2, axis = [0,0,1.0], spring = spring, damper = damper), Rotational2{T}(body1, body2; axis = SVector{3,T}(0,1.0,0), qoffset = qoffset, spring = spring, damper = damper)

# t1r2
"""
    PlanarAxis(body1, body2, axis; p1, p2, qoffset, spring, damper)

A planar joint between two bodies with a rotation axis perpendicular to the plane (turtle bot).
"""
PlanarAxis(body1::AbstractBody{T}, body2, axis; p1 = szeros(T, 3), p2 = szeros(T, 3), qoffset = one(UnitQuaternion{T}), spring = zero(T), damper = zero(T)) where T =
    Translational1{T}(body1, body2; p1, p2, axis, spring, damper), Rotational2{T}(body1, body2; axis, qoffset, spring, damper)

# t0r2
"""
    FreeRevolute(body1, body2, axis; p1, p2, qoffset, spring, damper)

A joint between two bodies with free translation and rotation along one axis.
"""
FreeRevolute(body1::AbstractBody{T}, body2, axis; p1 = szeros(T, 3), p2 = szeros(T, 3), qoffset = one(UnitQuaternion{T}), spring = zero(T), damper = zero(T)) where T =
    Translational0{T}(body1, body2; spring, damper), Rotational2{T}(body1, body2; axis, qoffset, spring, damper)

#t3r1
"""
    Orbital(body1, body2, axis; p1, p2, qoffset, spring, damper)

A rotational between two bodies with a 2 rotational degrees of freedom (skull-eye joint).
"""
Orbital(body1::AbstractBody{T}, body2, axis; p1 = szeros(T, 3), p2 = szeros(T, 3), qoffset = one(UnitQuaternion{T}), spring = zero(T), damper = zero(T)) where T =
    Translational3{T}(body1, body2; p1, p2, axis, spring, damper), Rotational1{T}(body1, body2; axis, qoffset, spring, damper)

#t2r1
"""
    PrismaticOrbital(body1, body2, axis; p1, p2, qoffset, spring, damper)

A prismatic joint between two bodies with a 2 rotational degrees of freedom (skull-eye joint).
"""
PrismaticOrbital(body1::AbstractBody{T}, body2, axis; p1 = szeros(T, 3), p2 = szeros(T, 3), qoffset = one(UnitQuaternion{T}), spring = zero(T), damper = zero(T)) where T =
    Translational2{T}(body1, body2; p1, p2, axis, spring, damper), Rotational1{T}(body1, body2; axis, qoffset, spring, damper)

#t1r1
"""
    PlanarOrbital(body1, body2, axis; p1, p2, qoffset, spring, damper)

A planar joint between two bodies with a 2 rotational degrees of freedom (skull-eye joint).
"""
PlanarOrbital(body1::AbstractBody{T}, body2, axis; p1 = szeros(T, 3), p2 = szeros(T, 3), qoffset = one(UnitQuaternion{T}), spring = zero(T), damper = zero(T)) where T =
    Translational1{T}(body1, body2; p1, p2, axis, spring, damper), Rotational1{T}(body1, body2; axis, qoffset, spring, damper)

#t0r1
"""
    FreeOrbital(body1, body2, axis; p1, p2, qoffset, spring, damper)

A free joint between two bodies with a 2 rotational degrees of freedom (skull-eye joint).
"""
FreeOrbital(body1::AbstractBody{T}, body2, axis; p1 = szeros(T, 3), p2 = szeros(T, 3), qoffset = one(UnitQuaternion{T}), spring = zero(T), damper = zero(T)) where T =
    Translational0{T}(body1, body2; spring, damper), Rotational1{T}(body1, body2; axis, qoffset, spring, damper)

#t3r0
"""
    Spherical(body1, body2; p1, p2, spring, damper)

A spherical joint between two bodies (ball-and-socket joint).
"""
Spherical(body1::AbstractBody{T}, body2; p1 = szeros(T, 3), p2 = szeros(T, 3), qoffset = one(UnitQuaternion{T}), spring = zero(T), damper = zero(T)) where T =
    Translational3{T}(body1, body2; p1, p2, spring, damper), Rotational0{T}(body1, body2; qoffset, spring, damper)

# t2r0
"""
    CylindricalFree(body1, body2, axis; p1, p2, spring, damper)

A cylindrical joint between two bodies with unconstrained orientation (point-on-line).
"""
CylindricalFree(body1::AbstractBody{T}, body2, axis; p1 = szeros(T, 3), p2 = szeros(T, 3), spring = zero(T), damper = zero(T)) where T =
    Translational2{T}(body1, body2; p1, p2, axis, spring, damper), Rotational0{T}(body1, body2; spring, damper)

# t1r0
"""
    PlanarFree(body1, body2, axis; p1, p2, spring, damper)

A planar joint between two bodies with unconstrained orientation.
"""
PlanarFree(body1::AbstractBody{T}, body2, axis; p1 = szeros(T, 3), p2 = szeros(T, 3), spring = zero(T), damper = zero(T)) where T =
    Translational1{T}(body1, body2; p1, p2, axis, spring, damper), Rotational0{T}(body1, body2; spring, damper)

# t0r0
"""
    Floating(body1, body2; spring, damper)

An unconstrained connection between two bodies (connection between floating base and origin).
"""
Floating(body1::AbstractBody{T}, body2; spring = zero(T), damper = zero(T)) where T =
    Translational0{T}(body1, body2; spring, damper), Rotational0{T}(body1, body2; spring, damper)



function Prototype(jointtype::Symbol, body1::AbstractBody{T}, body2, axis; p1 = szeros(T, 3), p2 = szeros(T, 3),
        qoffset = one(UnitQuaternion{T}), spring = zero(T), damper = zero(T)) where T
    (jointtype == :Fixed)            && (return            Fixed(body1, body2;       p1 = p1, p2 = p2, qoffset = qoffset))
    (jointtype == :Prismatic)        && (return        Prismatic(body1, body2, axis; p1 = p1, p2 = p2, qoffset = qoffset, spring = spring, damper = damper))
    (jointtype == :Planar)           && (return           Planar(body1, body2, axis; p1 = p1, p2 = p2, qoffset = qoffset, spring = spring, damper = damper))
    (jointtype == :FixedOrientation) && (return FixedOrientation(body1, body2;                         qoffset = qoffset, spring = spring, damper = damper))
    (jointtype == :Revolute)         && (return         Revolute(body1, body2, axis; p1 = p1, p2 = p2, qoffset = qoffset, spring = spring, damper = damper))
    (jointtype == :Cylindrical)      && (return      Cylindrical(body1, body2, axis; p1 = p1, p2 = p2, qoffset = qoffset, spring = spring, damper = damper))
    (jointtype == :PlanarAxis)       && (return       PlanarAxis(body1, body2, axis; p1 = p1, p2 = p2, qoffset = qoffset, spring = spring, damper = damper))
    (jointtype == :FreeRevolute)     && (return     FreeRevolute(body1, body2, axis; p1 = p1, p2 = p2, qoffset = qoffset, spring = spring, damper = damper))
    (jointtype == :Orbital)          && (return          Orbital(body1, body2, axis; p1 = p1, p2 = p2, qoffset = qoffset, spring = spring, damper = damper))
    (jointtype == :PrismaticOrbital) && (return PrismaticOrbital(body1, body2, axis; p1 = p1, p2 = p2, qoffset = qoffset, spring = spring, damper = damper))
    (jointtype == :PlanarOrbital)    && (return    PlanarOrbital(body1, body2, axis; p1 = p1, p2 = p2, qoffset = qoffset, spring = spring, damper = damper))
    (jointtype == :FreeOrbital)      && (return      FreeOrbital(body1, body2, axis; p1 = p1, p2 = p2, qoffset = qoffset, spring = spring, damper = damper))
    (jointtype == :Spherical)        && (return        Spherical(body1, body2;       p1 = p1, p2 = p2, qoffset = qoffset, spring = spring, damper = damper))
    (jointtype == :CylindricalFree)  && (return  CylindricalFree(body1, body2, axis; p1 = p1, p2 = p2,                    spring = spring, damper = damper))
    (jointtype == :PlanarFree)       && (return       PlanarFree(body1, body2, axis; p1 = p1, p2 = p2,                    spring = spring, damper = damper))
    (jointtype == :Floating)         && (return         Floating(body1, body2;                                            spring = spring, damper = damper))
end
