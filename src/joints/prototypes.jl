# fixed connection between two bodies.
Fixed(body1::Component{T}, body2::Component{T}; p1=szeros(T, 3), p2=szeros(T, 3),
    qoffset=one(UnitQuaternion{T})) where T =
    Translational3{T}(body1, body2; p1, p2), Rotational3{T}(body1, body2; qoffset)

# prismatic joint between two bodies.
Prismatic(body1::Component{T}, body2::Component{T}, axis; p1=szeros(T, 3), p2=szeros(T, 3),
    qoffset=one(UnitQuaternion{T}), spring=zero(T), damper=zero(T),
    tra_spring_offset=szeros(T,1)) where T =
    Translational2{T}(body1, body2; p1, p2, axis, spring, damper, spring_offset=tra_spring_offset),
    Rotational3{T}(body1, body2; qoffset, spring, damper)

# planar joint between two bodies.
Planar(body1::Component{T}, body2::Component{T}, axis; p1=szeros(T, 3), p2=szeros(T, 3),
    qoffset=one(UnitQuaternion{T}), spring=zero(T), damper=zero(T),
    tra_spring_offset=szeros(T,2)) where T =
    Translational1{T}(body1, body2; p1, p2, axis, spring, damper, spring_offset=tra_spring_offset),
    Rotational3{T}(body1, body2; qoffset, spring, damper)

# fixed orientation between two bodies (chicken's head).
FixedOrientation(body1::Component{T}, body2::Component{T}; qoffset=one(UnitQuaternion{T}),
    spring=zero(T), damper=zero(T),
    tra_spring_offset=szeros(T,3)) where T =
    Translational0{T}(body1, body2; spring, damper, spring_offset=tra_spring_offset),
    Rotational3{T}(body1, body2; qoffset, spring, damper)

# revolute joint between two bodies (pin, continuous, hinge joint).
Revolute(body1::Component{T}, body2::Component{T}, axis; p1=szeros(T, 3), p2=szeros(T, 3),
    qoffset=one(UnitQuaternion{T}), spring=zero(T), damper=zero(T),
    rot_spring_offset=szeros(T,1)) where T =
    Translational3{T}(body1, body2; p1, p2, spring, damper),
    Rotational2{T}(body1, body2; axis, qoffset, spring, damper, spring_offset=rot_spring_offset)

# cylindrical joint between two bodies.
Cylindrical(body1::Component{T}, body2::Component{T}, axis; p1=szeros(T, 3), p2=szeros(T, 3),
    qoffset=one(UnitQuaternion{T}), spring=zero(T), damper=zero(T),
    tra_spring_offset=szeros(T,1), rot_spring_offset=szeros(T,1)) where T =
    Translational2{T}(body1, body2; p1, p2, axis, spring, damper, spring_offset=tra_spring_offset),
    Rotational2{T}(body1, body2; axis, qoffset, spring, damper, spring_offset=rot_spring_offset)

# planar joint between two bodies with a rotation axis perpendicular to the plane (turtle bot).
PlanarAxis(body1::Component{T}, body2::Component{T}, axis; p1=szeros(T, 3), p2=szeros(T, 3),
    qoffset = one(UnitQuaternion{T}), spring=zero(T), damper=zero(T),
    tra_spring_offset=szeros(T,2), rot_spring_offset=szeros(T,1)) where T =
    Translational1{T}(body1, body2; p1, p2, axis, spring, damper, spring_offset=tra_spring_offset),
    Rotational2{T}(body1, body2; axis, qoffset, spring, damper, spring_offset=rot_spring_offset)

# joint between two bodies with free translation and rotation along one axis.
FreeRevolute(body1::Component{T}, body2::Component{T}, axis; p1=szeros(T, 3), p2=szeros(T, 3),
    qoffset=one(UnitQuaternion{T}), spring=zero(T), damper=zero(T),
    tra_spring_offset=szeros(T,3), rot_spring_offset=szeros(T,1)) where T =
    Translational0{T}(body1, body2; spring, damper, spring_offset=tra_spring_offset),
    Rotational2{T}(body1, body2; axis, qoffset, spring, damper, spring_offset=rot_spring_offset)

# rotational between two bodies with a 2 rotational degrees of freedom (skull-eye joint).
Orbital(body1::Component{T}, body2::Component{T}, axis; p1=szeros(T, 3), p2=szeros(T, 3),
    qoffset=one(UnitQuaternion{T}), spring=zero(T), damper=zero(T),
    rot_spring_offset=szeros(T,2)) where T =
    Translational3{T}(body1, body2; p1, p2, axis, spring, damper),
    Rotational1{T}(body1, body2; axis, qoffset, spring, damper, spring_offset=rot_spring_offset)

# prismatic joint between two bodies with a 2 rotational degrees of freedom (skull-eye joint).
PrismaticOrbital(body1::Component{T}, body2::Component{T}, axis; p1=szeros(T, 3), p2=szeros(T, 3),
    qoffset = one(UnitQuaternion{T}), spring=zero(T), damper=zero(T),
    tra_spring_offset=szeros(T,1), rot_spring_offset=szeros(T,2)) where T =
    Translational2{T}(body1, body2; p1, p2, axis, spring, damper, spring_offset=tra_spring_offset),
    Rotational1{T}(body1, body2; axis, qoffset, spring, damper, spring_offset=rot_spring_offset)

# planar joint between two bodies with a 2 rotational degrees of freedom (skull-eye joint).
PlanarOrbital(body1::Component{T}, body2::Component{T}, axis; p1=szeros(T, 3), p2=szeros(T, 3),
    qoffset=one(UnitQuaternion{T}), spring=zero(T), damper=zero(T),
    tra_spring_offset=szeros(T,2), rot_spring_offset=szeros(T,2)) where T =
    Translational1{T}(body1, body2; p1, p2, axis, spring, damper, spring_offset=tra_spring_offset),
    Rotational1{T}(body1, body2; axis, qoffset, spring, damper, spring_offset=rot_spring_offset)

# free joint between two bodies with a 2 rotational degrees of freedom (skull-eye joint).
FreeOrbital(body1::Component{T}, body2::Component{T}, axis; p1=szeros(T, 3), p2=szeros(T, 3),
    qoffset=one(UnitQuaternion{T}), spring=zero(T), damper=zero(T),
    tra_spring_offset=szeros(T,3), rot_spring_offset=szeros(T,2)) where T =
    Translational0{T}(body1, body2; spring, damper, spring_offset=tra_spring_offset),
    Rotational1{T}(body1, body2; axis, qoffset, spring, damper, spring_offset=rot_spring_offset)

# spherical joint between two bodies (ball-and-socket joint).
Spherical(body1::Component{T}, body2::Component{T}; p1=szeros(T, 3), p2=szeros(T, 3),
    qoffset = one(UnitQuaternion{T}), spring=zero(T), damper=zero(T),
    rot_spring_offset=szeros(T,3)) where T =
    Translational3{T}(body1, body2; p1, p2, spring, damper),
    Rotational0{T}(body1, body2; qoffset, spring, damper, spring_offset=rot_spring_offset)

# cylindrical joint between two bodies with unconstrained orientation (point-on-line).
CylindricalFree(body1::Component{T}, body2::Component{T}, axis; p1=szeros(T, 3), p2=szeros(T, 3),
    spring=zero(T), damper=zero(T),
    tra_spring_offset=szeros(T,1), rot_spring_offset=szeros(T,3)) where T =
    Translational2{T}(body1, body2; p1, p2, axis, spring, damper, spring_offset=tra_spring_offset),
    Rotational0{T}(body1, body2; spring, damper, spring_offset=rot_spring_offset)

# planar joint between two bodies with unconstrained orientation.
PlanarFree(body1::Component{T}, body2::Component{T}, axis; p1=szeros(T, 3), p2=szeros(T, 3),
    spring = zero(T), damper= zero(T),
    tra_spring_offset=szeros(T,2), rot_spring_offset=szeros(T,3)) where T =
    Translational1{T}(body1, body2; p1, p2, axis, spring, damper, spring_offset=tra_spring_offset),
    Rotational0{T}(body1, body2; spring, damper, spring_offset=rot_spring_offset)

# unconstrained connection between two bodies (connection between floating base and origin).
Floating(body1::Component{T}, body2::Component{T}; spring=zero(T), damper=zero(T),
    tra_spring_offset=szeros(T,3), rot_spring_offset=szeros(T,3)) where T =
    Translational0{T}(body1, body2; spring, damper, spring_offset=tra_spring_offset),
    Rotational0{T}(body1, body2; spring, damper, spring_offset=rot_spring_offset)

function Prototype(jointtype::Symbol, body1::Component{T}, body2::Component{T}, axis; p1=szeros(T, 3), p2=szeros(T, 3),
        qoffset=one(UnitQuaternion{T}), spring=zero(T), damper=zero(T),
        tra_spring_offset=nothing, rot_spring_offset=nothing) where T
    N̄tra, N̄rot = nullspacedims(jointtype)
    (tra_spring_offset == nothing) && (tra_spring_offset = szeros(T,N̄tra))
    (rot_spring_offset == nothing) && (rot_spring_offset = szeros(T,N̄rot))
    (jointtype == :Fixed)            && (return            Fixed(body1, body2;       p1=p1, p2=p2, qoffset=qoffset))
    (jointtype == :Prismatic)        && (return        Prismatic(body1, body2, axis; p1=p1, p2=p2, qoffset=qoffset, spring=spring, damper=damper, tra_spring_offset=tra_spring_offset))
    (jointtype == :Planar)           && (return           Planar(body1, body2, axis; p1=p1, p2=p2, qoffset=qoffset, spring=spring, damper=damper, tra_spring_offset=tra_spring_offset))
    (jointtype == :FixedOrientation) && (return FixedOrientation(body1, body2;                     qoffset=qoffset, spring=spring, damper=damper, tra_spring_offset=tra_spring_offset))
    (jointtype == :Revolute)         && (return         Revolute(body1, body2, axis; p1=p1, p2=p2, qoffset=qoffset, spring=spring, damper=damper,                                      rot_spring_offset=rot_spring_offset))
    (jointtype == :Cylindrical)      && (return      Cylindrical(body1, body2, axis; p1=p1, p2=p2, qoffset=qoffset, spring=spring, damper=damper, tra_spring_offset=tra_spring_offset, rot_spring_offset=rot_spring_offset))
    (jointtype == :PlanarAxis)       && (return       PlanarAxis(body1, body2, axis; p1=p1, p2=p2, qoffset=qoffset, spring=spring, damper=damper, tra_spring_offset=tra_spring_offset, rot_spring_offset=rot_spring_offset))
    (jointtype == :FreeRevolute)     && (return     FreeRevolute(body1, body2, axis; p1=p1, p2=p2, qoffset=qoffset, spring=spring, damper=damper, tra_spring_offset=tra_spring_offset, rot_spring_offset=rot_spring_offset))
    (jointtype == :Orbital)          && (return          Orbital(body1, body2, axis; p1=p1, p2=p2, qoffset=qoffset, spring=spring, damper=damper,                                      rot_spring_offset=rot_spring_offset))
    (jointtype == :PrismaticOrbital) && (return PrismaticOrbital(body1, body2, axis; p1=p1, p2=p2, qoffset=qoffset, spring=spring, damper=damper, tra_spring_offset=tra_spring_offset, rot_spring_offset=rot_spring_offset))
    (jointtype == :PlanarOrbital)    && (return    PlanarOrbital(body1, body2, axis; p1=p1, p2=p2, qoffset=qoffset, spring=spring, damper=damper, tra_spring_offset=tra_spring_offset, rot_spring_offset=rot_spring_offset))
    (jointtype == :FreeOrbital)      && (return      FreeOrbital(body1, body2, axis; p1=p1, p2=p2, qoffset=qoffset, spring=spring, damper=damper, tra_spring_offset=tra_spring_offset, rot_spring_offset=rot_spring_offset))
    (jointtype == :Spherical)        && (return        Spherical(body1, body2;       p1=p1, p2=p2, qoffset=qoffset, spring=spring, damper=damper,                                      rot_spring_offset=rot_spring_offset))
    (jointtype == :CylindricalFree)  && (return  CylindricalFree(body1, body2, axis; p1=p1, p2=p2,                  spring=spring, damper=damper, tra_spring_offset=tra_spring_offset, rot_spring_offset=rot_spring_offset))
    (jointtype == :PlanarFree)       && (return       PlanarFree(body1, body2, axis; p1=p1, p2=p2,                  spring=spring, damper=damper, tra_spring_offset=tra_spring_offset, rot_spring_offset=rot_spring_offset))
    (jointtype == :Floating)         && (return         Floating(body1, body2;                                      spring=spring, damper=damper, tra_spring_offset=tra_spring_offset, rot_spring_offset=rot_spring_offset))
end

function nullspacedims(jointtype::Symbol)
    (jointtype == :Fixed)            && (return 0, 0)
    (jointtype == :Prismatic)        && (return 1, 0)
    (jointtype == :Planar)           && (return 2, 0)
    (jointtype == :FixedOrientation) && (return 3, 0)
    (jointtype == :Revolute)         && (return 0, 1)
    (jointtype == :Cylindrical)      && (return 1, 1)
    (jointtype == :PlanarAxis)       && (return 2, 1)
    (jointtype == :FreeRevolute)     && (return 3, 1)
    (jointtype == :Orbital)          && (return 0, 2)
    (jointtype == :PrismaticOrbital) && (return 1, 2)
    (jointtype == :PlanarOrbital)    && (return 2, 2)
    (jointtype == :FreeOrbital)      && (return 3, 2)
    (jointtype == :Spherical)        && (return 0, 3)
    (jointtype == :CylindricalFree)  && (return 1, 3)
    (jointtype == :PlanarFree)       && (return 2, 3)
    (jointtype == :Floating)         && (return 3, 3)
end
