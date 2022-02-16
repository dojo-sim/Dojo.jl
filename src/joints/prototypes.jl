# fixed connection between two bodies.
Fixed(body1::Node{T}, body2::Node{T}; p1=szeros(T, 3), p2=szeros(T, 3),
    qoffset=one(UnitQuaternion{T})) where T =
    Translational{T,3}(body1, body2; p1, p2),
    Rotational{T,3}(body1, body2; qoffset)

# prismatic joint between two bodies.
Prismatic(body1::Node{T}, body2::Node{T}, axis; p1=szeros(T, 3), p2=szeros(T, 3),
    qoffset=one(UnitQuaternion{T}), spring=zero(T), damper=zero(T),
    tra_spring_offset=szeros(T,1),
    tra_joint_limits=[szeros(T,0), szeros(T,0)]) where T =
    Translational{T,2}(body1, body2; p1, p2, axis, spring, damper,
        spring_offset=tra_spring_offset, joint_limits=tra_joint_limits),
    Rotational{T,3}(body1, body2; qoffset, spring, damper)

# planar joint between two bodies.
Planar(body1::Node{T}, body2::Node{T}, axis; p1=szeros(T, 3), p2=szeros(T, 3),
    qoffset=one(UnitQuaternion{T}), spring=zero(T), damper=zero(T),
    tra_spring_offset=szeros(T,2),
    tra_joint_limits=[szeros(T,0), szeros(T,0)]) where T =
    Translational{T,1}(body1, body2; p1, p2, axis, spring, damper,
        spring_offset=tra_spring_offset, joint_limits=tra_joint_limits),
    Rotational{T,3}(body1, body2; qoffset, spring, damper)

# fixed orientation between two bodies (chicken's head).
FixedOrientation(body1::Node{T}, body2::Node{T}; qoffset=one(UnitQuaternion{T}),
    spring=zero(T), damper=zero(T),
    tra_spring_offset=szeros(T,3),
    tra_joint_limits=[szeros(T,0), szeros(T,0)]) where T =
    Translational{T,0}(body1, body2; spring, damper,
        spring_offset=tra_spring_offset, joint_limits=tra_joint_limits),
    Rotational{T,3}(body1, body2; qoffset, spring, damper)

# revolute joint between two bodies (pin, continuous, hinge joint).
Revolute(body1::Node{T}, body2::Node{T}, axis; p1=szeros(T, 3), p2=szeros(T, 3),
    qoffset=one(UnitQuaternion{T}), spring=zero(T), damper=zero(T),
    rot_spring_offset=szeros(T,1),
    rot_joint_limits=[szeros(T,0), szeros(T,0)],
    spring_type=:sinusoidal) where T =
    Translational{T,3}(body1, body2; p1, p2, spring, damper),
    Rotational{T,2}(body1, body2; axis, qoffset, spring, damper,
        spring_offset=rot_spring_offset, joint_limits=rot_joint_limits, spring_type=spring_type)

# cylindrical joint between two bodies.
Cylindrical(body1::Node{T}, body2::Node{T}, axis; p1=szeros(T, 3), p2=szeros(T, 3),
    qoffset=one(UnitQuaternion{T}), spring=zero(T), damper=zero(T),
    tra_spring_offset=szeros(T,1), rot_spring_offset=szeros(T,1),
    rot_joint_limits=[szeros(T,0), szeros(T,0)], tra_joint_limits=[szeros(T,0), szeros(T,0)],
    spring_type=:sinusoidal) where T =
    Translational{T,2}(body1, body2; p1, p2, axis, spring, damper,
        spring_offset=tra_spring_offset, joint_limits=tra_joint_limits),
    Rotational{T,2}(body1, body2; axis, qoffset, spring, damper,
        spring_offset=rot_spring_offset, joint_limits=rot_joint_limits, spring_type=spring_type)

# planar joint between two bodies with a rotation axis perpendicular to the plane (turtle bot).
PlanarAxis(body1::Node{T}, body2::Node{T}, axis; p1=szeros(T, 3), p2=szeros(T, 3),
    qoffset=one(UnitQuaternion{T}), spring=zero(T), damper=zero(T),
    tra_spring_offset=szeros(T,2), rot_spring_offset=szeros(T,1),
    rot_joint_limits=[szeros(T,0), szeros(T,0)], tra_joint_limits=[szeros(T,0), szeros(T,0)],
    spring_type=:sinusoidal) where T =
    Translational{T,1}(body1, body2; p1, p2, axis, spring, damper,
        spring_offset=tra_spring_offset, joint_limits=tra_joint_limits),
    Rotational{T,2}(body1, body2; axis, qoffset, spring, damper,
        spring_offset=rot_spring_offset, joint_limits=rot_joint_limits, spring_type=spring_type)

# joint between two bodies with free translation and rotation along one axis.
FreeRevolute(body1::Node{T}, body2::Node{T}, axis; p1=szeros(T, 3), p2=szeros(T, 3),
    qoffset=one(UnitQuaternion{T}), spring=zero(T), damper=zero(T),
    tra_spring_offset=szeros(T,3), rot_spring_offset=szeros(T,1),
    rot_joint_limits=[szeros(T,0), szeros(T,0)], tra_joint_limits=[szeros(T,0), szeros(T,0)],
    spring_type=:sinusoidal) where T =
    Translational{T,0}(body1, body2; spring, damper,
        spring_offset=tra_spring_offset, joint_limits=tra_joint_limits),
    Rotational{T,2}(body1, body2; axis, qoffset, spring, damper,
        spring_offset=rot_spring_offset, joint_limits=rot_joint_limits, spring_type=spring_type)

# rotational between two bodies with a 2 rotational degrees of freedom (skull-eye joint).
Orbital(body1::Node{T}, body2::Node{T}, axis; p1=szeros(T, 3), p2=szeros(T, 3),
    qoffset=one(UnitQuaternion{T}), spring=zero(T), damper=zero(T),
    rot_spring_offset=szeros(T,2),
    rot_joint_limits=[szeros(T,0), szeros(T,0)],
    spring_type=:sinusoidal) where T =
    Translational{T,3}(body1, body2; p1, p2, axis, spring, damper),
    Rotational{T,1}(body1, body2; axis, qoffset, spring, damper,
        spring_offset=rot_spring_offset, joint_limits=rot_joint_limits, spring_type=spring_type)

# prismatic joint between two bodies with a 2 rotational degrees of freedom (skull-eye joint).
PrismaticOrbital(body1::Node{T}, body2::Node{T}, axis; p1=szeros(T, 3), p2=szeros(T, 3),
    qoffset=one(UnitQuaternion{T}), spring=zero(T), damper=zero(T),
    tra_spring_offset=szeros(T,1), rot_spring_offset=szeros(T,2),
    rot_joint_limits=[szeros(T,0), szeros(T,0)], tra_joint_limits=[szeros(T,0), szeros(T,0)],
    spring_type=:sinusoidal) where T =
    Translational{T,2}(body1, body2; p1, p2, axis, spring, damper,
        spring_offset=tra_spring_offset, joint_limits=tra_joint_limits),
    Rotational{T,1}(body1, body2; axis, qoffset, spring, damper,
        spring_offset=rot_spring_offset, joint_limits=rot_joint_limits, spring_type=spring_type)

# planar joint between two bodies with a 2 rotational degrees of freedom (skull-eye joint).
PlanarOrbital(body1::Node{T}, body2::Node{T}, axis; p1=szeros(T, 3), p2=szeros(T, 3),
    qoffset=one(UnitQuaternion{T}), spring=zero(T), damper=zero(T),
    tra_spring_offset=szeros(T,2), rot_spring_offset=szeros(T,2),
    rot_joint_limits=[szeros(T,0), szeros(T,0)], tra_joint_limits=[szeros(T,0), szeros(T,0)],
    spring_type=:sinusoidal) where T =
    Translational{T,1}(body1, body2; p1, p2, axis, spring, damper,
        spring_offset=tra_spring_offset, joint_limits=tra_joint_limits),
    Rotational{T,1}(body1, body2; axis, qoffset, spring, damper,
        spring_offset=rot_spring_offset, joint_limits=rot_joint_limits, spring_type=spring_type)

# free joint between two bodies with a 2 rotational degrees of freedom (skull-eye joint).
FreeOrbital(body1::Node{T}, body2::Node{T}, axis; p1=szeros(T, 3), p2=szeros(T, 3),
    qoffset=one(UnitQuaternion{T}), spring=zero(T), damper=zero(T),
    tra_spring_offset=szeros(T,3), rot_spring_offset=szeros(T,2),
    rot_joint_limits=[szeros(T,0), szeros(T,0)], tra_joint_limits=[szeros(T,0), szeros(T,0)],
    spring_type=:sinusoidal) where T =
    Translational{T,0}(body1, body2; spring, damper,
        spring_offset=tra_spring_offset, joint_limits=tra_joint_limits),
    Rotational{T,1}(body1, body2; axis, qoffset, spring, damper,
        spring_offset=rot_spring_offset, joint_limits=rot_joint_limits, spring_type=spring_type)

# spherical joint between two bodies (ball-and-socket joint).
Spherical(body1::Node{T}, body2::Node{T}; p1=szeros(T, 3), p2=szeros(T, 3),
    qoffset=one(UnitQuaternion{T}), spring=zero(T), damper=zero(T),
    rot_spring_offset=szeros(T,3), rot_joint_limits=[szeros(T,0), szeros(T,0)],
    spring_type=:sinusoidal) where T =
    Translational{T,3}(body1, body2; p1, p2, spring, damper),
    Rotational{T,0}(body1, body2; qoffset, spring, damper,
        spring_offset=rot_spring_offset, joint_limits=rot_joint_limits, spring_type=spring_type)

# cylindrical joint between two bodies with unconstrained orientation (point-on-line).
CylindricalFree(body1::Node{T}, body2::Node{T}, axis; p1=szeros(T, 3), p2=szeros(T, 3),
    spring=zero(T), damper=zero(T),
    tra_spring_offset=szeros(T,1), rot_spring_offset=szeros(T,3),
    rot_joint_limits=[szeros(T,0), szeros(T,0)], tra_joint_limits=[szeros(T,0), szeros(T,0)],
    spring_type=:sinusoidal) where T =
    Translational{T,2}(body1, body2; p1, p2, axis, spring, damper,
        spring_offset=tra_spring_offset, joint_limits=tra_joint_limits),
    Rotational{T,0}(body1, body2; spring, damper,
        spring_offset=rot_spring_offset, joint_limits=rot_joint_limits, spring_type=spring_type)

# planar joint between two bodies with unconstrained orientation.
PlanarFree(body1::Node{T}, body2::Node{T}, axis; p1=szeros(T, 3), p2=szeros(T, 3),
    spring=zero(T), damper= zero(T),
    tra_spring_offset=szeros(T,2), rot_spring_offset=szeros(T,3),
    rot_joint_limits=[szeros(T,0), szeros(T,0)], tra_joint_limits=[szeros(T,0), szeros(T,0)],
    spring_type=:sinusoidal) where T =
    Translational{T,1}(body1, body2; p1, p2, axis,
        spring, damper, spring_offset=tra_spring_offset, joint_limits=tra_joint_limits),
    Rotational{T,0}(body1, body2; spring, damper,
        spring_offset=rot_spring_offset, joint_limits=rot_joint_limits, spring_type=spring_type)

# unconstrained connection between two bodies (connection between floating base and origin).
Floating(body1::Node{T}, body2::Node{T}; spring=zero(T), damper=zero(T),
    tra_spring_offset=szeros(T,3), rot_spring_offset=szeros(T,3),
    rot_joint_limits=[szeros(T,0), szeros(T,0)], tra_joint_limits=[szeros(T,0), szeros(T,0)],
    spring_type=:sinusoidal) where T =
    Translational{T,0}(body1, body2; spring, damper,
        spring_offset=tra_spring_offset, joint_limits=tra_joint_limits),
    Rotational{T,0}(body1, body2; spring, damper,
        spring_offset=rot_spring_offset, joint_limits=rot_joint_limits, spring_type=spring_type)

function Prototype(jointtype::Symbol, body1::Node{T}, body2::Node{T}, axis; p1=szeros(T, 3), p2=szeros(T, 3),
        qoffset=one(UnitQuaternion{T}), spring=zero(T), damper=zero(T),
        tra_spring_offset=nothing, rot_spring_offset=nothing,
        tra_joint_limits=[szeros(T,0), szeros(T,0)], rot_joint_limits=[szeros(T,0), szeros(T,0)],
        spring_type=:sinusoidal) where T

    N̄tra, N̄rot = nullspace_dimension(jointtype)
    (tra_spring_offset == nothing) && (tra_spring_offset = szeros(T,N̄tra))
    (rot_spring_offset == nothing) && (rot_spring_offset = szeros(T,N̄rot))
    (jointtype == :Fixed)            && (return            Fixed(body1, body2;       p1=p1, p2=p2, qoffset=qoffset))
    (jointtype == :Prismatic)        && (return        Prismatic(body1, body2, axis; p1=p1, p2=p2, qoffset=qoffset, spring=spring, damper=damper, tra_spring_offset=tra_spring_offset,                                      tra_joint_limits=tra_joint_limits))
    (jointtype == :Planar)           && (return           Planar(body1, body2, axis; p1=p1, p2=p2, qoffset=qoffset, spring=spring, damper=damper, tra_spring_offset=tra_spring_offset,                                      tra_joint_limits=tra_joint_limits))
    (jointtype == :FixedOrientation) && (return FixedOrientation(body1, body2;                     qoffset=qoffset, spring=spring, damper=damper, tra_spring_offset=tra_spring_offset,                                      tra_joint_limits=tra_joint_limits))
    (jointtype == :Revolute)         && (return         Revolute(body1, body2, axis; p1=p1, p2=p2, qoffset=qoffset, spring=spring, damper=damper,                                      rot_spring_offset=rot_spring_offset,                                    rot_joint_limits=rot_joint_limits, spring_type=spring_type))
    (jointtype == :Cylindrical)      && (return      Cylindrical(body1, body2, axis; p1=p1, p2=p2, qoffset=qoffset, spring=spring, damper=damper, tra_spring_offset=tra_spring_offset, rot_spring_offset=rot_spring_offset, tra_joint_limits=tra_joint_limits, rot_joint_limits=rot_joint_limits, spring_type=spring_type))
    (jointtype == :PlanarAxis)       && (return       PlanarAxis(body1, body2, axis; p1=p1, p2=p2, qoffset=qoffset, spring=spring, damper=damper, tra_spring_offset=tra_spring_offset, rot_spring_offset=rot_spring_offset, tra_joint_limits=tra_joint_limits, rot_joint_limits=rot_joint_limits, spring_type=spring_type))
    (jointtype == :FreeRevolute)     && (return     FreeRevolute(body1, body2, axis; p1=p1, p2=p2, qoffset=qoffset, spring=spring, damper=damper, tra_spring_offset=tra_spring_offset, rot_spring_offset=rot_spring_offset, tra_joint_limits=tra_joint_limits, rot_joint_limits=rot_joint_limits, spring_type=spring_type))
    (jointtype == :Orbital)          && (return          Orbital(body1, body2, axis; p1=p1, p2=p2, qoffset=qoffset, spring=spring, damper=damper,                                      rot_spring_offset=rot_spring_offset,                                    rot_joint_limits=rot_joint_limits, spring_type=spring_type))
    (jointtype == :PrismaticOrbital) && (return PrismaticOrbital(body1, body2, axis; p1=p1, p2=p2, qoffset=qoffset, spring=spring, damper=damper, tra_spring_offset=tra_spring_offset, rot_spring_offset=rot_spring_offset, tra_joint_limits=tra_joint_limits, rot_joint_limits=rot_joint_limits, spring_type=spring_type))
    (jointtype == :PlanarOrbital)    && (return    PlanarOrbital(body1, body2, axis; p1=p1, p2=p2, qoffset=qoffset, spring=spring, damper=damper, tra_spring_offset=tra_spring_offset, rot_spring_offset=rot_spring_offset, tra_joint_limits=tra_joint_limits, rot_joint_limits=rot_joint_limits, spring_type=spring_type))
    (jointtype == :FreeOrbital)      && (return      FreeOrbital(body1, body2, axis; p1=p1, p2=p2, qoffset=qoffset, spring=spring, damper=damper, tra_spring_offset=tra_spring_offset, rot_spring_offset=rot_spring_offset, tra_joint_limits=tra_joint_limits, rot_joint_limits=rot_joint_limits, spring_type=spring_type))
    (jointtype == :Spherical)        && (return        Spherical(body1, body2;       p1=p1, p2=p2, qoffset=qoffset, spring=spring, damper=damper,                                      rot_spring_offset=rot_spring_offset,                                    rot_joint_limits=rot_joint_limits, spring_type=spring_type))
    (jointtype == :CylindricalFree)  && (return  CylindricalFree(body1, body2, axis; p1=p1, p2=p2,                  spring=spring, damper=damper, tra_spring_offset=tra_spring_offset, rot_spring_offset=rot_spring_offset, tra_joint_limits=tra_joint_limits, rot_joint_limits=rot_joint_limits, spring_type=spring_type))
    (jointtype == :PlanarFree)       && (return       PlanarFree(body1, body2, axis; p1=p1, p2=p2,                  spring=spring, damper=damper, tra_spring_offset=tra_spring_offset, rot_spring_offset=rot_spring_offset, tra_joint_limits=tra_joint_limits, rot_joint_limits=rot_joint_limits, spring_type=spring_type))
    (jointtype == :Floating)         && (return         Floating(body1, body2;                                      spring=spring, damper=damper, tra_spring_offset=tra_spring_offset, rot_spring_offset=rot_spring_offset, tra_joint_limits=tra_joint_limits, rot_joint_limits=rot_joint_limits, spring_type=spring_type))
end

function nullspace_dimension(jointtype::Symbol)
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
