"""
    Fixed{T} <: JointConstraint{T}

    fixed connection between two bodies
"""
Fixed(pbody::Node{T}, cbody::Node{T}; 
    parent_vertex=szeros(T, 3), 
    child_vertex=szeros(T, 3),
    orientation_offset=one(Quaternion{T})) where T =
    Translational{T,3}(pbody, cbody; 
        parent_vertex, 
        child_vertex),
    Rotational{T,3}(pbody, cbody; 
        orientation_offset)

"""
    Prismatic{T} <: JointConstraint{T}

    one translational degree of freedom between two bodies
"""
Prismatic(pbody::Node{T}, cbody::Node{T}, axis; 
    parent_vertex=szeros(T, 3), 
    child_vertex=szeros(T, 3),
    orientation_offset=one(Quaternion{T}), 
    spring=zero(T), 
    damper=zero(T),
    tra_spring_offset=szeros(T,1),
    tra_joint_limits=[szeros(T,0), szeros(T,0)]) where T =
    Translational{T,2}(pbody, cbody; 
        parent_vertex, 
        child_vertex, 
        axis, 
        spring, 
        damper,
        spring_offset=tra_spring_offset, 
        joint_limits=tra_joint_limits),
    Rotational{T,3}(pbody, cbody; 
        orientation_offset, 
        spring, 
        damper)

"""
    Planar{T} <: JointConstraint{T} 

    two translational degree of freedom between two bodies
"""
Planar(pbody::Node{T}, cbody::Node{T}, axis; 
    parent_vertex=szeros(T, 3), 
    child_vertex=szeros(T, 3),
    orientation_offset=one(Quaternion{T}), 
    spring=zero(T), 
    damper=zero(T),
    tra_spring_offset=szeros(T,2),
    tra_joint_limits=[szeros(T,0), szeros(T,0)]) where T =
    Translational{T,1}(pbody, cbody; 
        parent_vertex, 
        child_vertex, 
        axis, 
        spring, 
        damper,
        spring_offset=tra_spring_offset, 
        joint_limits=tra_joint_limits),
    Rotational{T,3}(pbody, cbody; 
        orientation_offset, 
        spring, 
        damper)

"""
    FixedOrientation{T} <: JointConstraint{T} 

    three translational degree of freedom between two bodies
"""
FixedOrientation(pbody::Node{T}, cbody::Node{T}; 
    orientation_offset=one(Quaternion{T}),
    spring=zero(T), 
    damper=zero(T),
    tra_spring_offset=szeros(T,3),
    tra_joint_limits=[szeros(T,0), szeros(T,0)]) where T =
    Translational{T,0}(pbody, cbody; 
        spring, 
        damper,
        spring_offset=tra_spring_offset, 
        joint_limits=tra_joint_limits),
    Rotational{T,3}(pbody, cbody; 
        orientation_offset, 
        spring, 
        damper)

"""
    Revolute{T} <: JointConstraint{T} 

    one rotational degree of freedom between two bodies
"""
Revolute(pbody::Node{T}, cbody::Node{T}, axis; 
    parent_vertex=szeros(T, 3), 
    child_vertex=szeros(T, 3),
    orientation_offset=one(Quaternion{T}), 
    spring=zero(T), 
    damper=zero(T),
    rot_spring_offset=szeros(T,1),
    rot_joint_limits=[szeros(T,0), szeros(T,0)],
    spring_type=:linear) where T =
    Translational{T,3}(pbody, cbody; 
        parent_vertex, 
        child_vertex, 
        spring, 
        damper),
    Rotational{T,2}(pbody, cbody; 
        axis, 
        orientation_offset, 
        spring, 
        damper,
        spring_offset=rot_spring_offset, 
        joint_limits=rot_joint_limits, 
        spring_type=spring_type)

"""
    Cylindrical{T} <: JointConstraint{T} 

    one translational and one rotational degree of freedom between two bodies
"""
Cylindrical(pbody::Node{T}, cbody::Node{T}, axis; 
    parent_vertex=szeros(T, 3), 
    child_vertex=szeros(T, 3),
    orientation_offset=one(Quaternion{T}), 
    spring=zero(T), 
    damper=zero(T),
    tra_spring_offset=szeros(T,1), 
    rot_spring_offset=szeros(T,1),
    rot_joint_limits=[szeros(T,0), szeros(T,0)], 
    tra_joint_limits=[szeros(T,0), szeros(T,0)],
    spring_type=:linear) where T =
    Translational{T,2}(pbody, cbody; 
        parent_vertex, 
        child_vertex, 
        axis, 
        spring, 
        damper,
        spring_offset=tra_spring_offset, 
        joint_limits=tra_joint_limits),
    Rotational{T,2}(pbody, cbody; 
        axis, 
        orientation_offset, 
        spring, 
        damper,
        spring_offset=rot_spring_offset, 
        joint_limits=rot_joint_limits, 
        spring_type=spring_type)

"""
    PlanarAxis{T} <: JointConstraint{T} 

    two translational and one rotational degree of freedom between two bodies
"""
PlanarAxis(pbody::Node{T}, cbody::Node{T}, axis; 
    parent_vertex=szeros(T, 3), 
    child_vertex=szeros(T, 3),
    orientation_offset=one(Quaternion{T}), 
    spring=zero(T), 
    damper=zero(T),
    tra_spring_offset=szeros(T,2), 
    rot_spring_offset=szeros(T,1),
    rot_joint_limits=[szeros(T,0), szeros(T,0)], 
    tra_joint_limits=[szeros(T,0), szeros(T,0)],
    spring_type=:linear) where T =
    Translational{T,1}(pbody, cbody; 
        parent_vertex, 
        child_vertex, 
        axis, 
        spring, 
        damper,
        spring_offset=tra_spring_offset, 
        joint_limits=tra_joint_limits),
    Rotational{T,2}(pbody, cbody; 
        axis, 
        orientation_offset, 
        spring, 
        damper,
        spring_offset=rot_spring_offset, 
        joint_limits=rot_joint_limits, 
        spring_type=spring_type)

"""
    FreeRevolute{T} <: JointConstraint{T} 

    free translation with rotation along one axis
"""
FreeRevolute(pbody::Node{T}, cbody::Node{T}, axis; 
    parent_vertex=szeros(T, 3), 
    child_vertex=szeros(T, 3),
    orientation_offset=one(Quaternion{T}), 
    spring=zero(T), 
    damper=zero(T),
    tra_spring_offset=szeros(T,3), 
    rot_spring_offset=szeros(T,1),
    rot_joint_limits=[szeros(T,0), szeros(T,0)], 
    tra_joint_limits=[szeros(T,0), szeros(T,0)],
    spring_type=:linear) where T =
    Translational{T,0}(pbody, cbody; 
        spring, 
        damper,
        spring_offset=tra_spring_offset, 
        joint_limits=tra_joint_limits),
    Rotational{T,2}(pbody, cbody; 
        axis, 
        orientation_offset, 
        spring, 
        damper,
        spring_offset=rot_spring_offset, 
        joint_limits=rot_joint_limits, 
        spring_type=spring_type)

"""
    Orbital{T} <: JointConstraint{T} 

    two rotational degrees of freedom between two bodies
"""
Orbital(pbody::Node{T}, cbody::Node{T}, axis; 
    parent_vertex=szeros(T, 3), 
    child_vertex=szeros(T, 3),
    orientation_offset=one(Quaternion{T}), 
    spring=zero(T), 
    damper=zero(T),
    rot_spring_offset=szeros(T,2),
    rot_joint_limits=[szeros(T,0), szeros(T,0)],
    spring_type=:linear) where T =
    Translational{T,3}(pbody, cbody; 
        parent_vertex, 
        child_vertex, 
        axis, 
        spring, 
        damper),
    Rotational{T,1}(pbody, cbody; 
        axis, 
        orientation_offset, 
        spring, 
        damper,
        spring_offset=rot_spring_offset, 
        joint_limits=rot_joint_limits, 
        spring_type=spring_type)

"""
    PrismaticOrbital{T} <: JointConstraint{T} 

    one translational and two rotational degrees of freedom between two bodies
"""
PrismaticOrbital(pbody::Node{T}, cbody::Node{T}, axis; 
    parent_vertex=szeros(T, 3), 
    child_vertex=szeros(T, 3),
    orientation_offset=one(Quaternion{T}), 
    spring=zero(T), 
    damper=zero(T),
    tra_spring_offset=szeros(T,1), 
    rot_spring_offset=szeros(T,2),
    rot_joint_limits=[szeros(T,0), szeros(T,0)], 
    tra_joint_limits=[szeros(T,0), szeros(T,0)],
    spring_type=:linear) where T =
    Translational{T,2}(pbody, cbody; 
        parent_vertex, 
        child_vertex, 
        axis, 
        spring, 
        damper,
        spring_offset=tra_spring_offset, 
        joint_limits=tra_joint_limits),
    Rotational{T,1}(pbody, cbody; 
        axis, 
        orientation_offset, 
        spring, 
        damper,
        spring_offset=rot_spring_offset, 
        joint_limits=rot_joint_limits, 
        spring_type=spring_type)

"""
    PlanarOrbital{T} <: JointConstraint{T} 

    two translational and two rotational degrees of freedom between two bodies
"""
PlanarOrbital(pbody::Node{T}, cbody::Node{T}, axis; 
    parent_vertex=szeros(T, 3), 
    child_vertex=szeros(T, 3),
    orientation_offset=one(Quaternion{T}), 
    spring=zero(T), 
    damper=zero(T),
    tra_spring_offset=szeros(T,2), 
    rot_spring_offset=szeros(T,2),
    rot_joint_limits=[szeros(T,0), szeros(T,0)], 
    tra_joint_limits=[szeros(T,0), szeros(T,0)],
    spring_type=:linear) where T =
    Translational{T,1}(pbody, cbody; 
        parent_vertex, 
        child_vertex, 
        axis, 
        spring, 
        damper,
        spring_offset=tra_spring_offset, 
        joint_limits=tra_joint_limits),
    Rotational{T,1}(pbody, cbody; 
        axis, 
        orientation_offset, 
        spring, 
        damper,
        spring_offset=rot_spring_offset, 
        joint_limits=rot_joint_limits, 
        spring_type=spring_type)

"""
    FreeOrbital{T} <: JointConstraint{T} 

    three translational and two rotational degrees of freedom between two bodies
"""
FreeOrbital(pbody::Node{T}, cbody::Node{T}, axis; 
    parent_vertex=szeros(T, 3), 
    child_vertex=szeros(T, 3),
    orientation_offset=one(Quaternion{T}), 
    spring=zero(T), 
    damper=zero(T),
    tra_spring_offset=szeros(T,3), 
    rot_spring_offset=szeros(T,2),
    rot_joint_limits=[szeros(T,0), szeros(T,0)], 
    tra_joint_limits=[szeros(T,0), szeros(T,0)],
    spring_type=:linear) where T =
    Translational{T,0}(pbody, cbody; 
        spring, 
        damper,
        spring_offset=tra_spring_offset, 
        joint_limits=tra_joint_limits),
    Rotational{T,1}(pbody, cbody; 
        axis, 
        orientation_offset, 
        spring, 
        damper,
        spring_offset=rot_spring_offset, 
        joint_limits=rot_joint_limits, 
        spring_type=spring_type)

"""
    Spherical{T} <: JointConstraint{T} 

    three rotational degrees of freedom between two bodies
"""
Spherical(pbody::Node{T}, cbody::Node{T}; 
    parent_vertex=szeros(T, 3), 
    child_vertex=szeros(T, 3),
    orientation_offset=one(Quaternion{T}), 
    spring=zero(T), 
    damper=zero(T),
    rot_spring_offset=szeros(T,3), 
    rot_joint_limits=[szeros(T,0), szeros(T,0)],
    spring_type=:linear) where T =
    Translational{T,3}(pbody, cbody; 
        parent_vertex, 
        child_vertex, 
        spring, 
        damper),
    Rotational{T,0}(pbody, cbody; 
        orientation_offset, 
        spring, 
        damper,
        spring_offset=rot_spring_offset, 
        joint_limits=rot_joint_limits, 
        spring_type=spring_type)

"""
    CylindricalFree{T} <: JointConstraint{T} 

    one translational and three rotational degrees of freedom between two bodies
"""
CylindricalFree(pbody::Node{T}, cbody::Node{T}, axis; 
    parent_vertex=szeros(T, 3), 
    child_vertex=szeros(T, 3),
    spring=zero(T), damper=zero(T),
    tra_spring_offset=szeros(T,1), 
    rot_spring_offset=szeros(T,3),
    rot_joint_limits=[szeros(T,0), szeros(T,0)], 
    tra_joint_limits=[szeros(T,0), szeros(T,0)],
    spring_type=:linear) where T =
    Translational{T,2}(pbody, cbody; 
        parent_vertex, 
        child_vertex, 
        axis, 
        spring, 
        damper,
        spring_offset=tra_spring_offset, 
        joint_limits=tra_joint_limits),
    Rotational{T,0}(pbody, cbody; 
        spring, 
        damper,
        spring_offset=rot_spring_offset, 
        joint_limits=rot_joint_limits, 
        spring_type=spring_type)

"""
    PlanarFree{T} <: JointConstraint{T} 

    two translational and three rotational degrees of freedom between two bodies
"""
PlanarFree(pbody::Node{T}, cbody::Node{T}, axis; 
    parent_vertex=szeros(T, 3), 
    child_vertex=szeros(T, 3),
    spring=zero(T), 
    damper= zero(T),
    tra_spring_offset=szeros(T,2), 
    rot_spring_offset=szeros(T,3),
    rot_joint_limits=[szeros(T,0), szeros(T,0)], 
    tra_joint_limits=[szeros(T,0), szeros(T,0)],
    spring_type=:linear) where T =
    Translational{T,1}(pbody, cbody; 
        parent_vertex, 
        child_vertex, 
        axis,
        spring, 
        damper, 
        spring_offset=tra_spring_offset, 
        joint_limits=tra_joint_limits),
    Rotational{T,0}(pbody, cbody; 
        spring, 
        damper,
        spring_offset=rot_spring_offset, 
        joint_limits=rot_joint_limits, 
        spring_type=spring_type)

"""
    Floating{T} <: JointConstraint{T} 

    no restricted degrees of freedom between two bodies
"""
Floating(pbody::Node{T}, cbody::Node{T}; 
    spring=zero(T), 
    damper=zero(T),
    tra_spring_offset=szeros(T,3), 
    rot_spring_offset=szeros(T,3),
    rot_joint_limits=[szeros(T,0), szeros(T,0)], 
    tra_joint_limits=[szeros(T,0), szeros(T,0)],
    spring_type=:linear) where T =
    Translational{T,0}(pbody, cbody; 
        spring, 
        damper,
        spring_offset=tra_spring_offset, 
        joint_limits=tra_joint_limits),
    Rotational{T,0}(pbody, cbody; 
        spring, 
        damper,
        spring_offset=rot_spring_offset, 
        joint_limits=rot_joint_limits, 
        spring_type=spring_type)

function Prototype(joint_type::Symbol, pbody::Node{T}, cbody::Node{T}, axis; 
        parent_vertex=szeros(T, 3), 
        child_vertex=szeros(T, 3),
        orientation_offset=one(Quaternion{T}), 
        spring=zero(T), 
        damper=zero(T),
        tra_spring_offset=nothing, 
        rot_spring_offset=nothing,
        tra_joint_limits=[szeros(T,0), szeros(T,0)], 
        rot_joint_limits=[szeros(T,0), szeros(T,0)],
        spring_type=:linear) where T

    N̄tra, N̄rot = nullspace_dimension(joint_type)
    (tra_spring_offset == nothing) && (tra_spring_offset = szeros(T,N̄tra))
    (rot_spring_offset == nothing) && (rot_spring_offset = szeros(T,N̄rot))
    (joint_type == :Fixed)            && (return            Fixed(pbody, cbody;       parent_vertex=parent_vertex, child_vertex=child_vertex, orientation_offset=orientation_offset))
    (joint_type == :Prismatic)        && (return        Prismatic(pbody, cbody, axis; parent_vertex=parent_vertex, child_vertex=child_vertex, orientation_offset=orientation_offset, spring=spring, damper=damper, tra_spring_offset=tra_spring_offset,                                      tra_joint_limits=tra_joint_limits))
    (joint_type == :Planar)           && (return           Planar(pbody, cbody, axis; parent_vertex=parent_vertex, child_vertex=child_vertex, orientation_offset=orientation_offset, spring=spring, damper=damper, tra_spring_offset=tra_spring_offset,                                      tra_joint_limits=tra_joint_limits))
    (joint_type == :FixedOrientation) && (return FixedOrientation(pbody, cbody;                     orientation_offset=orientation_offset, spring=spring, damper=damper, tra_spring_offset=tra_spring_offset,                                      tra_joint_limits=tra_joint_limits))
    (joint_type == :Revolute)         && (return         Revolute(pbody, cbody, axis; parent_vertex=parent_vertex, child_vertex=child_vertex, orientation_offset=orientation_offset, spring=spring, damper=damper,                                      rot_spring_offset=rot_spring_offset,                                    rot_joint_limits=rot_joint_limits, spring_type=spring_type))
    (joint_type == :Cylindrical)      && (return      Cylindrical(pbody, cbody, axis; parent_vertex=parent_vertex, child_vertex=child_vertex, orientation_offset=orientation_offset, spring=spring, damper=damper, tra_spring_offset=tra_spring_offset, rot_spring_offset=rot_spring_offset, tra_joint_limits=tra_joint_limits, rot_joint_limits=rot_joint_limits, spring_type=spring_type))
    (joint_type == :PlanarAxis)       && (return       PlanarAxis(pbody, cbody, axis; parent_vertex=parent_vertex, child_vertex=child_vertex, orientation_offset=orientation_offset, spring=spring, damper=damper, tra_spring_offset=tra_spring_offset, rot_spring_offset=rot_spring_offset, tra_joint_limits=tra_joint_limits, rot_joint_limits=rot_joint_limits, spring_type=spring_type))
    (joint_type == :FreeRevolute)     && (return     FreeRevolute(pbody, cbody, axis; parent_vertex=parent_vertex, child_vertex=child_vertex, orientation_offset=orientation_offset, spring=spring, damper=damper, tra_spring_offset=tra_spring_offset, rot_spring_offset=rot_spring_offset, tra_joint_limits=tra_joint_limits, rot_joint_limits=rot_joint_limits, spring_type=spring_type))
    (joint_type == :Orbital)          && (return          Orbital(pbody, cbody, axis; parent_vertex=parent_vertex, child_vertex=child_vertex, orientation_offset=orientation_offset, spring=spring, damper=damper,                                      rot_spring_offset=rot_spring_offset,                                    rot_joint_limits=rot_joint_limits, spring_type=spring_type))
    (joint_type == :PrismaticOrbital) && (return PrismaticOrbital(pbody, cbody, axis; parent_vertex=parent_vertex, child_vertex=child_vertex, orientation_offset=orientation_offset, spring=spring, damper=damper, tra_spring_offset=tra_spring_offset, rot_spring_offset=rot_spring_offset, tra_joint_limits=tra_joint_limits, rot_joint_limits=rot_joint_limits, spring_type=spring_type))
    (joint_type == :PlanarOrbital)    && (return    PlanarOrbital(pbody, cbody, axis; parent_vertex=parent_vertex, child_vertex=child_vertex, orientation_offset=orientation_offset, spring=spring, damper=damper, tra_spring_offset=tra_spring_offset, rot_spring_offset=rot_spring_offset, tra_joint_limits=tra_joint_limits, rot_joint_limits=rot_joint_limits, spring_type=spring_type))
    (joint_type == :FreeOrbital)      && (return      FreeOrbital(pbody, cbody, axis; parent_vertex=parent_vertex, child_vertex=child_vertex, orientation_offset=orientation_offset, spring=spring, damper=damper, tra_spring_offset=tra_spring_offset, rot_spring_offset=rot_spring_offset, tra_joint_limits=tra_joint_limits, rot_joint_limits=rot_joint_limits, spring_type=spring_type))
    (joint_type == :Spherical)        && (return        Spherical(pbody, cbody;       parent_vertex=parent_vertex, child_vertex=child_vertex, orientation_offset=orientation_offset, spring=spring, damper=damper,                                      rot_spring_offset=rot_spring_offset,                                    rot_joint_limits=rot_joint_limits, spring_type=spring_type))
    (joint_type == :CylindricalFree)  && (return  CylindricalFree(pbody, cbody, axis; parent_vertex=parent_vertex, child_vertex=child_vertex,                  spring=spring, damper=damper, tra_spring_offset=tra_spring_offset, rot_spring_offset=rot_spring_offset, tra_joint_limits=tra_joint_limits, rot_joint_limits=rot_joint_limits, spring_type=spring_type))
    (joint_type == :PlanarFree)       && (return       PlanarFree(pbody, cbody, axis; parent_vertex=parent_vertex, child_vertex=child_vertex,                  spring=spring, damper=damper, tra_spring_offset=tra_spring_offset, rot_spring_offset=rot_spring_offset, tra_joint_limits=tra_joint_limits, rot_joint_limits=rot_joint_limits, spring_type=spring_type))
    (joint_type == :Floating)         && (return         Floating(pbody, cbody;                                      spring=spring, damper=damper, tra_spring_offset=tra_spring_offset, rot_spring_offset=rot_spring_offset, tra_joint_limits=tra_joint_limits, rot_joint_limits=rot_joint_limits, spring_type=spring_type))
end

function nullspace_dimension(joint_type::Symbol)
    (joint_type == :Fixed)            && (return 0, 0)
    (joint_type == :Prismatic)        && (return 1, 0)
    (joint_type == :Planar)           && (return 2, 0)
    (joint_type == :FixedOrientation) && (return 3, 0)
    (joint_type == :Revolute)         && (return 0, 1)
    (joint_type == :Cylindrical)      && (return 1, 1)
    (joint_type == :PlanarAxis)       && (return 2, 1)
    (joint_type == :FreeRevolute)     && (return 3, 1)
    (joint_type == :Orbital)          && (return 0, 2)
    (joint_type == :PrismaticOrbital) && (return 1, 2)
    (joint_type == :PlanarOrbital)    && (return 2, 2)
    (joint_type == :FreeOrbital)      && (return 3, 2)
    (joint_type == :Spherical)        && (return 0, 3)
    (joint_type == :CylindricalFree)  && (return 1, 3)
    (joint_type == :PlanarFree)       && (return 2, 3)
    (joint_type == :Floating)         && (return 3, 3)
end
