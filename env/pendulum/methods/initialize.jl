function get_pendulum(; timestep::T = 0.01, gravity = -9.81, m::T = 1.0, l::T = 1.0,
        spring = 0.0, damper = 0.0, spring_offset = szeros(1), qoffset=one(UnitQuaternion)) where T
    # Parameters
    joint_axis = [1.0; 0; 0]
    width, depth = 0.1, 0.1
    p2 = [0; 0; l/2] # joint connection point

    # Links
    origin = Origin{T}()
    body1 = Box(width, depth, l, m)

    # Constraints
    joint_between_origin_and_body1 = JointConstraint(Revolute(origin, body1,
        joint_axis; 
        p2=p2, 
        spring=spring, 
        damper=damper, 
        rot_spring_offset=spring_offset,
        qoffset=qoffset))
    bodies = [body1]
    joints = [joint_between_origin_and_body1]

    mech = Mechanism(origin, bodies, joints, gravity=gravity, timestep=timestep, spring=spring, damper=damper)
    return mech
end

function initialize_pendulum!(mechanism::Mechanism; ϕ1::T = 0.7, ω1::T = 0.0) where T
    body = collect(mechanism.bodies)[1]
    joint = collect(mechanism.joints)[1]
    set_minimal_coordinates_velocities!(mechanism, joint; xmin=[ϕ1, ω1])
end

function get_npendulum(; timestep::T = 0.01, gravity = -9.81, m::T = 1.0, l::T = 1.0,
        spring = 0.0, damper::T = 0.0, Nb::Int = 5,
        basetype::Symbol = :Revolute, jointtype::Symbol = :Revolute) where T
    # Parameters
    ex = [1.; 0; 0]
    r = 0.05
    vert11 = [0; 0; l/2]
    vert12 = -vert11

    # Links
    origin = Origin{T}()
    bodies = [Box(r, r, l, m, color = RGBA(1., 0., 0.)) for i = 1:Nb]

    # Constraints
    jointb1 = JointConstraint(Prototype(basetype, origin, bodies[1], ex; p2 = vert11, spring=spring, damper=damper))
    if Nb > 1
        joints = [JointConstraint(Prototype(jointtype, bodies[i - 1], bodies[i], ex; p1 = vert12, p2 = vert11, spring=spring, damper=damper)) for i = 2:Nb]
        joints = [jointb1; joints]
    else
        joints = [jointb1]
    end

    mech = Mechanism(origin, bodies, joints, gravity=gravity, timestep=timestep)
    return mech
end

function initialize_npendulum!(mechanism::Mechanism; ϕ1::T = pi/4, ω = [0.0, 0.0, 0.0],
    Δv::AbstractVector{T} = [0, 0, 0.], Δω::AbstractVector{T} = [0, 0, 0.]) where T

    body1 = mechanism.bodies[1]
    joint = mechanism.joints[1]
    vert11 = joint.translational.vertices[2]
    vert12 = - vert11

    # set position and velocities
    set_position!(mechanism.origin, body1, p2 = vert11, Δq = UnitQuaternion(RotX(ϕ1)))
    set_velocity!(body1, ω = ω)

    previd = body1.id
    for (i,body) in enumerate(Iterators.drop(mechanism.bodies, 1))
        set_position!(get_body(mechanism, previd), body, p1 = vert12, p2 = vert11)
        set_velocity!(get_body(mechanism, previd), body, p1 = vert12, p2 = vert11,
                Δv = Δv, Δω = 1/i*Δω)
        previd = body.id
    end
end
