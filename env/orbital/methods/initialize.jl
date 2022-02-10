function get_orbital(; timestep::T = 0.01, gravity = -9.81, spring = 0.0, damper = 0.0, Nb::Int = 5) where T
    # Parameters
    ex = [0; 0; 1]
    h = 1.
    r = 0.05
    vert11 = [0; 0; h/2]
    vert12 = -vert11

    # Links
    origin = Origin{T}()

    bodies = [Box(r, r, h, h, color = RGBA(1., 0., 0.)) for i = 1:Nb]

    # Constraints
    jointb1 = JointConstraint(Fixed(origin, bodies[1]; p2 = vert11))
    if Nb > 1
        joints = [
            jointb1;
            [JointConstraint(Orbital(bodies[i - 1], bodies[i], ex; p1=vert12, p2=vert11, spring=spring, damper=damper)) for i = 2:Nb]
            ]
    else
        joints = [jointb1]
    end
    mech = Mechanism(origin, bodies, joints, gravity=gravity, timestep=timestep, spring=spring, damper=damper)
    return mech
end

function initialize_orbital!(mechanism::Mechanism; ϕx::T = pi/4, ϕy::T = pi/8) where T
    body1 = collect(mechanism.bodies)[1]
    joint = collect(mechanism.joints)[1]
    vert11 = joint.translational.vertices[2]
    vert12 = - vert11

    # set position and velocities
    set_position!(mechanism.origin, body1, p2 = vert11, Δq = UnitQuaternion(RotX(0.0)))

    previd = body1.id
    set_position!(mechanism, collect(mechanism.joints)[2], [ϕx, ϕy])

    zero_velocity!(mechanism)
end
