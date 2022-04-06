function get_orbital(; 
    timestep=0.01, 
    gravity=-9.81, 
    spring=0.0, 
    damper=0.0, 
    num_bodies=5,
    T=Float64)

    # Parameters
    ex = [0.0; 0.0; 1.0]
    h = 1.0
    r = 0.05
    vert11 = [0.0; 0.0; h / 2.0]
    vert12 = -vert11

    # Links
    origin = Origin{T}()

    bodies = [Box(r, r, h, h, color=RGBA(1.0, 0.0, 0.0)) for i = 1:num_bodies]

    # Constraints
    jointb1 = JointConstraint(Fixed(origin, bodies[1]; 
        child_vertex=vert11))
    if num_bodies > 1
        joints = [
            jointb1;
            [
                JointConstraint(Orbital(bodies[i - 1], bodies[i], ex; 
                parent_vertex=vert12, 
                child_vertex=vert11, 
                spring=spring, 
                damper=damper)) for i = 2:num_bodies]
            ]
    else
        joints = [jointb1]
    end

    mech = Mechanism(origin, bodies, joints, 
        gravity=gravity, 
        timestep=timestep, 
        spring=spring, 
        damper=damper)

    return mech
end

function initialize_orbital!(mechanism::Mechanism{T}; 
    orientation=[pi / 4.0, pi / 8.0]) where T

    pbody = mechanism.bodies[1]
    joint = mechanism.joints[1]
    vert11 = joint.translational.vertices[2]
    vert12 = -vert11

    # set position and velocities
    set_maximal_configurations!(mechanism.origin, pbody, 
        child_vertex=vert11, 
        Δq=RotX(0.0))

    set_minimal_coordinates!(mechanism, mechanism.joints[2], orientation)

    zero_velocity!(mechanism)
end
