function get_slider(;
    timestep=0.01,
    gravity=-9.81,
    spring=0.0,
    damper=0.0,
    T=Float64)

    # Parameters
    joint_axis = [0.0; 0.0; 1.0]
    len = 1.0
    width, depth = 0.1, 0.1
    child_vertex = [0.0; 0.0; len / 2.0] # joint connection point

    # Links
    origin = Origin{T}()
    pbody = Box(width, depth, len, len, name=:slider)

    # Constraints
    joint_between_origin_and_pbody = JointConstraint(Prismatic(origin, pbody, joint_axis;
        child_vertex=child_vertex,
        spring=spring,
        damper=damper))

    bodies = [pbody]
    joints = [joint_between_origin_and_pbody]

    mech = Mechanism(origin, bodies, joints,
        gravity=gravity,
        timestep=timestep,
        spring=spring,
        damper=damper)

    return mech
end

function initialize_slider!(mechanism::Mechanism{T};
        position=0.0) where T

    zero_velocity!(mechanism)
    body = mechanism.bodies[1]
    set_maximal_configurations!(body, x=[0,0,position], q=Quaternion(1,0,0,0.0))
end

function get_nslider(;
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
    vert11 = [0.0; r; 0.0]
    vert12 = -vert11

    # Links
    origin = Origin{T}()
    bodies = [Cylinder(r, h, h, color=RGBA(1.0, 0.0, 0.0)) for i = 1:num_bodies]

    # Constraints
    jointb1 = JointConstraint(Prismatic(origin, bodies[1], ex;
        child_vertex = 0 * vert11))

    if num_bodies > 1
        joints = [
            jointb1;
            [JointConstraint(Prismatic(bodies[i - 1], bodies[i], ex;
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
        timestep=timestep)

    return mech
end

function initialize_nslider!(mechanism::Mechanism{T};
    position=0.2,
    relative_position=0.0) where T

    pbody = mechanism.bodies[1]

    # set position and velocities
    set_maximal_configurations!(mechanism.origin, pbody,
        parent_vertex=[0.0, 0.0, position])

    # set relative positions
    previd = pbody.id
    for body in mechanism.bodies[2:end]
        set_maximal_configurations!(get_body(mechanism, previd), body,
            parent_vertex=[0.0, -0.1, relative_position])
        previd = body.id
    end

    zero_velocity!(mechanism)
end
