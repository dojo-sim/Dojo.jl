using Dojo

vis = Visualizer()
open(vis)

function get_tugbot(;
    timestep=0.01,
    gravity=[0.0; 0.0; -9.81],
    friction_coefficient=0.8,
    radius=0.2,
    object_dimension=[0.4, 0.4, 0.1],
    contact=true,
    contact_type=:nonlinear,
    T=Float64)

    origin = Origin{T}(name=:origin)
    drone_mass = 1.0
    object_mass = 10.0
    bodies = [
        Sphere(radius, drone_mass, name=:drone),
        Box(object_dimension..., object_mass, name=:object)]
    joints = [
        JointConstraint(Floating(origin, bodies[1]), name=:drone_joint)
        JointConstraint(Floating(origin, bodies[2]), name=:object_joint)]
    mechanism = Mechanism(origin, bodies, joints,
        timestep=timestep,
        gravity=gravity)

    if contact
        contacts = [[0.2, -0.2, -0.05], [-0.2, 0.2, -0.05], [0.2, 0.2, -0.05], [-0.2, -0.2, -0.05]]
        normal = [0.0, 0.0, 1.0]
        contacts = [contact_constraint(get_body(mechanism, :object), normal,
            friction_coefficient=friction_coefficient,
            contact_origin=contacts[i],
            contact_radius=0.0,
            contact_type=contact_type,
            name=Symbol(:contact,i)) for i=1:length(contacts)]
        mechanism = Mechanism(origin, bodies, joints, contacts,
            gravity=gravity,
            timestep=timestep)
    end
    return mechanism
end

function initialize_tugbot!(mechanism::Mechanism{T};
    drone_position=[0,1,1.0],
    drone_orientation=one(Quaternion),
    drone_velocity=zeros(3),
    drone_angular_velocity=zeros(3),
    object_position=zeros(3),
    object_orientation=one(Quaternion),
    object_velocity=zeros(3),
    object_angular_velocity=zeros(3),
    ) where T

    r_drone = get_body(mechanism, :drone).shape.r
    r_object = get_body(mechanism, :object).shape.xyz[3]/2
    joint = get_joint(mechanism, :drone_joint)
    zero_velocity!(mechanism)
    z_drone = [drone_position + [0,0,r_drone]; drone_velocity; vector(drone_orientation); drone_angular_velocity]
    z_object = [object_position + [0,0,r_object]; object_velocity; vector(object_orientation); object_angular_velocity]
    set_maximal_state!(mechanism, [z_drone; z_object])
end


mech = get_tugbot(gravity=-9.81, timestep=0.1)

function ctrl!(mechanism::Mechanism, k::Int; x_target=SVector{3}(1,1,1.0), kp=2.0, kd=2.0)
    dt = mechanism.timestep
    drone_body = get_body(mechanism, :drone)
    drone_joint = get_joint(mechanism, :drone_joint)
    x = current_position(drone_body.state)
    v = current_velocity(drone_body.state)[1]
    u_gravity = -drone_body.mass * mechanism.gravity * dt
    u_tra = u_gravity + kp*(x_target - x)*dt - kd * dt .* v
    u_rot = szeros(3)
    set_input!(drone_joint, [u_tra; u_rot])
    return nothing
end

initialize!(mech, :tugbot)
storage = simulate!(mech, 5.0, ctrl!, record=true, verbose=true)
visualize(mech, storage, vis=vis, show_contact=true)
