using Dojo

vis = Visualizer()
open(vis)

function get_tugbot(;
    timestep=0.01,
    gravity=[0.0; 0.0; -9.81],
    friction_coefficient=0.8,
    radius=0.1,
    object_dimension=[0.2, 0.2, 0.05],
    contact=true,
    contact_type=:nonlinear,
    T=Float64)

    origin = Origin{T}(name=:origin)
    drone_mass = 0.3
    object_mass = 1.0
    bodies = [
        Sphere(radius, drone_mass, name=:drone),
        Box(object_dimension..., object_mass, name=:object)]
    joints = [
        JointConstraint(Floating(origin, bodies[1]), name=:drone_joint)
        JointConstraint(Floating(origin, bodies[2]), name=:object_joint)]
    mechanism = Mechanism(origin, bodies, joints,
        timestep=timestep,
        gravity=gravity)

    # if contact
    #     contact = [0.0, 0.0, 0.0]
    #     normal = [0.0, 0.0, 1.0]
    #     contacts = [contact_constraint(get_body(mechanism, :tugbot), normal,
    #         friction_coefficient=friction_coefficient,
    #         contact_origin=contact,
    #         contact_radius=radius,
    #         contact_type=contact_type)]
    #     set_minimal_coordinates!(mechanism, get_joint(mechanism, :floating_joint), [0.0; 0.0; radius; zeros(3)])
    #     mechanism = Mechanism(origin, bodies, joints, contacts,
    #         gravity=gravity,
    #         timestep=timestep)
    # end
    return mechanism
end

function initialize_tugbot!(mechanism::Mechanism{T};
    drone_position=zeros(3),
    drone_orientation=one(Quaternion),
    drone_velocity=zeros(3),
    drone_angular_velocity=zeros(3),
    object_position=zeros(3),
    object_orientation=one(Quaternion),
    object_velocity=zeros(3),
    object_angular_velocity=zeros(3),
    ) where T

    r = get_body(mechanism, :drone).shape.r
    joint = get_joint(mechanism, :drone_joint)
    zero_velocity!(mechanism)
    z_drone = [drone_position + [0,0,r]; vector(drone_orientation); drone_velocity; drone_angular_velocity]
    z_object = [object_position + [0,0,r]; vector(object_orientation); object_velocity; object_angular_velocity]
    set_maximal_state!(mechanism, [z_drone; z_object])
end


mech = get_tugbot(gravity=0.0, timestep=0.01)
initialize!(mech, :tugbot)
storage = simulate!(mech, 0.01, record=true, verbose=true)
visualize(mech, storage, vis=vis)
