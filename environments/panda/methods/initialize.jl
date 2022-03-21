
function get_panda(;
        timestep=0.01,
        gravity=-9.81,
        friction_coefficient=0.8,
        spring=0.0,
        damper=0.01, # this value comes from the official URDF https://github.com/frankaemika/franka_ros/blob/develop/franka_gazebo/test/launch/panda-gazebo.urdf
        contact=false,
        T=Float64)

    path = joinpath(@__DIR__, "../deps/panda_simple.urdf")
    mech = Mechanism(path, false, T,
        gravity=gravity,
        timestep=timestep,
        spring=spring,
        damper=damper)

    # Adding springs and dampers
    for joint in mech.joints
        joint.damper = true
        joint.spring = true
        joint.translational.spring=spring
        joint.translational.damper=damper
        joint.rotational.spring=spring
        joint.rotational.damper=damper
    end

    origin = Origin{T}()
    bodies = mech.bodies
    eqs = mech.joints

    contacts = ContactConstraint{T}[]

    # if contact
    #     # Foot contact
    #     locations = [
    #         [-0.08; -0.04; 0.015],
    #         [+0.12; -0.02; 0.015],
    #         [-0.08; +0.04; 0.015],
    #         [+0.12; +0.02; 0.015],
    #         ]
    #     n = length(locations)
    #     normal = [[0.0; 0.0; 1.0] for i = 1:n]
    #     offset = [[0.0; 0.0; 0.025] for i = 1:n]
    #     friction_coefficients = friction_coefficient * ones(T, n)
    #
    #     names = ["RR", "FR", "RL", "RR"]
    #
    #     foot_contacts1 = contact_constraint(get_body(mech, :l_foot), normal,
    #         friction_coefficient=friction_coefficients,
    #         contact_points=locations,
    #         offset=offset,
    #         names=[Symbol("l_" .* name) for name in names])
    #     foot_contacts2 = contact_constraint(get_body(mech, :r_foot), normal,
    #         friction_coefficient=friction_coefficients,
    #         contact_points=locations,
    #         offset=offset,
    #         names=[Symbol("r_" .* name) for name in names])
    #     contacts = [contacts..., foot_contacts1..., foot_contacts2...]
    # end

    # set_minimal_coordinates!(mech, get_joint(mech, :panda_joint0), [0.0; 0.0])

    # mech = Mechanism(origin, bodies, eqs, contacts,
    #     gravity=gravity,
    #     timestep=timestep,
    #     spring=spring,
    #     damper=damper)

    return mech
end

function initialize_panda!(mechanism::Mechanism;
    model_type=:simple,
    body_position=[0.0, 0.0, 0.2],
    body_orientation=[0.0, 0.0, 0.0],
    link_linear_velocity=[zeros(3)  for i=1:length(mechanism.bodies)],
    link_angular_velocity=[zeros(3) for i=1:length(mechanism.bodies)],
    hip_orientation=0.0,
    knee_orientation=0.0) where T

    body_position += (model_type == :armless) ? [0.0, 0.0, 0.9385 + 0.14853] : [0.0, 0.0, 0.9385]

    # positions
    try
        set_minimal_coordinates!(mechanism,
                get_joint(mechanism, :floating_base),
                [body_position; body_orientation])
        set_minimal_coordinates!(mechanism, get_joint(mechanism, :back_bkxyz), [0.0, 0.0, 0.0])
        set_minimal_coordinates!(mechanism, get_joint(mechanism, :l_leg_hpxyz), [0.0, -hip_orientation, 0.0])
        set_minimal_coordinates!(mechanism, get_joint(mechanism, :r_leg_hpxyz), [0.0, -hip_orientation, 0.0])
        set_minimal_coordinates!(mechanism, get_joint(mechanism, :l_leg_kny), [knee_orienation])
        set_minimal_coordinates!(mechanism, get_joint(mechanism, :r_leg_kny), [knee_orienation])
        set_minimal_coordinates!(mechanism, get_joint(mechanism, :l_leg_akxy), [hip_orientation-knee_orienation, 0.0])
        set_minimal_coordinates!(mechanism, get_joint(mechanism, :r_leg_akxy), [hip_orientation-knee_orienation, 0.0])
    catch
        nothing
    end

    zero_velocity!(mechanism)

    return nothing
end
