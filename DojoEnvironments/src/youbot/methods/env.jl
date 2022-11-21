"""
    Youbot <: Environment

    Youbot robot designed by Kuka (URDF and mesh files from https://github.com/mas-group/youbot_description)
"""
struct Youbot end

using Dojo


# path = "DojoEnvironments/src/youbot/deps/youbot_arm_only.urdf"
# path = "DojoEnvironments/src/youbot/deps/youbot_arm_only_fixed_gripper.urdf"
# path = "DojoEnvironments/src/youbot/deps/youbot_base_only.urdf"
path = "DojoEnvironments/src/youbot/deps/youbot.urdf"
# path = "DojoEnvironments/src/youbot/deps/youbot_dual_arm.urdf"
# path = "DojoEnvironments/src/youbot/deps/youbot_with_cam3d.urdf"
mech = Mechanism(path; gravity = -9.81, timestep = 0.01)

joint1 = get_joint(mech,:base_footprint_joint)

set_minimal_coordinates!(mech,joint1,[0;0;1.5])

function controller!(mechanism,k)
    # set_input!(joint1, [0.1;0.1;0.4])
end

storage = simulate!(mech, 5.0, controller!, record = true)
visualize(mech, storage)[1]


 