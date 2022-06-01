"""
    Youbot <: Environment

    Youbot robot designed by Kuka (URDF and mesh files from https://github.com/mas-group/youbot_description)
"""
struct Youbot end

using Dojo


# path = "DojoEnvironments/src/youbot/deps/youbot_arm_only.urdf"
path = "DojoEnvironments/src/youbot/deps/youbot_arm_only_fixed_gripper.urdf"
# path = "DojoEnvironments/src/youbot/deps/youbot_base_only.urdf"
# path = "DojoEnvironments/src/youbot/deps/youbot.urdf"
# path = "DojoEnvironments/src/youbot/deps/youbot_dual_arm.urdf"
# path = "DojoEnvironments/src/youbot/deps/youbot_with_cam3d.urdf"
mech = Mechanism(path; gravity = -9.81, timestep = 0.01)

joint1 = get_joint(mech,:arm_joint_1)
joint2 = get_joint(mech,:arm_joint_2)
joint3 = get_joint(mech,:arm_joint_3)
joint4 = get_joint(mech,:arm_joint_4)
joint5 = get_joint(mech,:arm_joint_5)

set_minimal_coordinates!(mech,joint1,[0])
set_minimal_coordinates!(mech,joint2,[0])
set_minimal_coordinates!(mech,joint3,[0])
set_minimal_coordinates!(mech,joint4,[0])
set_minimal_coordinates!(mech,joint5,[0])

storage = simulate!(mech, 5.0, record = true)
visualize(mech, storage)
 