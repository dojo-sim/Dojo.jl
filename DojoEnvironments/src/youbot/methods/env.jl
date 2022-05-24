"""
    Youbot <: Environment

    Youbot robot designed by Kuka (URDF and mesh files from https://github.com/mas-group/youbot_description)
"""
struct Youbot end

using Dojo


path = "DojoEnvironments/src/youbot/deps/youbot_arm_only.urdf"
# path = "DojoEnvironments/src/youbot/deps/youbot_base_only.urdf"
# path = "DojoEnvironments/src/youbot/deps/youbot.urdf"
# path = "DojoEnvironments/src/youbot/deps/youbot_dual_arm.urdf"
# path = "DojoEnvironments/src/youbot/deps/youbot_with_cam3d.urdf"
mech = Mechanism(path, gravity = -9.81, timestep = 0.001)

# joint1 = get_joint(mech,:link1_joint)
# joint2 = get_joint(mech,:link2_joint)
# joint3 = get_joint(mech,:link3_joint)
# joint4 = get_joint(mech,:link4_joint)

# set_minimal_coordinates!(mech,joint1,[pi/2])
# set_minimal_coordinates!(mech,joint2,[0])
# set_minimal_coordinates!(mech,joint3,[-pi/2])
# set_minimal_coordinates!(mech,joint4,[pi/2])

storage = simulate!(mech, 5.0, record = true)
visualize(mech, storage)
 