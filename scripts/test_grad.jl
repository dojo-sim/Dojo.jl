using Dojo, StaticArrays

m = 10.0

pad_r = 0.5
pad_m = m * 4/3 * pi * pad_r ^ 3

box_x = 20.0
box_y = 1.0
box_z = 1.0
box_m = m * box_x * box_y * box_z
μ = 0.50

spring = 2000.0
damper = 0.0

origin = Origin()
pad = Dojo.Sphere(pad_r, pad_m)
box = Dojo.Box(box_x, box_y, box_z, box_m)
bodies = [pad, box]

# pad_joint = JointConstraint(Planar(
#     origin, pad, Y_AXIS, 
#     spring = spring,
#     damper = damper,
#     child_vertex = - Z_AXIS * (pad_r + box_z),
#     ))

pad_joint = JointConstraint(Prismatic(
    origin, pad, Dojo.X_AXIS, 
    spring = spring,
    damper = damper,
    child_vertex = - Dojo.Z_AXIS * (pad_r + box_z),
    ))

box_joint = JointConstraint(Prismatic(
    origin, box, Dojo.X_AXIS, 
    child_vertex = -Dojo.Z_AXIS * box_z / 2 .+ Dojo.X_AXIS * box_x / 2 .- 2.0 * Dojo.X_AXIS * pad_r))

joints = [pad_joint, box_joint]

collision = SphereBoxCollision{Float64, 2, 3, 6}(szeros(3), box_x, box_y, box_z, pad_r)
friction_parameterization = SA{Float64}[
    1.0  0.0
    0.0  1.0]
contact = NonlinearContact{Float64, 8}(μ, friction_parameterization, collision)
contacts = [ContactConstraint((contact, pad.id, box.id), name = :pad_belt)]

mech = Mechanism(origin, bodies, joints, contacts)

zero_coordinates!(mech)
zero_velocities!(mech)
# set_maximal_configurations!(pad, x = Z_AXIS * pad_r * 5.0)

u = zeros(input_dimension(mech))
u[1] = 0.0 #-100#-100#0.0 ## z axis force on pad
# u[2] = 0.0 ## x axis force on pad
# u[3] = 26.0 ## x axis force on box
u[2] = 380.0

K = 1000
storage = Storage(K, length(mech.bodies))

for k in 1:K
    Jx, Ju = get_minimal_gradients!(mech, get_minimal_state(mech), u)
    Dojo.save_to_storage!(mech, storage, k)
end

vis = Visualizer()
open(vis)
visualize(mech, storage, vis = vis)