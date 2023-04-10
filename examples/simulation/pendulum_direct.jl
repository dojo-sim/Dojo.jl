using Dojo

radius = 0.1
length = 1
mass = 1
axis = [0;0;1]
connection = [0;0;length/2]

origin = Origin{T}()
body = Cylinder(length, radius, mass)
joint = JointConstraint(Revolute(origin, body, axis; child_vertex=connection))

mechanism = Mechanism(origin, [body], [joint])

set_minimal_coordinates!(mechanism, joint, [pi/4])

storage = simulate!(mechanism, 5.0, record=true)

vis = visualize(mechanism, storage; visualize_floor=false)
render(vis)