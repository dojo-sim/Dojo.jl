using ConstrainedDynamics


# Parameters
joint_axis = [.0;1.0;1.0]
v1 = -[0.5;0;0.]
v2 = -[0.;0.5;0.]
v3 = -[0.;0;0.5]

length1 = 1.
width, depth = 0.1, 0.1

# Links
origin = Origin{Float64}()
link1 = Box(length1, width, depth, length1, color = RGBA(1., 1., 0.))
link2 = Box(width, length1, depth, length1, color = RGBA(1., 0., 0.))
link3 = Box(width, depth, length1, length1, color = RGBA(1., 0., 1.))

# Constraints
# joint1 = EqualityConstraint(Spherical(origin, link1; p2=v1))
# joint2 = EqualityConstraint(Spherical(origin, link2; p2=v2))
# joint3 = EqualityConstraint(Spherical(origin, link3; p2=v3))
joint1 = EqualityConstraint(Revolute(origin, link1, joint_axis; p2=v1))
joint2 = EqualityConstraint(Revolute(origin, link2, joint_axis; p2=v2))
joint3 = EqualityConstraint(Revolute(origin, link3, joint_axis; p2=v3))

links = [link1;link2;link3]
constraints = [joint1;joint2;joint3]


mech = Mechanism(origin, links, constraints, g = 0.)
setPosition!(origin,link1,Δq = UnitQuaternion(RotZ(0.)),p2=v1)
setPosition!(origin,link2,Δq = UnitQuaternion(RotX(0.)),p2=v2)
setPosition!(origin,link3,Δq = UnitQuaternion(RotX(0.)),p2=v3)
