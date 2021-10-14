using ConstrainedDynamics


# Parameters
joint_axis = [1.0;0.0;0.0]

length1 = 0.5
width, depth = 0.5, 0.5

# Links
origin = Origin{Float64}()
link1 = Box(width, depth, length1, 1., color = RGBA(1., 1., 0.))

# Constraints
joint1 = InequalityConstraint(Impact(link1, [0;-0.1;1.0]))
joint2 = InequalityConstraint(Impact(link1, [0.1;0;1.0]))
joint3 = InequalityConstraint(Impact(link1, [-0.1;0;1.0]))
joint4 = InequalityConstraint(Impact(link1, [0;0.1;1.0]))

links = [link1]
ineqs = [joint1;joint2;joint3;joint4]


mech = Mechanism(origin, links, ineqs)
setVelocity!(link1,v = [1;0.5;5])

