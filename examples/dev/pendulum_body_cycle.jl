# Parameters
joint_axis = [1.0; 0; 0]
length1 = 1.0
width, depth = 0.1, 0.1
p2 = [0; 0; length1/2] # joint connection point

# Links
origin = Origin{T}()
link1 = Box(width, depth, length1, length1)

# Constraints
joint_between_origin_and_link1 = EqualityConstraint(Revolute(origin, link1, joint_axis; p2=p2))
links = [link1]
eqcs = [joint_between_origin_and_link1]

links[1].id

mech = Mechanism(origin, links, eqcs, g = g, Δt = Δt)
return mech