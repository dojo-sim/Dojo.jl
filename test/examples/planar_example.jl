using ConstrainedDynamics


# Parameters
ex = [1.;0.;0.]
ey = [0.;1.;0.]
ez = [0.;0.;1.]

l1 = 1.0
x, y = .1, .1
b1 = Box(x, l1, y, l1, color = RGBA(1., 1., 0.))

vert11 = [0.;l1 / 2;0.]
vert12 = -vert11

# Initial orientation
phi1, phi2, phi3 = 0., 2pi / 3, 4pi / 3
q1, q2, q3 = UnitQuaternion(RotZ(phi1)), UnitQuaternion(RotZ(phi2)), UnitQuaternion(RotZ(phi3))

# Links
origin = Origin{Float64}()
link1 = Box(x, l1, y, l1, color = RGBA(1., 1., 0.))
link2 = deepcopy(link1)

# Constraints
joint0to123 = EqualityConstraint(PlanarFree(origin, link1, ez; p2=vert12), PlanarFree(origin, link2, ez; p2=vert12))
joint1to23 = EqualityConstraint(Spherical(link1, link2; p1=vert11, p2=vert11))



links = [link1; link2]
constraints = [joint0to123;joint1to23]

mech = Mechanism(origin, links, constraints)
setPosition!(origin,link1,p2 = vert11,Δq = q1)
setPosition!(origin,link2,p2 = vert11,Δq = q2)