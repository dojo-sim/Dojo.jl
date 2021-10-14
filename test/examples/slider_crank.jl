using StaticArrays
using ConstrainedDynamics


# Parameters
ex = [1.;0.;0.]
ey = [0.;1.;0.]
ez = [0.;0.;1.]

l1 = 2.0
x, y = .1, .1
r = l1 / 4
h = r / 10
box = Box(x, l1, y, l1, color = RGBA(1., 1., 0.))
cyl = Cylinder(r, h, 2 * r, color = RGBA(1., 0., 0.))

p0 = [0;l1 / 2;0]
p1 = -p0
p2 = [0;r;0]
p3 = -p2

# Initial orientation
q1 = UnitQuaternion(RotY(pi / 2))

# Links
origin = Origin{Float64}()
link1 = Box(x, l1, y, l1, color = RGBA(1., 1., 0.))
link2 = Cylinder(r, h, 2 * r, color = RGBA(1., 0., 0.))

# Constraints
joint0to12 = EqualityConstraint(CylindricalFree(origin, link1, ey; p2=p0), Revolute(origin, link2, ez; p1=p3 + p1))
joint1to2 = EqualityConstraint(Cylindrical(link1, link2, ez; p1=p1, p2=p3))


links = [link1; link2]
constraints = [joint0to12;joint1to2]

function slider_crank_control!(mechanism, k)
    τ = SA[0;0.3]
    setForce!(mechanism, joint1to2, τ * 0)
    return
end


mech = Mechanism(origin, links, constraints)
setPosition!(link1,x = -p0)
setPosition!(link2,x = p3 + p1)
