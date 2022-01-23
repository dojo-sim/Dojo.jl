using Dojo
using MeshCat

vis = Visualizer()
open(vis)


# Parameters
ex = [1.;0.;0.]
ey = [0.;1.;0.]
ez = [0.;0.;1.]
h = 1.0
r = 0.1
vert11 = [0.;0.;h / 2]
vert12 = -vert11

# Links
origin = Origin()
bodies = [Box(3r, 2r, h, h, color = RGBA(1., 0., 0.)) for i = 1:2]

# Constraints
eqc1 = EqualityConstraint(Floating(origin, bodies[1], spring = 0.0, damper = 0.0))
eqc2 = EqualityConstraint(Revolute(bodies[1], bodies[2], ez, spring = 0.0, damper = 0.0))
eqcs = [eqc1, eqc2]

mech = Mechanism(origin, bodies, eqcs, g = g, Δt = Δt, spring=spring, damper=damper)

bthigh = geteqconstraint(mech, "bthigh")
eqcs[bthigh.id] = add_limits(mech, bthigh, rot_limits=[SVector{1}(joint_limits[1][1]), SVector{1}(joint_limits[2][1])])
eqc2_limited = add_limits(mech, bthigh, rot_limits=[SVector{1}(joint_limits[1][1]), SVector{1}(joint_limits[2][1])])
