include("examples/doublependulum_3d.jl")

using ConstrainedDynamics: currentasknot!, getbody, geteqconstraint, getcomponent

dq = minimalCoordinates(mech)
@test true
setPosition!(mech,dq)
@test true

dF = ConstrainedDynamics.UnitDict(dq.keys,[[0.0;0;0] for i=1:2])
setForce!(mech,dF)
@test true

currentasknot!(mech)
@test true

body1 = getbody(mech, "body1")
body2 = getbody(mech, "body2")
eqc1 = geteqconstraint(mech, "jointb1")
eqc2 = geteqconstraint(mech, "joint12")
@test body1 === mech.bodies[3]
@test body2 === mech.bodies[4]
@test eqc1 === mech.eqconstraints[1]
@test eqc2 === mech.eqconstraints[2]

body1 = getcomponent(mech, "body1")
body2 = getcomponent(mech, "body2")
eqc1 = getcomponent(mech, "jointb1")
eqc2 = getcomponent(mech, "joint12")
@test body1 === mech.bodies[3]
@test body2 === mech.bodies[4]
@test eqc1 === mech.eqconstraints[1]
@test eqc2 === mech.eqconstraints[2]