using ConstrainedDynamics
using LinearAlgebra
using StaticArrays


# Parameters
h = .1
r = 1.

# Links
origin = Origin{Float64}()
link1 = Cylinder(r, h, h*r, color = RGBA(1., 0., 0.))

# Constraints
joint1 = EqualityConstraint(Floating(origin, link1))

links = [link1]
constraints = [joint1]


mech = Mechanism(origin, links, constraints, g = 0.)

axis = [0;0;1.]
speed = 2pi #*0
setVelocity!(link1, ω = speed*axis)

function nutation_control!(mechanism, k)
    if k==1
        setForce!(link1, F = SA[0;0;0.2] * 0, p = SA[0;1.;0])
    else
        setForce!(link1, F = szeros(3), p = szeros(3))
    end
    return
end
