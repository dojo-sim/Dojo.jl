using ConstrainedDynamics
using LinearAlgebra
using StaticArrays

# Parameters
h = .1
r = 1.

# length1 = 0.1
# width, depth = 2., 2.
# b1 = Box(width, depth, length1, length1, color = RGBA(1., 1., 0.))

# Links
origin = Origin{Float64}()
link1 = Cylinder(r, h, h*r, color = RGBA(1., 0., 0.))

# Constraints
joint1 = EqualityConstraint(Floating(origin, link1))

links = [link1]
constraints = [joint1]


mech = Mechanism(origin, links, constraints, g = 0., Δt=0.005)

axis = [0;0;1.]
speed = 50pi #*0
setVelocity!(link1, ω = speed*axis)

function control!(mechanism, k)
    if k==1
        setForce!(link1, F = SA[0;0;2.], p = SA[0;1.;0])
    else
        setForce!(link1, F = szeros(3), p = szeros(3))
    end
    return
end

storage = simulate!(mech, 10., control!, record = true)
# visualize(mech, storage)

p = [storage.x[1][i] + ConstrainedDynamics.vrotate([0;1.0;0],storage.q[1][i]) for i=1:1000]
@test maximum(getindex.(p-[[0;0;(i-1)*0.002] for i=1:1000],3)) < 0.01