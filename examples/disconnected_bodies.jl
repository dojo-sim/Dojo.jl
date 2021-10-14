using ConstrainedDynamics
using ConstrainedDynamicsVis


# Parameters
joint_axis = [1.0;0.0;0.0]

length1 = 0.5
width, depth = 0.5, 0.5

# Links
origin = Origin{Float64}()
link1 = Box(width, depth, length1, length1, color = RGBA(1., 1., 0.))
link2 = Box(width, depth, length1, length1, color = RGBA(1., 0., 0.))

links = [link1;link2]


mech = Mechanism(origin, links, g = -5.)
setPosition!(link1,x = [0.;1.;0.])
setPosition!(link2,x = [0.;-1.;0.])
setVelocity!(link1,v = [-5.;-5;5])
setVelocity!(link2,v = [-2.5;2.5;5],ω = [2;4;10.])


storage = simulate!(mech, 10., record = true)
visualize(mech, storage)
