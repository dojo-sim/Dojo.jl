using ConstrainedDynamics
using ConstrainedDynamicsVis

# Initial orientation
ϕ1 = 0;
q1 = UnitQuaternion(RotX(ϕ1))

# Links
origin = Origin{Float64}()
link1 = Cylinder(0.5,0.1, 1., color = RGBA(1., 1., 0.))

# Constraints
N = 5

angle = [i*pi/N for i=1:N]
corners1 = [0.5*[sin(angle[i]);cos(angle[i]);0] for i=1:N]
corners2 = [0.5*[sin(angle[i]+pi);cos(angle[i]+pi);0] for i=1:N]
corners = [corners1;corners2]

fricsandineqs = [Friction(link1, [0;0;1.0], 0.2; p = corners[i]) for i=1:2*N]
frics = getindex.(fricsandineqs,1)
ineqcs = vcat(getindex.(fricsandineqs,2)...)


joint0to1 = EqualityConstraint(Floating(origin, link1))

links = [link1]
eqcs = [joint0to1]

mech = Mechanism(origin, links, eqcs, ineqcs, frics)

setPosition!(link1,x = [0.;0;0.0])
storage = simulate!(mech, 5, record = true, debug=false, ε=1e-3)
visualize(mech, storage)



