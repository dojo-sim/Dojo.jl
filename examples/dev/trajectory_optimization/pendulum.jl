# Utils
function module_dir()
    return joinpath(@__DIR__, "..", "..")
end

# Activate package
using Pkg
Pkg.activate(module_dir())

# Load packages
using Plots
using Random
using MeshCat

# Open visualizer
vis = Visualizer()
open(vis)

# Include new files
include(joinpath(module_dir(), "loader.jl"))

# System 
gravity = -9.81 
Δt = 0.1

# Parameters
pendulum_axis = [1.0; 0.0; 0.0]
pendulum_length = 1.0
width, depth, height = 0.1, 0.1, 0.1
pendulum_mass = 1.0 

# Float type 
T = Float64 

# Links
origin = Origin{T}()
pendulum = Box(width, depth, pendulum_length, pendulum_mass)
links = [pendulum]

# Joint Constraints
joint_slider_pendulum = EqualityConstraint(Revolute(origin, pendulum, pendulum_axis; p1=zeros(3), p2=[0.0; 0.0; 0.5 * pendulum_length]))
eqcs = [joint_slider_pendulum]

# Mechanism
mech = Mechanism(origin, links, eqcs, g=gravity, Δt=Δt)

# origin to pendulum
setPosition!(mech.origin, mech.bodies[2], Δx=[0.0; 0.0; -0.5 * pendulum_length])
setVelocity!(mech.bodies[2], v = [0.0; 0.0; 0.0], ω = zeros(3))
mech.bodies[2].state.xc 

setPosition!(mech.origin, mech.bodies[2], Δx=[0.0; 0.0; 0.5 * pendulum_length], Δq=UnitQuaternion(RotX(π)))
setVelocity!(mech.bodies[2], v = zeros(3), ω = zeros(3))

# controller
function controller!(mech, k)
    j1 = geteqconstraint(mech, mech.eqconstraints[1].id)

    u1 = 0.1

    setForce!(mech, j1, SA[u1])

    return
end 

# simulate
storage = simulate!(mech, 5 * mech.Δt, controller!, record = true, solver = :mehrotra!)

plot(hcat(storage.x[1]...)')

# visualize
visualize(mech, storage, vis = vis)

## state space 
n = 13 * length(mech.bodies)
m = isempty(mech.eqconstraints) ? 0 : sum(getcontroldim.(mech.eqconstraints))

function step!(mech::Mechanism, x1::Vector{T}, u1::Vector{T}) 
    # set data
    data = [x1; u1] 
    setdata!(mech, data) 

    function controller!(mech, k)
        j1 = geteqconstraint(mech, mech.eqconstraints[1].id)
        setForce!(mech, j1, SA[u1[1]])
        return
    end 

    # simulate  
    simulate!(mech, mech.Δt, 
        controller!, 
        record=true, solver=:mehrotra!)

    # next state
    x2 = getnextstate(mech) 

    return x2
end

# initial state 
x1 = [0.0; 0.0; -0.5 * pendulum_length] 
v1 = [0.0; 0.0; 0.0] 
q1 = [1.0; 0.0; 0.0; 0.0]
ω1 = [0.0; 0.0; 0.0] 
z1 = [x1; v1; q1; ω1]

# target state 
xT = [0.0; 0.0; 0.5 * pendulum_length]
vT = [0.0; 0.0; 0.0] 
qT = [0.0; 1.0; 0.0; 0.0]
ωT = [0.0; 0.0; 0.0]
zT = [xT; vT; qT; ωT]

z = [copy(z1)]
for t = 1:5 
    znext = step!(mech, z[end], [0.1]) 
    push!(z, znext)
end

plot(hcat(z...)[1:3, :]')

