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
include(joinpath(module_dir(), "examples", "loader.jl"))

# System 
gravity = -9.81 
Δt = 0.1

# Parameters
slider_axis = [0.0; 1.0; 0.0]
pendulum_axis = [1.0; 0.0; 0.0]
slider_length = 1.0
pendulum_length = 1.0
width, depth, height = 0.1, 0.1, 0.1
slider_mass = 1.0 
pendulum_mass = 1.0 

# Float type 
T = Float64 

# Links
origin = Origin{T}()
slider = Box(width, slider_length, height, slider_mass)
pendulum = Box(width, depth, pendulum_length, pendulum_mass)
links = [slider, pendulum]

# Joint Constraints
joint_origin_slider = EqualityConstraint(Prismatic(origin, slider, slider_axis; p1=szeros(T, 3), p2=szeros(T, 3)))
joint_slider_pendulum = EqualityConstraint(Revolute(slider, pendulum, pendulum_axis; p1=szeros(T, 3), p2=[0.0; 0.0; 0.5 * pendulum_length]))
eqcs = [joint_origin_slider, joint_slider_pendulum]

# Mechanism
mech = Mechanism(origin, links, eqcs, g=gravity, Δt=Δt)

# origin to slider
setPosition!(mech.origin, mech.bodies[3])
setVelocity!(mech.bodies[3], v = [0.0; 0.0; 0.0], ω = zeros(3))

# slider to pendulum
setPosition!(mech.bodies[3], mech.bodies[4], Δx=[0.0; 0.0; -0.5 * pendulum_length])
setVelocity!(mech.bodies[4], v = zeros(3), ω = zeros(3))

setPosition!(mech.bodies[3], mech.bodies[4], Δx=[0.0; 0.0; 0.5 * pendulum_length], Δq=UnitQuaternion(RotX(π)))
setVelocity!(mech.bodies[4], v = zeros(3), ω = zeros(3))

# controller 
function controller!(mech, k)
    j1 = geteqconstraint(mech, mech.eqconstraints[1].id)
    j2 = geteqconstraint(mech, mech.eqconstraints[2].id)

    u1 = 0.2
    u2 = 0.0

    setForce!(mech, j1, SA[u1])
    setForce!(mech, j2, SA[u2])

    return
end 

mech.bodies[3].state

mech.eqconstraints[1].constraints[1].Fτ

# simulate
storage = simulate!(mech, Δt, controller!, record = true, solver = :mehrotra!)

# visualize
visualize(mech, storage, vis = vis)

## state space 
n = 13 * 2 
m = isempty(mech.eqconstraints) ? 0 : sum(getcontroldim.(mech.eqconstraints))

function step!(mech::Mechanism, x1::Vector{T}, u1::Vector{T}) 
    # set data
    data = [x1; u1] 
    setdata!(mech, data) 

    # simulate  
    simulate!(mech, mech.Δt, record=false, solver=:mehrotra!)

    # next state
    x2 = getnextstate(mech) 

    return x2
end






