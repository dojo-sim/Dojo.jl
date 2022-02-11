using ForwardDiff
using FiniteDiff
using StaticArrays
using LinearAlgebra
using Random
using SparseArrays
using BenchmarkTools
using Dojo
using Test

include("integrator.jl")
include("minimal.jl")
include("jacobian.jl")
include("momentum.jl")
include("energy.jl")
include("behaviors.jl")
include("joint_limits.jl")
include("joint_impulse_map.jl")
include("data.jl")






mech = get_pendulum()
joint = mech.joints[1]
xa = rand(3)
qa = UnitQuaternion(rand(4)...)
xb = rand(3)
qb = UnitQuaternion(rand(4)...)
displacement_jacobian_configuration(:parent, joint.translational, xa, qa, xb, qb, attjac=true)
joint


a = 1
a = 1
a = 1
a = 1
a = 1
