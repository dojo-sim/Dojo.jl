using Pkg
Pkg.develop(path=joinpath(@__DIR__, "../../DojoEnvironments"))
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

# ## visualizer
vis = Visualizer()
open(vis)

# ## setup
using Dojo
using IterativeLQR
using LinearAlgebra
using FiniteDiff
using DojoEnvironments

using BenchmarkTools

mech = get_mechanism(:pendulum)
initialize!(mech, :pendulum, angle=0.1)
Main.@profiler storage = simulate!(mech, 10.0)
storage = simulate!(mech, 1.0)

visualize(mech, storage, vis=vis)
#
# function joint_impulse_index(joint::JointConstraint{T,N,Nc}, i::Int) where {T,N,Nc}
#     s = 0
#     for j = 1:i-1
#         element = (joint.translational, joint.rotational)[j]
#         s += impulses_length(element)
#     end
#     joint_impulse_index((joint.translational, joint.rotational)[i], s) # to be allocation free
#     # joint_impulse_index((joint.translational, joint.rotational)[i], s) # to be allocation free
#     # joint_impulse_index((joint.translational, joint.rotational)[i], s) # to be allocation free
# end
#
# function joint_impulse_index(joint::Joint{T,Nλ,Nb,N}, s::Int) where {T,Nλ,Nb,N}
#     ind = SVector{N,Int}(s+1:s+N)
#     return ind
# end


joint_con0 = mech.joints[1]
joint0 = joint_con0.translational
joint0 = joint_con0.rotational
i0 = 1
s0 = 1
@benchmark $joint_impulse_index($joint_con0, $i0)
@benchmark $joint_impulse_index($joint0, $s0)


for (i,el) in enumerate(3,4,5)
    @show el
end



rosrun unitree_controller unitree_servo # let the robot stretch legs
rosrun unitree_controller unitree_move_kinetic # place the robot back to origin

roslaunch unitree_gazebo normal.launch rname:=a1 wname:=stairs_new
