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

mech = get_mechanism(:pendulum, limits=true)
initialize!(mech, :pendulum, angle=0.1)
# Main.@profiler
storage = simulate!(mech, 1.0)
storage = simulate!(mech, 1.0)

visualize(mech, storage, vis=vis)
#


mech = get_mechanism(:quadruped)
initialize!(mech, :quadruped)
# Main.@profiler
storage = simulate!(mech, 1.0)
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









rotation_matrix0(q::Quaternion) = VRᵀmat(q) * LVᵀmat(q)
function ∂rotation_matrix∂q0(q::Quaternion, p::AbstractVector; attjac::Bool=false)
    M = ∂VRᵀmat∂q(LVᵀmat(q) * p) + VRᵀmat(q) * ∂LVᵀmat∂q(p)
    attjac && (M *= LVᵀmat(q))
    return M
end
function ∂rotation_matrix_inv∂q0(q::Quaternion, p::AbstractVector; attjac::Bool=false)
    M = ∂rotation_matrix∂q0(inv(q), p) * Tmat()
    attjac && (M *= LVᵀmat(q))
    return M
end

using FiniteDiff

p = srand(3)
q = Quaternion(normalize(srand(4))...)
J0 = FiniteDiff.finite_difference_jacobian(q -> rotation_matrix0(Quaternion(q...)) * p, vector(q)) * LVᵀmat(q)
J1 = ∂rotation_matrix∂q0(q, p, attjac=true)
J0 - J1

p = srand(3)
q = Quaternion(normalize(srand(4))...)
J0 = FiniteDiff.finite_difference_jacobian(q -> rotation_matrix0(inv(Quaternion(q...))) * p, vector(q)) * LVᵀmat(q)
J1 = ∂rotation_matrix_inv∂q0(q, p, attjac=true)
J0 - J1





p = srand(3)
q = Quaternion(normalize(srand(4))..., true)
inv(q)
# q = Quaternion(1,0,0,0.0)
# J0 = FiniteDiff.finite_difference_jacobian(q -> rotation_matrix(inv(Quaternion(q..., true))) * p, vector(q)) * LVᵀmat(q)
J0 = FiniteDiff.finite_difference_jacobian(q -> rotation_matrix(inv(Quaternion(q..., true))) * p, vector(q))# * LVᵀmat(q)
J1 = ∂rotation_matrix_inv∂q(q, p)
J0 - J1


rotation_matrix(q) = VRᵀmat(q) * LVᵀmat(q)
