using Pkg
# Pkg.develop(path=joinpath(@__DIR__, "../../DojoEnvironments"))
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
Main.@profiler storage = simulate!(mech, 100.0)
@btime storage = simulate!(mech, 10.0)

visualize(mech, storage, vis=vis)



function reset2!(joint::JointConstraint{T,N,Nc}; scale::T=1.0) where {T,N,Nc}
    Nλ_tra = joint_length(joint.translational)
    Nb_tra = limits_length(joint.translational)
    Nλ_rot = joint_length(joint.rotational)
    Nb_rot = limits_length(joint.rotational)
    joint.impulses[1] = [scale * sones(2Nb_tra); szeros(Nλ_tra); scale * sones(2Nb_rot); szeros(Nλ_rot)]
    joint.impulses[2] = [scale * sones(2Nb_tra); szeros(Nλ_tra); scale * sones(2Nb_rot); szeros(Nλ_rot)]
    return
end


joint = mech.joints[1]

reset!(joint)
Main.@code_warntype reset!(joint)
@benchmark $reset!($joint)



a = 10
a = 10
a = 10
a = 10
a = 10











function constraint_jacobian_configuration(mechanism::Mechanism{T,Nn,Ne,Nb}, body::Body{T}; reg::T=Dojo.REG) where {T,Nn,Ne,Nb}
    state = body.state
    timestep = mechanism.timestep
    mass = body.mass
    inertia = body.inertia

    x2, q2 = current_configuration(state)
    x3, q3 = next_configuration(state, timestep)

    I3 = SMatrix{3,3,T,9}(Diagonal(sones(T,3)))
    Z33 = szeros(T, 3, 3)
    Z34 = szeros(T, 3, 4)

    # dynamics
    dynT = I3 * mass / timestep
    dynR = -2.0 / timestep * LVᵀmat(q2)' * Tmat() * (∂Rᵀmat∂q(Vᵀmat() * inertia * Vmat() * Lmat(q2)' * vector(q3)) + Rmat(q3)' * Vᵀmat() * inertia * Vmat() * Lmat(q2)')

    state.D = [[dynT; Z33] [Z34; dynR]] * integrator_jacobian_velocity(body, timestep)
    state.D += [[reg * I3; Z33] [Z33; reg * I3]]

    # inputs
    nothing

    # impulses
    for id in connections(mechanism.system, body.id)
        Ne < id <= Ne + Nb && continue # body
        impulses_jacobian_velocity!(mechanism, body, get_node(mechanism, id))
    end

    return state.D
end

body = mech.bodies[1]
constraint_jacobian_configuration(mech, body)
Main.@code_warntype constraint_jacobian_configuration(mech, body)
@benchmark $constraint_jacobian_configuration($mech, $body)


a = 10
a = 10
a = 10
a = 10
a = 10
a = 10




















mech = get_mechanism(:quadruped)
initialize!(mech, :quadruped)
# Main.@profiler
storage = simulate!(mech, 1.0)
storage = simulate!(mech, 1.0)

visualize(mech, storage, vis=vis)
#


using BenchmarkTools

# function constraint_jacobian_configuration(relative::Symbol, joint::Joint{T,Nλ,Nb},
#         xa::AbstractVector, qa::Quaternion,
#         xb::AbstractVector, qb::Quaternion,
#         η) where {T,Nλ,Nb}
#
#     ∇comp = szeros(T,Nb,7)
#     ∇mincoord = minimal_coordinates_jacobian_configuration(relative, joint, xa, qa, xb, qb, attjac=false)
#     ∇unlim = joint_constraint_jacobian_configuration(relative, joint, xa, qa, xb, qb, η)
#
#     return [∇comp; ∇mincoord; -∇mincoord; ∇unlim]
# end

relative = :parent
joint = mech.joints[5].rotational
xa = srand(3)
xb = srand(3)
qa = Quaternion(normalize(rand(4))...)
qb = Quaternion(normalize(rand(4))...)
η = 0.0

constraint_jacobian_configuration(relative, joint, xa, qa, xb, qb, η)
@benchmark $constraint_jacobian_configuration($relative, $joint, $xa, $qa, $xb, $qb, $η)


a = 10
a = 10
a = 10
a = 10













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
