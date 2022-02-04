set_position
# set_velocity
minimal_velocities
minimal_coordinates
# minimal_to_maximal
# maximal_to_minimal
# set_joint_position!

mech = get_atlas()
for j in mech.joints
    @show vector(j.constraints[2].qoffset)
end

mech = get_humanoid()
for j in mech.joints
    @show vector(j.constraints[2].qoffset), j.name
end


include(joinpath(module_dir(), "src", "joints", "rotational", "minimal.jl"))
include(joinpath(module_dir(), "src", "joints", "translational", "minimal.jl"))

include(joinpath(module_dir(), "src", "joints", "rotational", "minimal_dev.jl"))
include(joinpath(module_dir(), "src", "joints", "translational", "minimal_dev.jl"))
################################################################################
# Test humanoid
################################################################################
mech = Dojo.get_mechanism(:humanoid)
Random.seed!(100)
nx = Dojo.minimal_dimension(mech)
x0 = [rand(3); Dojo.vector(UnitQuaternion(rand(4)...)); rand(3); rand(3); rand(abs(nx - 13))]
z0 = Dojo.minimal_to_maximal(mech, x0)
x1 = Dojo.maximal_to_minimal(mech, z0)
@test norm(x0[1:3] - x1[1:3], Inf) < 1e-10
@test quaterror(x0[4:7], x1[4:7]) < 1e-10
@test norm(x0[8:10] - x1[8:10], Inf) < 1e-10
@test norm(x0[11:nx] - x1[11:nx], Inf) < 1e-10


mech = get_snake(Nb=2)
Random.seed!(100)
nx = Dojo.minimal_dimension(mech)
x0 = [rand(3); Dojo.vector(UnitQuaternion(rand(4)...)); rand(3); rand(3); rand(abs(nx - 13))]
z0 = Dojo.minimal_to_maximal(mech, x0)
x1 = Dojo.maximal_to_minimal(mech, z0)
@test norm(x0[1:3] - x1[1:3], Inf) < 1e-10
@test quaterror(x0[4:7], x1[4:7]) < 1e-10
@test norm(x0[8:10] - x1[8:10], Inf) < 1e-10
@test norm(x0[11:nx] - x1[11:nx], Inf) < 1e-10

# vis = Visualizer()
# open(vis)


################################################################################
# qoffset = 0, zeros
################################################################################
mech = get_snake(Nb=1)
initialize!(mech, :snake, x=[0,0,0.], v=[0,0,0.], q1=one(UnitQuaternion), ω=[0,0,0.])
storage = generate_storage(mech, [get_maximal_state(mech)])
visualize(mech, storage, vis=vis)

z1 = [zeros(3); zeros(3); 1;0;0;0; zeros(3)]
x1 = [zeros(3); zeros(3); zeros(3); zeros(3)]
z0 = get_maximal_state(mech)
x0 = get_minimal_state(mech)
norm(z0 - z1, Inf) < 1e-10
norm(x0 - x1, Inf) < 1e-10

################################################################################
# qoffset = 0, Translational
################################################################################
mech = get_snake(Nb=1)
xi = rand(3)
vi = rand(3)
initialize!(mech, :snake, x=xi, v=vi, q1=one(UnitQuaternion), ω=[0,0,0.])
storage = generate_storage(mech, [get_maximal_state(mech)])
visualize(mech, storage, vis=vis)

z1 = [-xi; vi; 1;0;0;0; zeros(3)]
x1 = [-xi; zeros(3); vi; zeros(3)]
z0 = get_maximal_state(mech)
x0 = get_minimal_state(mech)
norm(z0 - z1, Inf) < 1e-10
norm(x0 - x1, Inf) < 1e-10

################################################################################
# qoffset = 0, Rotational
################################################################################
mech = get_snake(Nb=1)
qi = UnitQuaternion(rand(4)...)
ωi = rand(3)
initialize!(mech, :snake, x=zeros(3), v=zeros(3), q1=qi, ω=ωi)
storage = generate_storage(mech, [get_maximal_state(mech)])
visualize(mech, storage, vis=vis)

z1 = [zeros(3); zeros(3); vector(qi); ωi]
x1 = [zeros(3); rotation_vector(qi); zeros(3); rotation_matrix(qi)*ωi]
z0 = get_maximal_state(mech)
x0 = get_minimal_state(mech)
norm(z0 - z1, Inf) < 1e-10
norm(x0 - x1, Inf) < 1e-10
norm(x0 - x1, Inf)

################################################################################
# qoffset = 0, Translatioanl & Rotational
################################################################################
mech = get_snake(Nb=1)
xi = [0,1,0.]#zeros(3)#rand(3)
# vi = rand(3)
vi = [1,2,3.]
# qi = UnitQuaternion(rand(4)...)
qi = UnitQuaternion(1,1,0,0.)
ωi = [3,4,5.]
initialize!(mech, :snake, x=xi, v=vi, q1=qi, ω=ωi)
storage = generate_storage(mech, [get_maximal_state(mech)])
visualize(mech, storage, vis=vis)

z1 = [-xi; vi; vector(qi); ωi]
x1 = [-xi; rotation_vector(qi); vi; rotation_matrix(qi)*ωi]
z0 = get_maximal_state(mech)
x0 = get_minimal_state(mech)
norm(z0 - z1, Inf) < 1e-10
norm(x0 - x1, Inf) < 1e-10
norm(x0 - x1, Inf)

scatter(z0 - z1)










storage = generate_storage(mech, [get_maximal_state(mech)])
visualize(mech, storage, vis=vis)

z1 = [zeros(3); zeros(3); 1;0;0;0; zeros(3)]
x1 = [zeros(3); zeros(3); zeros(3); zeros(3)]
z0 = get_maximal_state(mech)
x0 = get_minimal_state(mech)
norm(z0 - z1, Inf) < 1e-10
norm(x0 - x1, Inf) < 1e-10




joint0 = mech.joints[1]
tra0 = joint0.constraints[1]
rot0 = joint0.constraints[2]
nullspace_mask(tra0)
nullspace_mask(rot0)
joint0.constraints[2].qoffset

q0 = UnitQuaternion(1,0,0,0.)
rotation_vector(q0)

origin0 = mech.origin
body0 = mech.bodies[1]
minimal_coordinates(tra0, origin0, body0)
minimal_coordinates(rot0, origin0, body0)














qa = UnitQuaternion(rand(4)...)
qb = UnitQuaternion(rand(4)...)
qoffset = UnitQuaternion(rand(4)...)
vector(qa \ qb) - vector(inv(qa) * qb)
vector(qb / qoffset) - vector(qb * inv(qoffset))
vector(qa \ qb / qoffset) - vector(inv(qa) * qb * inv(qoffset))









nullspace_mask(mech.joints[1].constraints[1])
constraint_mask(mech.joints[1].constraints[1])

mech = get_slider()

mech.joints[1].constraints[1].vertices
mech.joints[1].minimal_index
