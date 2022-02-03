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
