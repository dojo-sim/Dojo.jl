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
vis=visualizer()
open(vis)

# Include new files
include(joinpath(module_dir(), "examples", "loader.jl"))


# Build mechanism
include("mechanism_zoo.jl")


################################################################################
# DEVELOPMENT NEW EQUALITY CONSTRAINT
################################################################################

# Parameters
joint_axis = [1.0;0.0;0.0]
length1 = 0.5
width, depth = 0.5, 0.5

origin = Origin{Float64}()
pbody = Box(width, depth, length1, 1., color = RGBA(1., 1., 0.))
joint0to1 = JointConstraint(Floating(origin, pbody, spring = 100.0, damper = 0.0))
bodies = [pbody]
eqcs = [joint0to1]
mech = Mechanism(origin, bodies, eqcs, g = 0.0, timestep=0.01)

pbody = mech.bodies[1]
set_position(pbody, x = [1, 0, 0.])
set_maximal_velocities!(pbody, v = [0, 0, 0.], ω = [0, 0, 0.])

storage = simulate!(mech, 5.0, record=true, solver = :mehrotra!)
visualize(mech, storage, vis=vis)

plot([x[1] for x in storage.x[1]])


################################################################################
# Differentiation
################################################################################

include(joinpath(module_dir(), "examples", "diff_tools.jl"))
# Set data
Nb = length(mech.bodies)
data = get_data(mech)
set_data!(mech, data)

mehrotra!(mech, opts = SolverOptions(rtol = 1e-6, btol = 1e-1, undercut=1.2, verbose=true))
sol = get_solution(mech)
attjac = attitude_jacobian(data, Nb)

# IFT
set_entries!(mech)
datamat = full_data_matrix(mech)
solmat = full_matrix(mech.system)
sensi = - (solmat \ datamat)
@show norm(full_vector(mech.system), Inf)

# finite diff
fd_datamat = finitediff_data_matrix(mech, data, sol) * attjac
@test norm(fd_datamat + datamat, Inf) < 1e-7
@test norm((fd_datamat + datamat)[:, 1:end-2], Inf) < 1e-7

plot(Gray.(abs.(datamat)))
plot(Gray.(abs.(fd_datamat)))

norm(fd_datamat + datamat, Inf)
norm((fd_datamat + datamat)[1:6,1:25], Inf)
norm((fd_datamat + datamat)[7:8,1:25], Inf)
norm((fd_datamat + datamat)[9:9,1:25], Inf)
fd_datamat[9:9,1:6]
fd_datamat[9:9,7:12]
fd_datamat[9:9,13:18]
fd_datamat[9:9,19:24]

datamat[9:9,1:6]
datamat[9:9,7:12]
datamat[9:9,13:18]
datamat[9:9,19:24]

(fd_datamat + datamat)[9:9,1:6] # pbody x2z
(fd_datamat + datamat)[9:9,7:12] # pbody q2x
(fd_datamat + datamat)[9:9,13:18] # cbody x2z
(fd_datamat + datamat)[9:9,19:24] # cbody q2x

norm((fd_datamat + datamat)[10:12,1:25], Inf)
norm((fd_datamat + datamat)[13:18,1:25], Inf)
norm((fd_datamat + datamat)[19:20,1:25], Inf)
norm((fd_datamat + datamat)[21:21,1:25], Inf)
norm((fd_datamat + datamat)[22:24,1:25], Inf)

fd_solmat = finitediff_sol_matrix(mech, data, sol)
@test norm(fd_solmat + solmat, Inf) < 1e-7
plot(Gray.(abs.(solmat)))
plot(Gray.(abs.(fd_solmat)))

fd_sensi = finitediff_sensitivity(mech, data) * attjac
@test norm(fd_sensi - sensi) / norm(fd_sensi) < 5e-3
plot(Gray.(sensi))
plot(Gray.(fd_sensi))
norm(fd_sensi - sensi, Inf)
norm(fd_sensi, Inf)

################################################################################
# Finite Diff
################################################################################
include(joinpath(@__DIR__, "finite_diff.jl"))

timestep=0.01
tra1 = mech.joints[1].constraints[1]
tra2 = mech.joints[2].constraints[1]
origin = mech.origin
pbody = mech.bodies[1]
cbody = mech.bodies[2]


jac0, jac1 = finitediff_vel(tra2, pbody, cbody, timestep, spring_parent, spring_parent_jacobian_velocity_parent, diff_body = :parent)
norm(jac0 - jac1, Inf)
jac0, jac1 = finitediff_vel(tra2, pbody, cbody, timestep, damper_parent, damper_parent_jacobian_velocity_parent, diff_body = :parent)
norm(jac0 - jac1, Inf)
jac0, jac1 = finitediff_vel(tra2, pbody, cbody, timestep, spring_parent, spring_parent_jacobian_velocity_child, diff_body = :child)
norm(jac0 - jac1, Inf)
jac0, jac1 = finitediff_vel(tra2, pbody, cbody, timestep, damper_parent, damper_parent_jacobian_velocity_child, diff_body = :child)
norm(jac0 - jac1, Inf)
jac0, jac1 = finitediff_vel(tra2, pbody, cbody, timestep, spring_child, spring_child_jacobian_velocity_child, diff_body = :child)
norm(jac0 - jac1, Inf)
jac0, jac1 = finitediff_vel(tra2, pbody, cbody, timestep, damper_child, damper_child_configuration_velocity_child, diff_body = :child)
norm(jac0 - jac1, Inf)
jac0, jac1 = finitediff_vel(tra2, pbody, cbody, timestep, spring_child, spring_child_configuration_velocity_parent, diff_body = :parent)
norm(jac0 - jac1, Inf)
jac0, jac1 = finitediff_vel(tra2, pbody, cbody, timestep, damper_child, damper_child_configuration_velocity_parent, diff_body = :parent)
norm(jac0 - jac1, Inf)
jac0, jac1 = finitediff_vel(tra1, origin, pbody, timestep, spring_child, spring_child_jacobian_velocity_child, diff_body = :child)
norm(jac0 - jac1, Inf)
jac0, jac1 = finitediff_vel(tra1, origin, pbody, timestep, damper_child, damper_child_configuration_velocity_child, diff_body = :child)
norm(jac0 - jac1, Inf)

jac0, jac1 = finitediff_pos(tra2, pbody, cbody, timestep, spring_parent, spring_parent_jacobian_configuration_parent, diff_body = :parent)
@test norm(jac0 - jac1, Inf) < 1e-8
jac0, jac1 = finitediff_pos(tra2, pbody, cbody, timestep, damper_parent, damper_parent_jacobian_configuration_parent, diff_body = :parent)
@test norm(jac0 - jac1, Inf) < 1e-8
jac0, jac1 = finitediff_pos(tra2, pbody, cbody, timestep, spring_parent, spring_parent_jacobian_configuration_child, diff_body = :child)
@test norm(jac0 - jac1, Inf) < 1e-8
jac0, jac1 = finitediff_pos(tra2, pbody, cbody, timestep, damper_parent, damper_parent_jacobian_configuration_child, diff_body = :child)
@test norm(jac0 - jac1, Inf) < 1e-8
jac0, jac1 = finitediff_pos(tra2, pbody, cbody, timestep, spring_child, spring_child_jacobian_configuration_child, diff_body = :child)
@test norm(jac0 - jac1, Inf) < 1e-8
jac0, jac1 = finitediff_pos(tra2, pbody, cbody, timestep, damper_child, damper_child_jacobian_configuration_child, diff_body = :child)
@test norm(jac0 - jac1, Inf) < 1e-8
jac0, jac1 = finitediff_pos(tra2, pbody, cbody, timestep, spring_child, spring_child_jacobian_configuration_parent, diff_body = :parent)
@test norm(jac0 - jac1, Inf) < 1e-8
jac0, jac1 = finitediff_pos(tra2, pbody, cbody, timestep, damper_child, damper_child_jacobian_configuration_parent, diff_body = :parent)
@test norm(jac0 - jac1, Inf) < 1e-8
jac0, jac1 = finitediff_pos(tra1, origin, pbody, timestep, spring_child, spring_child_jacobian_configuration_child, diff_body = :child)
@test norm(jac0 - jac1, Inf) < 1e-8
jac0, jac1 = finitediff_pos(tra1, origin, pbody, timestep, damper_child, damper_child_jacobian_configuration_child, diff_body = :child)
@test norm(jac0 - jac1, Inf) < 1e-8
