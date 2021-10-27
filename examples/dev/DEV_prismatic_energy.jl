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
include(joinpath(module_dir(), "examples", "dev", "loader.jl"))


# Build mechanism
include("mechanism_zoo.jl")


################################################################################
# DEVELOPMENT NEW EQUALITY CONSTRAINT
################################################################################

# Parameters
ex = [0; 0; 1.0]
h = 1.
r = .05
vert11 = [0; r; 0.0]
vert12 = -vert11
Nlink = 2

# Links
origin = Origin{Float64}()
links = [Cylinder(r, h, h, color = RGBA(1., 0., 0.)) for i = 1:Nlink]

# Constraints
jointb1 = EqualityConstraint(Fixed(origin, links[1]; p1 = zeros(3), p2 = zeros(3)))
if Nlink > 1
    eqcs = [
        jointb1;
        [EqualityConstraint(Prismatic(links[i - 1], links[i], ex; p1=vert12, p2=vert11, spring = 100.0, damper = 0.0)) for i = 2:Nlink]
        ]
else
    eqcs = [jointb1]
end
mech = Mechanism(origin, links, eqcs, g = -9.81, Δt = 0.1)

# mech = getmechanism(:nslider, Nlink = 5)
initialize!(mech, :nslider)
storage = simulate!(mech, 5.0, record = true, solver = :mehrotra!)
visualize(mech, storage, vis = vis)

################################################################################
# Differentiation
################################################################################

include(joinpath(module_dir(), "examples", "dev", "diff_tools.jl"))
# Set data
Nb = length(mech.bodies)
data = getdata(mech)
setdata!(mech, data)

mehrotra!(mech, opts = InteriorPointOptions(rtol = 1e-6, btol = 1e-1, undercut=1.2, verbose=true))
sol = getsolution(mech)
attjac = attitudejacobian(data, Nb)

# IFT
setentries!(mech)
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

(fd_datamat + datamat)[9:9,1:6] # body1 x2z
(fd_datamat + datamat)[9:9,7:12] # body1 q2x
(fd_datamat + datamat)[9:9,13:18] # body2 x2z
(fd_datamat + datamat)[9:9,19:24] # body2 q2x

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

Δt = 0.01
tra1 = mech.eqconstraints[1].constraints[1]
tra2 = mech.eqconstraints[2].constraints[1]
origin = mech.origin
body1 = collect(mech.bodies)[1]
body2 = collect(mech.bodies)[2]


jac0, jac1 = finitediff_vel(tra2, body1, body2, Δt, springforcea, ∂springforcea∂vela, diff_body = :parent)
norm(jac0 - jac1, Inf)
jac0, jac1 = finitediff_vel(tra2, body1, body2, Δt, damperforcea, ∂damperforcea∂vela, diff_body = :parent)
norm(jac0 - jac1, Inf)
jac0, jac1 = finitediff_vel(tra2, body1, body2, Δt, springforcea, ∂springforcea∂velb, diff_body = :child)
norm(jac0 - jac1, Inf)
jac0, jac1 = finitediff_vel(tra2, body1, body2, Δt, damperforcea, ∂damperforcea∂velb, diff_body = :child)
norm(jac0 - jac1, Inf)
jac0, jac1 = finitediff_vel(tra2, body1, body2, Δt, springforceb, ∂springforceb∂velb, diff_body = :child)
norm(jac0 - jac1, Inf)
jac0, jac1 = finitediff_vel(tra2, body1, body2, Δt, damperforceb, ∂damperforceb∂velb, diff_body = :child)
norm(jac0 - jac1, Inf)
jac0, jac1 = finitediff_vel(tra2, body1, body2, Δt, springforceb, ∂springforceb∂vela, diff_body = :parent)
norm(jac0 - jac1, Inf)
jac0, jac1 = finitediff_vel(tra2, body1, body2, Δt, damperforceb, ∂damperforceb∂vela, diff_body = :parent)
norm(jac0 - jac1, Inf)
jac0, jac1 = finitediff_vel(tra1, origin, body1, Δt, springforceb, ∂springforceb∂velb, diff_body = :child)
norm(jac0 - jac1, Inf)
jac0, jac1 = finitediff_vel(tra1, origin, body1, Δt, damperforceb, ∂damperforceb∂velb, diff_body = :child)
norm(jac0 - jac1, Inf)

jac0, jac1 = finitediff_pos(tra2, body1, body2, Δt, springforcea, ∂springforcea∂posa, diff_body = :parent)
@test norm(jac0 - jac1, Inf) < 1e-8
jac0, jac1 = finitediff_pos(tra2, body1, body2, Δt, damperforcea, ∂damperforcea∂posa, diff_body = :parent)
@test norm(jac0 - jac1, Inf) < 1e-8
jac0, jac1 = finitediff_pos(tra2, body1, body2, Δt, springforcea, ∂springforcea∂posb, diff_body = :child)
@test norm(jac0 - jac1, Inf) < 1e-8
jac0, jac1 = finitediff_pos(tra2, body1, body2, Δt, damperforcea, ∂damperforcea∂posb, diff_body = :child)
@test norm(jac0 - jac1, Inf) < 1e-8
jac0, jac1 = finitediff_pos(tra2, body1, body2, Δt, springforceb, ∂springforceb∂posb, diff_body = :child)
@test norm(jac0 - jac1, Inf) < 1e-8
jac0, jac1 = finitediff_pos(tra2, body1, body2, Δt, damperforceb, ∂damperforceb∂posb, diff_body = :child)
@test norm(jac0 - jac1, Inf) < 1e-8
jac0, jac1 = finitediff_pos(tra2, body1, body2, Δt, springforceb, ∂springforceb∂posa, diff_body = :parent)
@test norm(jac0 - jac1, Inf) < 1e-8
jac0, jac1 = finitediff_pos(tra2, body1, body2, Δt, damperforceb, ∂damperforceb∂posa, diff_body = :parent)
@test norm(jac0 - jac1, Inf) < 1e-8
jac0, jac1 = finitediff_pos(tra1, origin, body1, Δt, springforceb, ∂springforceb∂posb, diff_body = :child)
@test norm(jac0 - jac1, Inf) < 1e-8
jac0, jac1 = finitediff_pos(tra1, origin, body1, Δt, damperforceb, ∂damperforceb∂posb, diff_body = :child)
@test norm(jac0 - jac1, Inf) < 1e-8
