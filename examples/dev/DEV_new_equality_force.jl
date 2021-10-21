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


###############################################################################
# METHODS NEW EQUALITY CONSTRAINT
################################################################################

# t2r3
function ForcePrismatic(body1::AbstractBody{T}, body2, axis;
        p1 = szeros(T, 3), p2 = szeros(T, 3), qoffset = one(UnitQuaternion{T}),
        spring = zero(T), damper = zero(T)) where T
    Translational2{T}(body1, body2; p1, p2, axis, spring, damper),
    Force121{T}(body1, body2; p1, p2, axis, spring, damper),
    Rotational3{T}(body1, body2; qoffset, spring, damper)
end

################################################################################
# DEVELOPMENT NEW EQUALITY CONSTRAINT
################################################################################

# Parameters
ex = [0; 0; 1.0]
h = 1.
r = .05
vert11 = [0; r; 0.0]
vert12 = -vert11
Nlink = 3

# Links
origin = Origin{Float64}()
links = [Cylinder(r, h, h, color = RGBA(1., 0., 0.)) for i = 1:Nlink]

# Constraints
jointb1 = EqualityConstraint(Fixed(origin, links[1]; p1 = zeros(3), p2 = zeros(3)))
if Nlink > 1
    eqcs = [
        jointb1;
        [EqualityConstraint(ForcePrismatic(links[i - 1], links[i], ex; p1=vert12, p2=vert11, spring = 1.0, damper = 1.0)) for i = 2:Nlink]
        ]
else
    eqcs = [jointb1]
end
mech = Mechanism(origin, links, eqcs, g = -9.81, Δt = 0.01)

# mech = getmechanism(:nslider, Nlink = 5)
initialize!(mech, :nslider)
# storage = simulate!(mech, 1.0, record = true, solver = :mehrotra!)
storage = simulate!(mech, 0.1, record = true, solver = :mehrotra!)

# visualize(mech, storage, vis = vis)

forc = eqcs[2].constraints[2]
xb = mech.bodies[4].state.xc
qb = mech.bodies[4].state.qc
ωb = mech.bodies[4].state.ωc
Δt = mech.Δt
∂g∂posb(forc, xb, qb)
fg = q -> ∂g∂posb(forc, xb, UnitQuaternion(q...))[2] * LVᵀmat(UnitQuaternion(q...))
fg(qb)
ForwardDiff.jacobian(fg, [qb.w; qb.x; qb.y; qb.z])


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




eqc1 = mech.eqconstraints[1]
eqc2 = mech.eqconstraints[2]
tra2 = eqc2.constraints[1]
force2 = eqc2.constraints[2]
rot2 = eqc2.constraints[3]
body1, body2 = mech.bodies


∂Fτ∂ua(mech, eqc1, getbody(mech, eqc2.parentid))
∂Fτ∂ua(mech, eqc2, getbody(mech, eqc2.parentid))
∂Fτ∂ua(force2, body1, body2, body2.id)

eqc1.λinds
∂g∂ʳposa(tra1, body1, body2, body2.id)


constraintmat(tra1)
A = constraintmat(force1)
Aᵀ = zerodimstaticadjoint(A)
Aᵀ * A * ([0.1] * A)

([0.1] * A)
∂g∂ʳvelb(mech, eqc1, body1)
∂g∂ʳvelb(tra1, body1, body2, body2.id, mech.Δt)
∂g∂ʳvelb(force1, body1, body2, body2.id, mech.Δt)
∂g∂ʳvelb(rot1, body1, body2, body2.id, mech.Δt)




∂g∂ʳself(mech, eqc1)
∂g∂ʳself(tra1)
∂g∂ʳself(force1)
∂g∂ʳself(rot1)
g(mech, eqc1)
g(mech, eqc2)
gc(mech, eqc2)

constraintmat(force1)
g(tra1, origin, body1, mech.Δt, rand(33))
g(rot1, origin, body1, mech.Δt, rand(33))
g(force1, origin, body1)
g(force1, body1, body2)


λindsf1 = Vector(eqc1.λinds[2][1]:eqc1.λinds[2][2])
eqc1.λsol[2][λindsf1]
g(force1, origin, body1)
g(force1, body1, body2)

springforce(force1, body1.state)
springforce(force1, body1.state, body2.state)

damperforce(force1, body1.state)
damperforce(force1, body1.state, body2.state)


eqc1.indss
eqc2.inds

Vector(SUnitRange(eqc1.inds[2][1], eqc1.inds[2][2]))



tra, ftr, rot = ForcePrismatic(origin, links[1], ex; p2 = vert11, spring = spring0, damper = damper0)
tra = tra[1]
ftr = ftr[1]
rot = rot[1]

constraintmat(tra)
constraintmat(ftr)
constraintmat(rot)
nullspacemat(tra)
nullspacemat(ftr)
nullspacemat(rot)




jointb1.inds
eqcs[1].inds
eqcs[2].inds