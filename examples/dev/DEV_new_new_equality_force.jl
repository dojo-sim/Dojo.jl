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
storage = simulate!(mech, 0.1, record = true, solver = :mehrotra!)

# visualize(mech, storage, vis = vis)

mech.eqconstraints[1].constraints[1]
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
function fdjac(f, x; δ = 1e-5)
    x = Vector(x)
    n = length(f(x))
    m = length(x)
    jac = zeros(n, m)
    for i = 1:m
        xp = deepcopy(x)
        xm = deepcopy(x)
        xp[i] += δ
        xm[i] -= δ
        jac[:,i] = (f(xp) - f(xm)) / (2δ)
    end
    return jac
end

function finitediff_helper(joint::AbstractJoint, pbody::AbstractBody, cbody::AbstractBody, Δt::T,
        evalf, jacf; diff_body::Symbol = :child) where {T}

    jac0 = jacf(joint, pbody, cbody, Δt)
    function f(x)
        v2 = x[1:3]
        ω2 = x[4:6]
        if diff_body == :parent
            cstate = deepcopy(cbody.state)
            pstate = deepcopy(pbody.state)
            pstate.vsol[2] = v2
            pstate.ωsol[2] = ω2
        elseif diff_body == :child
            cstate = deepcopy(cbody.state)
            cstate.vsol[2] = v2
            cstate.ωsol[2] = ω2
            if typeof(pbody) <: Origin
                return evalf(joint, cstate, Δt)
            else
                pstate = deepcopy(pbody.state)
            end
        end
        return evalf(joint, pstate, cstate, Δt)
    end

    if diff_body == :child
        v2 = cbody.state.vsol[2]
        ω2 = cbody.state.ωsol[2]
        x = [v2; ω2]
    elseif diff_body == :parent
        v2 = pbody.state.vsol[2]
        ω2 = pbody.state.ωsol[2]
        x = [v2; ω2]
    else
        error("invalid diff_body")
    end
    jac1 = fdjac(f, x)
    return jac0, jac1
end

Δt = 0.01
tra1 = mech.eqconstraints[1].constraints[1]
tra2 = mech.eqconstraints[2].constraints[1]
origin = mech.origin
body1 = collect(mech.bodies)[1]
body2 = collect(mech.bodies)[2]


jac0, jac1 = finitediff_helper(tra2, body1, body2, Δt, springforcea, ∂springforcea∂vela, diff_body = :parent)
norm(jac0 - jac1, Inf)
jac0, jac1 = finitediff_helper(tra2, body1, body2, Δt, damperforcea, ∂damperforcea∂vela, diff_body = :parent)
norm(jac0 - jac1, Inf)
jac0, jac1 = finitediff_helper(tra2, body1, body2, Δt, springforcea, ∂springforcea∂velb, diff_body = :child)
norm(jac0 - jac1, Inf)
jac0, jac1 = finitediff_helper(tra2, body1, body2, Δt, damperforcea, ∂damperforcea∂velb, diff_body = :child)
norm(jac0 - jac1, Inf)
jac0, jac1 = finitediff_helper(tra2, body1, body2, Δt, springforceb, ∂springforceb∂velb, diff_body = :child)
norm(jac0 - jac1, Inf)
jac0, jac1 = finitediff_helper(tra2, body1, body2, Δt, damperforceb, ∂damperforceb∂velb, diff_body = :child)
norm(jac0 - jac1, Inf)
jac0, jac1 = finitediff_helper(tra2, body1, body2, Δt, springforceb, ∂springforceb∂vela, diff_body = :parent)
norm(jac0 - jac1, Inf)
jac0, jac1 = finitediff_helper(tra2, body1, body2, Δt, damperforceb, ∂damperforceb∂vela, diff_body = :parent)
norm(jac0 - jac1, Inf)
jac0, jac1 = finitediff_helper(tra1, origin, body1, Δt, springforceb, ∂springforceb∂velb, diff_body = :child)
norm(jac0 - jac1, Inf)
jac0, jac1 = finitediff_helper(tra1, origin, body1, Δt, damperforceb, ∂damperforceb∂velb, diff_body = :child)
norm(jac0 - jac1, Inf)

damperforceb(tra1, body1.state, Δt)
damperforceb(tra1, origin, body1, Δt, body1.id)