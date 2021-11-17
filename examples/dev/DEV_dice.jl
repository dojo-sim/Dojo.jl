################################################################################
# Development
################################################################################
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
include(joinpath(module_dir(), "examples", "loader.jl"))

# mech = getmechanism(:dice, Δt = 0.01, g = -9.81, cf = 0.2, contact = true, mode=:box, conetype = :soc)
mech = getmechanism(:dice, Δt = 0.01, g = -9.81, cf = 0.2, contact = true, mode=:box, conetype = :linear)
# mech = getmechanism(:dice, Δt = 0.01, g = -9.81, contact = true, mode=:box, conetype = :impact)
Random.seed!(100)
ω = 5.0 * (rand(3) .- 0.5) * 1
x = [0, 0, 1.0]
v = 100.0 * [1, 0.3, 0.0]
initialize!(mech, :dice, x = x, v = v, ω = ω)
storage = simulate!(mech, 15.1, record = true, solver = :mehrotra!, verbose = false)
visualize(mech, storage, vis = vis)






ineqc1 = collect(mech.ineqconstraints)[1]
cont = ineqc1.constraints[1]
body = getbody(mech, ineqc1.parentid)
x, v, q, ω = fullargssol(body.state)
x3, q3 = posargsnext(body.state, mech.Δt)

# transforms the velocities of the origin of the link into velocities along all 4 axes of the friction pyramid
Bxmat = cont.Bx
Bqmat = Bxmat * ∂vrotate∂q(cont.p, q3) * LVᵀmat(q3)
γsol = rand(6)
ssol = rand(6)
function _g(x3, q3, γsol, ssol)
    γ = γsol[1]
    sγ = ssol[1]
    ψ = γsol[2]
    sψ = ssol[2]
    η = γsol[@SVector [3,4,5,6]]
    sη = ssol[@SVector [3,4,5,6]]
    return SVector{6,Float64}(
        cont.ainv3 * (x3 + vrotate(cont.p,q3) - cont.offset) - sγ,
        cont.cf * γ - sum(sη) - sψ,
        (Bxmat * v + Bqmat * ω + ψ * sones(4) - η)...)
end
function _dgdpos(x3, q3, γsol, ssol)
    _g_local(vars) = _g(vars[1:3], UnitQuaternion(vars[4:7]..., false), γsol, ssol)
    X = FiniteDiff.finite_difference_jacobian(_g_local, [x3; q3.w; q3.x; q3.y; q3.z])[:,1:3]
    Q = FiniteDiff.finite_difference_jacobian(_g_local, [x3; q3.w; q3.x; q3.y; q3.z])[:,4:7] * LVᵀmat(q3)
    return X, Q
end

g0 = _g(x3, q3, γsol, ssol)
X, Q = ∂g∂pos(cont, x3, q3)
X0, Q0 = _dgdpos(x3, q3, γsol, ssol)
X
X0
X - X0
norm(X - X0)
norm(Q - Q0)

# - XQ0








################################################################################
# Differentiation
################################################################################

include(joinpath(module_dir(), "examples", "diff_tools.jl"))
# Set data
data = getdata(mech)
setdata!(mech, data)
sol = getsolution(mech)
Nb = length(collect(mech.bodies))
attjac = attitudejacobian(data, Nb)

# IFT
setentries!(mech)
datamat = full_data_matrix(mech)
solmat = full_matrix(mech.system)
sensi = - (solmat \ datamat)
@show cond(solmat)
@show rank(solmat)
@show norm(full_vector(mech.system), Inf)

# finite diff
fd_datamat = finitediff_data_matrix(mech, data, sol) * attjac
@test norm(fd_datamat + datamat, Inf) < 1e-7
plot(Gray.(abs.(datamat)))
plot(Gray.(abs.(fd_datamat)))

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

###############################################################################
# plot
###############################################################################

plot(hcat(Vector.(storage.x[1])...)')
plot(hcat([[q.w, q.x, q.y, q.z] for q in storage.q[1]]...)')
plot(hcat(Vector.(storage.v[1])...)')
plot(hcat(Vector.(storage.ω[1])...)')

sdf = get_sdf(mech, storage)
plot(hcat(sdf[1]...)', ylims = (-0.01,0.01))
plot(hcat(sdf[2]...)', ylims = (-0.01,0.01))
plot(hcat(sdf[3]...)', ylims = (-0.01,0.01))
plot(hcat(sdf[4]...)', ylims = (-0.01,0.01))
plot(hcat(sdf[5]...)')
plot(hcat(sdf[6]...)')
plot(hcat(sdf[7]...)')
plot(hcat(sdf[8]...)')
