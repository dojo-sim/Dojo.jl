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
function TorqueRevolute(body1::AbstractBody{T}, body2, axis; p1 = szeros(T, 3), p2 = szeros(T, 3), qoffset = one(UnitQuaternion{T}), spring = zero(T), damper = zero(T)) where T
    return Translational3{T}(body1, body2, p1=p1, p2=p2),
        Torque1{T}(body1, body2; axis=axis, qoffset=qoffset, spring=spring, damper=damper),
        Rotational2{T}(body1, body2, axis=axis, qoffset=qoffset)
end

################################################################################
# DEVELOPMENT NEW EQUALITY CONSTRAINT
################################################################################
T = Float64

# Parameters
ex = [1.; 0; 0]
h = 1.
r = .05
vert11 = [0; 0; h/2]
vert12 = -vert11
Nlink = 10

# Links
origin = Origin{T}()
links = [Cylinder(r, h, h, color = RGBA(1., 0., 0.)) for i = 1:Nlink]

# Constraints
spring0 = 1.0 * 1e1
damper0 = 2.0 * 1e2
spring1 = 3.0 * 1e1
damper1 = 1.0 * 1e2
# jointb1 = EqualityConstraint(TorqueRevolute(origin, links[1], ex; spring=spring0, damper=damper0, p2 = vert11))
jointb1 = EqualityConstraint(Revolute(origin, links[1], ex; spring=spring0, damper=damper0, p2 = vert11))
if Nlink > 1
    eqcs = [
        jointb1;
        # [EqualityConstraint(TorqueRevolute(links[i - 1], links[i], ex; spring=spring1, damper=damper1, p1=vert12, p2=vert11)) for i = 2:Nlink]
        [EqualityConstraint(Revolute(links[i - 1], links[i], ex; spring=spring1, damper=damper1, p1=vert12, p2=vert11)) for i = 2:Nlink]
        ]
else
    eqcs = [jointb1]
end

mech = Mechanism(origin, links, eqcs, g = -9.81, Δt = 0.01)

initialize!(mech, :npendulum)
@elapsed storage = simulate!(mech, 0.3, record = true, solver = :mehrotra!, verbose = true)
@profiler storage = simulate!(mech, 5.0, record = true, solver = :mehrotra!, verbose = false)

visualize(mech, storage, vis = vis)






################################################################################
# Differentiation
################################################################################
include(joinpath(module_dir(), "examples", "dev", "diff_tools.jl"))
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
@show norm(full_vector(mech.system), Inf)

# finite diff
fd_datamat = finitediff_data_matrix(mech, data, sol) * attjac
@test norm(fd_datamat + datamat, Inf) < 1e-7
@test norm((fd_datamat + datamat)[:, 1:24], Inf) < 1e-7
# plot(Gray.(abs.(datamat)))
# plot(Gray.(abs.(fd_datamat)))
norm((fd_datamat + datamat)[1:6, 1:13], Inf)
norm((fd_datamat + datamat)[7:12, 1:13], Inf)
norm((fd_datamat + datamat)[7:12, 1:3], Inf)
norm((fd_datamat + datamat)[7:12, 4:6], Inf)
norm((fd_datamat + datamat)[7:12, 7:10], Inf)
norm((fd_datamat + datamat)[7:9, 7:10], Inf)
norm((fd_datamat + datamat)[10:12, 7:10], Inf)
norm((fd_datamat + datamat)[7:12, 11:13], Inf)

norm((fd_datamat + datamat)[1:12, 1:26], Inf)
norm((fd_datamat + datamat)[13:18, 1:13], Inf)
norm((fd_datamat + datamat)[13:18, 14:26], Inf)
datamat[13:18, 24:26]
-fd_datamat[13:18, 24:26]

norm((fd_datamat + datamat)[19:24, 1:26], Inf)
norm((fd_datamat + datamat)[19:24, 1:13], Inf)
norm((fd_datamat + datamat)[19:24, 14:26], Inf)
datamat[19:24, 24:26]
-fd_datamat[19:24, 24:26]

datamat[12:22, 24:26]
-fd_datamat[12:22, 24:26]

datamat

norm((fd_datamat + datamat)[1:24, 1:24], Inf)
norm((fd_datamat + datamat)[1:24, 25:26], Inf)
(fd_datamat + datamat)[1:12, 25:26]
(fd_datamat + datamat)[13:24, 25:26]
fd_datamat[13:24, 25:26]
datamat[13:24, 25:26]

fd_solmat = finitediff_sol_matrix(mech, data, sol)
@test norm(fd_solmat + solmat, Inf) < 1e-7
plot(Gray.(abs.(solmat)))
plot(Gray.(abs.(fd_solmat)))

fd_sensi = finitediff_sensitivity(mech, data) * attjac
@test norm(fd_sensi - sensi) / norm(fd_sensi) < 5e-3
# plot(Gray.(sensi))
# plot(Gray.(fd_sensi))
# norm(fd_sensi - sensi, Inf)
# norm(fd_sensi, Inf)
# norm(fd_sensi - sensi) / norm(fd_sensi)

################################################################################
# Finite Diff
################################################################################

function fdjac(f, x; δ = 1e-5)
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

begin
    eqc1 = collect(mech.eqconstraints)[1]
    eqc2 = collect(mech.eqconstraints)[2]
    tra1 = eqc1.constraints[1]
    rot1 = eqc1.constraints[3]
    torque1 = eqc1.constraints[2]
    tra2 = eqc2.constraints[1]
    rot2 = eqc2.constraints[3]
    torque2 = eqc2.constraints[2]
    A1 = constraintmat(torque1)
    A1ᵀ = zerodimstaticadjoint(A1)
    A2 = constraintmat(torque2)
    A2ᵀ = zerodimstaticadjoint(A2)
end
torque1
tra1
constraintmat(tra1)
rot1
zerodimstaticadjoint(constraintmat(rot1)) * srand(1)

qa = UnitQuaternion(rand(4)...)
ωa = rand(3)
qb = UnitQuaternion(rand(4)...)
ωb = rand(3)
Δt = 0.01

Xb, Qb = ∂g∂posb(torque1, qb, ωb, Δt)
Xb
Qb

Qb_fd = fdjac(
    qb -> A1ᵀ * A1 * (springtorque(torque1, UnitQuaternion(qb..., false)) + dampertorque(torque1, UnitQuaternion(qb..., false), ωb)),
    [qb.w, qb.x, qb.y, qb.z])
norm(Qb - Qb_fd, Inf) / (norm(Qb_fd) + 1) < 1e-7
norm(Qb - Qb_fd, Inf) / (norm(Qb_fd) + 1)


Xb, Qb = ∂g∂posb(torque2, qa, ωa, qb, ωb, Δt)
Xb
Qb

Qb_fd = fdjac(
    qb -> A2ᵀ * A2 * (springtorque(torque2, qa, UnitQuaternion(qb..., false)) + dampertorque(torque2, qa, ωa, UnitQuaternion(qb..., false), ωb)),
    [qb.w, qb.x, qb.y, qb.z])
norm(Qb - Qb_fd, Inf) / (norm(Qb_fd) + 1) < 1e-7
norm(Qb - Qb_fd, Inf) / (norm(Qb_fd) + 1)

function fm(qa)
    q = UnitQuaternion(qa..., false) * qb
    return [q.w, q.x, q.y, q.z]
end
ForwardDiff.jacobian(fm, [qa.w, qa.x, qa.y, qa.z]) - Rmat(qb)


Xa, Qa = ∂g∂posa(torque2, qa, ωa, qb, ωb, Δt)
Xa
Qa

Qa_fd = fdjac(
    qa -> A2ᵀ * A2 * (springtorque(torque2, UnitQuaternion(qa..., false), qb) + dampertorque(torque2, UnitQuaternion(qa..., false), ωa, qb, ωb)),
    [qa.w, qa.x, qa.y, qa.z])
norm(Qa - Qa_fd, Inf) / (norm(Qa_fd) + 1) < 1e-7
norm(Qa - Qa_fd, Inf) / (norm(Qa_fd) + 1)

Qa - Qa_fd
Qa_fd
torque2.spring
torque2.damper
torque1.qoffset
torque2.qoffset

A2ᵀ * A2

#
#
# x3a = rand(3)
# q3a = rand(4) #UnitQuaternion(rand(4)...)
# x3b = rand(3)
# q3b = rand(4) #UnitQuaternion(rand(4)...)
# qoff = rand(4) #UnitQuaternion(rand(4)...)
# λ = rand(3)
#
# Gtλ_rot(x3a, q3a, x3b, q3b, qoff, λ)
# rot1.qoffset = UnitQuaternion(qoff..., false)
# fd = ForwardDiff.jacobian(vars -> Gtλ_rot(vars[1:3], vars[4:7], x3b, q3b, qoff, λ), [x3a; q3a])
# sb = _dG(rot1, x3a, UnitQuaternion(q3a..., false), x3b, UnitQuaternion(q3b..., false), λ)
# norm(fd - sb)
#
#
# tra = 2 .* ones(3, 7)
# rot = zeros(2, 7)
# tor = ones(1, 7)
#
# vv = [tra, rot, tor]
# vcat(vv...)





#
# eqc2.parentid
# eqc2.constraints
#
#
# eqc2.childids
# #
# function _dGa(mechanism, pbody::Body, cbody::Body, eqc::EqualityConstraint{T,N,Nc}) where {T,N,Nc} # 6 x 7
#     Δt = mechanism.Δt
#     dG = zeros(6,7)
#
#     off = 0
#     for i = 1:Nc
#         joint = eqc.constraints[i]
#         Nj = length(joint)
#         dG += _dG(joint, pbody, cbody, eqc.λsol[2][off .+ (1:Nj)], Δt) # 6 * 7 = ∂(6 x d * d)/∂(7)
#         off += Nj
#     end
#     return dG
# end
#
# function _dG(joint::AbstactJoint, pbody::Body, cbody::Body, λ::AbtractVector, Δt)
#     A = constraintmat(joint)
#     Aᵀ = zerodimstaticadjoint(A)
#     x3a, q3a = posargsnext(pbody.state, Δt)
#     x3b, q3b = posargsnext(cbody.state, Δt)
#     _dG(joint, x3a, q3a, x3b, q3b, Aᵀ * λ)
#     return dG
# end
#
# @inline function _dG(mechanism, bodya::Body, bodyb::Body)
#     body.id == constraint.parentid ? (return _dGaa(mechanism, constraint, body)) : (return _dGbb(mechanism, constraint, body))
# end
#

#
# using Symbolics
#
# @variables x3a[1:3], q3a[1:4], x3b[1:3], q3b[1:4], λ[1:3], qoff[1:4]
#
# function Gatλ_rot(x3a, q3a, x3b, q3b, qoff, λ)
#     T = eltype(x3a)
#     q3a_ = UnitQuaternion(q3a..., false)
#     q3b_ = UnitQuaternion(q3b..., false)
#     qoff_ = UnitQuaternion(qoff..., false)
#
#     X = szeros(T, 3, 3)
#     Q = VRᵀmat(qoff_) * Rmat(q3b_) * Tmat() * LVᵀmat(q3a_)
#     return [X'; Q']* λ
# end
#
# function Gbtλ_rot(x3a, q3a, x3b, q3b, qoff, λ)
#     T = eltype(x3a)
#     q3a_ = UnitQuaternion(q3a..., false)
#     q3b_ = UnitQuaternion(q3b..., false)
#     qoff_ = UnitQuaternion(qoff..., false)
#
#     X = szeros(T, 3, 3)
#     Q = VRᵀmat(qoff_) * Lᵀmat(q3a_) * LVᵀmat(q3b_)
#
#     return [X'; Q']* λ
# end
#
# @show jac = Symbolics.jacobian(Gbtλ_rot(x3a, q3a, x3b, q3b, qoff, λ), [x3a; q3a])
# for i = 4:6
#     for j = 1:7
#         # @show i j
#         println("dG[$i, $j] = ", jac[i,j])
#     end
# end
#
# x3a = rand(3)
# q3a = rand(4) #UnitQuaternion(rand(4)...)
# x3b = rand(3)
# q3b = rand(4) #UnitQuaternion(rand(4)...)
# qoff = rand(4) #UnitQuaternion(rand(4)...)
# λ = rand(3)
#
# Gbtλ_rot(x3a, q3a, x3b, q3b, qoff, λ)
# rot1.qoffset = UnitQuaternion(qoff..., false)
# fd = ForwardDiff.jacobian(vars -> Gbtλ_rot(vars[1:3], vars[4:7], x3b, q3b, qoff, λ), [x3a; q3a])
# sb = _dGba(rot1, x3a, UnitQuaternion(q3a..., false), x3b, UnitQuaternion(q3b..., false), λ)
# norm(fd - sb)
#
# fd
#
# sb
# tra = 2 .* ones(3, 7)
# rot = zeros(2, 7)
# tor = ones(1, 7)
#
# vv = [tra, rot, tor]
# vcat(vv...)
# ∂gab∂ʳba














#
#
#
#
#
#
# using BenchmarkTools
# function ∂integration(q2::UnitQuaternion{T}, ω2::SVector{3,T}, Δt::T) where {T}
#     Δ = Δt * SMatrix{3,3,T,9}(Diagonal(sones(T,3)))
#     X = hcat(Δ, szeros(T,3,3))
#     Q = hcat(szeros(T,4,3), Lmat(q2)*derivωbar(ω2, Δt)*Δt/2)
#     return svcat(X, Q)
# end
#
# function NNN()
#     return ∂integration(q2_, ω2_, Δt_)
# end
#
# ω2_ = srand(3)
# q2_ = UnitQuaternion(rand(4)...)
# Δt_ = 0.1
#
# ∂integration(q2_, ω2_, Δt_)
# NNN()
# @benchmark ∂integration(q2_, ω2_, Δt_)
# @benchmark NNN()
#
# aa = sones(3,3)
# aa = SMatrix{}Diagonal(sones(3))
# bb = szeros(4,3)
# cc = szeros(3,4)
# svcat(aa, bb)
# hcat(aa, cc)
# svcat(aa, bb)
#
# a = 10
# a = 10
