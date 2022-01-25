# Utils
function module_dir()
    return joinpath(@__DIR__, "..", "..", "..")
end

# Activate package
using Pkg
Pkg.activate(module_dir())

using MeshCat
# Open visualizer
vis = Visualizer()
open(vis)

# Include new files
include(joinpath(module_dir(), "examples", "loader.jl"))

using IterativeLQR

# System
gravity = -9.81
Δt = 0.05
mech = getmechanism(:atlas, Δt = Δt, g = gravity, cf = 1.5, damper = 10.0, spring = 0.0, model_type = :armless)
initialize!(mech, :atlas, tran = [0,0,0.], rot = [0,0,0.])

@elapsed storage = simulate!(mech, 0.05, record = true, solver = :mehrotra!, verbose = false)
visualize(mech, storage, vis = vis)

T = 24
xref = atlas_trajectory(mech, r = 0.10, z = 0.88; N = Int(T/2), Ncycles = 1)
zref = [min2max(mech, x) for x in xref]
storage = generate_storage(mech, zref)
visualize(mech, storage, vis = vis)
zref = [max2min(mech, z) for z in zref]

z1 = zref[1]
visualizeMaxCoord(mech, min2max(mech, z1), vis)

function gravity_compensation(mechanism::Mechanism)
    # only works with revolute joints for now
    nu = control_dimension(mechanism)
    u = zeros(nu)
    off  = 0
    for eqc in mechanism.eqconstraints
        nu = control_dimension(eqc)
        if eqc.parentid != nothing
            body = get_body(mechanism, eqc.parentid)
            rot = eqc.constraints[2]
            A = Matrix(nullspacemat(rot))
            Fτ = springforce(mechanism, eqc, body)
            F = Fτ[1:3]
            τ = Fτ[4:6]
            u[off .+ (1:nu)] = -A * τ
        else
            @warn "need to treat the joint to origin"
        end
        off += nu
    end
    return u
end

mech = getmechanism(:atlas, Δt = Δt, g = gravity, cf = 1.5, damper = 1000.0, spring = 30.0, model_type = :armless)
initialize!(mech, :atlas)
set_state!(mech, min2max(mech, zref[1]))
@elapsed storage = simulate!(mech, 0.05, record = true, solver = :mehrotra!, verbose = false)
visualize(mech, storage, vis = vis)
ugc = 6.00 * gravity_compensation(mech)

mech = getmechanism(:atlas, Δt = Δt, g = gravity, cf = 1.5, damper = 30.0, spring = 0.0, model_type = :armless)
control_dimension(mech)
u_control = ugc[6 .+ (1:15)]
u_mask = [zeros(15,6) I(15)]

z = [copy(z1)]
for t = 1:5
    znext = max2min(mech, step!(mech, min2max(mech, z[end]), u_mask'*u_control))
    push!(z, znext)
end
storage = generate_storage(mech, [min2max(mech, zi) for zi in z])
visualize(mech, storage, vis = vis)


# Model
function fd(y, x, u, w)
	z = step!(mech, min2max(mech, x), u_mask'*u, ϵ = 3e-4, btol = 3e-4, undercut = 1.5, verbose = false)
	y .= copy(max2min(mech, z))
end

function fdx(fx, x, u, w)
	fx .= copy(getMinGradients!(mech, min2max(mech, x), u_mask'*u, ϵ = 3e-4, btol = 3e-4, undercut = 1.5, verbose = false)[1])
end

function fdu(fu, x, u, w)
	∇u = copy(getMinGradients!(mech, min2max(mech, x), u_mask'*u, ϵ = 3e-4, btol = 3e-4, undercut = 1.5, verbose = false)[2])
	fu .= ∇u * u_mask'
end


# Time
h = mech.Δt
n, m, d = minimal_dimension(mech), 15, 0
dyn = Dynamics(fd, fdx, fdu, n, n, m, d)
model = [dyn for t = 1:T-1]

X[1][1:12]
X[1][13:18]

# Initial conditions, controls, disturbances
ū = [u_control for t = 1:T-1]
w = [zeros(d) for t = 1:T-1]

# Rollout


x̄ = rollout(model, z1, ū, w)
# step!(model.mech, x, u_mask'*u_control, ϵ = 1e-6, btol = 1e-6, undercut = 1.5, verbose = false)
# getGradients!(model.mech, x, u_mask'*u_control, ϵ = 1e-6, btol = 1e-3, undercut = 1.5, verbose = false)
storage = generate_storage(mech, [min2max(mech, x) for x in x̄])
visualize(mech, storage; vis = vis)

# Objective
# qt1 = [0.1; 0.1; 1.0; 0.001 * ones(3); 0.01 * ones(4); 0.01 * ones(3)]
# qt2 = [0.1; 0.1; 1.0; 0.001 * ones(3); 0.01 * ones(4); 0.01 * ones(3)]
# body_scale = [1; 0.1ones(12)]
# qt = vcat([body_scale[i] * [0.1 * ones(3); 0.001 * ones(3); 0.1 * ones(4); 0.01 * ones(3)] for i = 1:Nb]...)
qt = [0.3; 0.05; 0.05; 0.01 * ones(3); 0.01 * ones(3); 0.01 * ones(3); fill([0.2, 0.2], 15)...]

# ot1 = (x, u, w) -> transpose(x - zM) * Diagonal(Δt * qt) * (x - zM) + transpose(u) * Diagonal(Δt * 0.01 * ones(m)) * u
# ot2 = (x, u, w) -> transpose(x - zT) * Diagonal(Δt * qt) * (x - zT) + transpose(u) * Diagonal(Δt * 0.01 * ones(m)) * u
# oT = (x, u, w) -> transpose(x - zT) * Diagonal(Δt * qt) * (x - zT)
ots = [(x, u, w) -> transpose(x - zref[t]) * Diagonal(Δt * qt) * (x - zref[t]) + transpose(u) * Diagonal(Δt * 0.01 * ones(m)) * u for t = 1:T-1]
oT = (x, u, w) -> transpose(x - zref[end]) * Diagonal(Δt * qt) * (x - zref[end])

# ct1 = Cost(ot1, n, m, d)
# ct2 = Cost(ot2, n, m, d)
# cT = Cost(oT, n, 0, 0)
cts = Cost.(ots, n, m, d)
cT = Cost(oT, n, 0, 0)
# obj = [[ct1 for t = 1:10]..., [ct2 for t = 1:10]..., cT]
obj = [cts..., cT]

# Constraints
function goal(x, u, w)
	# Δ = x - zT
    Δ = x - zref[end]
    return Δ[collect(1:42)]
end

cont = Constraint()
conT = Constraint(goal, n, 0)
cons = [[cont for t = 1:T-1]..., conT]

prob = problem_data(model, obj, cons)
initialize_controls!(prob, ū)
initialize_states!(prob, x̄)

using BenchmarkTools
# Solve
@profiler IterativeLQR.constrained_ilqr_solve!(prob,
    verbose = true,
	linesearch=:armijo,
    α_min=1.0e-5,
	con_tol = 1e-2,
    obj_tol=1.0e-3,
    grad_tol=1.0e-3,
	# max_iter=100,
    max_iter=1,
	# max_al_iter=5,
    max_al_iter=1,
    ρ_init=1.0,
    ρ_scale=10.0)

using BenchmarkTools
using InteractiveUtils
prob.m_data.obj.costs[1].val(prob.m_data.obj.costs[1].val_cache, x̄[1], ū[1], w[1])
@benchmark prob.m_data.obj.costs[1].hessxx($prob.m_data.obj.costs[1].hessxx_cache, $x̄[1], $ū[1], $w[1])
@code_warntype prob.m_data.obj.costs[1].hessxx(prob.m_data.obj.costs[1].hessxx_cache, x̄[1], ū[1], w[1])

x_sol, u_sol = get_trajectory(prob)
storage = generate_storage(mech, [min2max(mech, x) for x in x_sol])
visualize(mech, storage, vis = vis)



function my_eval_obj_hess!(hessxx::Vector{Matrix{S}}, hessuu::Vector{Matrix{S}},
    hessux::Vector{Matrix{S}}, obj::IterativeLQR.Objective{S},
	x::Vector{Vector{S}}, u::Vector{Vector{S}}, w::Vector{Vector{S}}) where {S}
    T = length(obj)
    for (t, cost) in enumerate(obj)
        cost.hessxx(cost.hessxx_cache, x[t], u[t], w[t])
        @views hessxx[t] .+= cost.hessxx_cache
        fill!(cost.hessxx_cache, 0.0) # TODO: confirm this is necessary
        t == T && continue
        cost.hessuu(cost.hessuu_cache, x[t], u[t], w[t])
        cost.hessux(cost.hessux_cache, x[t], u[t], w[t])
        @views hessuu[t] .+= cost.hessuu_cache
        @views hessux[t] .+= cost.hessux_cache
        fill!(cost.hessuu_cache, 0.0) # TODO: confirm this is necessary
        fill!(cost.hessux_cache, 0.0) # TODO: confirm this is necessary
    end
end


obj = prob.m_data.obj.costs
hessxx = prob.m_data.obj_deriv.gxx
hessuu = prob.m_data.obj_deriv.guu
hessux = prob.m_data.obj_deriv.gux

my_eval_obj_hess!(hessxx, hessuu, hessux, obj, x̄, [ū..., zeros(0)], [w..., zeros(0)])
@benchmark my_eval_obj_hess!($hessxx, $hessuu, $hessux, $obj, $x̄, $([ū..., zeros(0)]), $([w..., zeros(0)]))
@code_warntype my_eval_obj_hess!(hessxx, hessuu, hessux, obj, x̄, [ū..., zeros(0)], [w..., zeros(0)])

x̄
ū
w

[ū..., zeros(0)]
