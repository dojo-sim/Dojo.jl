using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
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
using JLD2

# ## scripts
include(joinpath(module_dir(), "examples/policy/methods/continuation.jl"))
include(joinpath(module_dir(), "examples/policy/methods/tvlqr.jl"))
include(joinpath(module_dir(), "DojoEnvironments/src",
    "quadruped/methods/template.jl"))


################################################################################
# ## system
################################################################################
gravity = -9.81
timestep = 0.02
friction_coefficient = 0.8
damper = 0.5
spring = 1.0
env = get_environment(:quadruped,
    representation=:minimal,
    timestep=timestep,
    contact_body=false,
    gravity=gravity,
    friction_coefficient=friction_coefficient,
    damper=damper,
    spring=spring,
    infeasible_control=true,
    vis=vis)

# ## dimensions
nx = env.num_states
nu = env.num_inputs
nu_infeasible = 6

################################################################################
# Load trajectory
################################################################################
file = JLD2.jldopen(joinpath(@__DIR__, "../data/trotting_forward.jld2"))
x_ref = file["x"]
u_ref = file["u"]
JLD2.close(file)

T = length(x_ref)

################################################################################
# TVLQR
################################################################################

K_tv, P_tv = tvlqr(x_ref, u_ref, env; q_tracking=q_tracking, r_tracking=r_tracking)

################################################################################
# tuning
################################################################################
function tracking_cost(x1, K; Q=Diagonal(ones(nx)))
    c = 0.0
    dcdx = [zeros(1, nx) for i = 1:T-1]
    dxdu = [zeros(nx, nu) for i = 1:T-1]
    dxpdx = [zeros(nx, nx) for i = 1:T-1]
    dudK = [spzeros(nu, nu*nx) for i = 1:T-1]
    dudx = [zeros(nu, nx) for i = 1:T-1]
    dcdK = [zeros(1, nu*nx) for i = 1:T-1]
    dxdK = [zeros(nx, nu*nx) for i = 1:T-1]
    dxdx = [zeros(nx, nx) for i = 1:T-1]
    dx = zeros(nx, nx)
    du = zeros(nx, nu)

    x = [zeros(nx) for i = 1:T]
    x[1] .= x1

    #Forward Pass
    for i = 1:T-1
        u_π = u_ref[i] + K[i] * (x_ref[i] - x[i])
        dynamics(x[i+1], env, x[i], u_π, zeros(0); gradients=true)
        dxpdx[i] .= env.dynamics_jacobian_state
        dxdu[i] .= env.dynamics_jacobian_input
        dudx[i] .= K[i]
        Δx = (x_ref[i] - x[i])
        for j = 1:nx
            dudK[i][:,(j-1)*nu .+ (1:nu)] = Δx[j] * I(nu)
        end
        dxdx[i] = dxdu[i] * dudx[i]
    end

    # Backward pass
    for i = T-1:-1:1
        c += 0.5 * (x[i+1] - x_ref[i+1])' * Q * (x[i+1] - x_ref[i+1])
        dcdx[i] .= (x[i+1] - x_ref[i+1])' * Q
        # for j = i:i
        for j = i:T-1
            X = I(nx)
            for k = j:-1:i+1
                X *= dxdx[k]
            end
            dcdK[i] .+= dcdx[j] * X *dxdu[j] * dudK[j]
        end
    end
    return x, c, dcdK
end

x, c, dcdK = tracking_cost(x_ref[1], K_tv, Q=Diagonal(ones(nx)))

# K_start = deepcopy(K_tv)
for i = 1:1
    c = 0.0
    dcdK = [zeros(1, nu*nx) for i=1:T-1]
    for j = 1:3
        x1 = deepcopy(x_ref[1])
        x1[j] += 0.05
        xi, ci, dcdKi = tracking_cost(x1, K_start, Q=Diagonal(q_tracking))
        c += ci / 3
        for k = 1:T-1
            dcdK[k] += dcdKi[k]/3
        end
    end
    # x, c, dcdK = tracking_cost(x_ref[1], K_start, Q=Diagonal(q_tracking))
    @show i, c
    for i = 1:T-1
        K_start[i] -= 0.001 * reshape(dcdK[i], (nu,nx))
    end
end

################################################################################
# Save
################################################################################
JLD2.jldsave(joinpath(@__DIR__, "../data/tvlqr_gains.jld2"), K=K_tv)
JLD2.jldsave(joinpath(@__DIR__, "../data/tuned_gains.jld2"), K=K_start)
file = JLD2.jldopen(joinpath(@__DIR__, "../data/tuned_gains.jld2"))
file["K"]
JLD2.close(file)




DojoEnvironments.visualize(env, x, build=false)
