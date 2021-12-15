using FiniteDiff 
using BenchmarkTools

struct MuJoCoModel{T,CX,CU}
    m::jlModel 
    d::jlData
    x::Vector{T}
    fx::Matrix{T} 
    fu::Matrix{T}
    nx::Int 
    nu::Int
    idx_u::Vector{Int} 
    idx_pos::Vector{Int} 
    idx_vel::Vector{Int} 
    idx_next::Vector{Int}
    cache_x::CX 
    cache_u::CU
end

function MuJoCoModel(path)
    m = jlModel(path)
    d = jlData(m)

    nq = length(d.qpos)
    nv = length(d.qvel)
    nx = nq + nv
    nu = length(d.ctrl)

    x = zeros(nx) 
    u = zeros(nu) 
    fx = zeros(nx, nx) 
    fu = zeros(nx, nu) 
   
    idx_u = collect(1:nu) 
    idx_pos = collect(1:nq) 
    idx_vel = collect(nq .+ (1:nv))
    idx_next = collect(1:nq) 

    cache_x = FiniteDiff.JacobianCache(x, x)
    cache_u = FiniteDiff.JacobianCache(u, x)

    MuJoCoModel(m, d, x, 
        fx, fu, 
        nx, nu, 
        idx_u, idx_pos, idx_vel, idx_next,
        cache_x, cache_u) 
end

# dynamics methods
function f!(y, model::MuJoCoModel, x, u) 
    model.d.ctrl .= @views u[model.idx_u]
    model.d.qpos[model.idx_next] .= @views x[model.idx_pos] 
    model.d.qvel[model.idx_next] .= @views x[model.idx_vel] 
    mj_step(model.m, model.d) 
    y[model.idx_pos] .= @views model.d.qpos[model.idx_next]
    y[model.idx_vel] .= @views model.d.qvel[model.idx_next]
    return y
end

f(model::MuJoCoModel, x, u) = f!(model.x, model, x, u)

function fx!(dx, model::MuJoCoModel, x, u) 
    FiniteDiff.finite_difference_jacobian!(dx, (a, b) -> f!(a, model, b, u), x, model.cache_x)
    return dx 
end

fx(model::MuJoCoModel, x, u) = fx!(model.fx, model, x, u)

function fu!(du, model::MuJoCoModel, x, u) 
    FiniteDiff.finite_difference_jacobian!(du, (a, b) -> f!(a, model, x, b), u, model.cache_u)
    return du
end

fu(model::MuJoCoModel, x, u) = fu!(model.fu, model, x, u)

# x2 = zeros(acrobot.nx)
# x1 = zeros(acrobot.nx)
# u1 = zeros(acrobot.nu)
# f(acrobot, x1, u1)
# @benchmark f($acrobot, $x1, $u1)
# f!(x2, acrobot, x1, u1)
# @benchmark f!($x2, $acrobot, $x1, $u1)
# fx(acrobot, x1, u1)
# @benchmark fx($acrobot, $x1, $u1)
# fu(acrobot, x1, u1)
# @benchmark fu($acrobot, $x1, $u1)
