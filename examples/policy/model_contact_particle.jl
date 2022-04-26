mutable struct ContactParticle19{T}
    nx::Int
    nu::Int
    nw::Int
    nv::Int
    nc::Int
    nvars::Int
    mass::T
    gravity::T
    side::T
    timestep::T
    friction_coefficient::T
    r_cache::Vector{T}
    r∂vars_cache::Matrix{T}
    r∂x_cache::Matrix{T}
    r∂u_cache::Matrix{T}
    fct_r::Function
    fct_r∂vars::Function
    fct_r∂x::Function
    fct_r∂u::Function
end

function ContactParticle19(;nx=6, nu=3, nw=0, nv=3, nc=1,
    mass = 1.0,
    gravity = 9.81,
    side = 0.5,
    timestep = 0.05,
    friction_coefficient = 0.2)
    nvars = 2*6nc + nv
    r_cache = zeros(nvars)
    r∂vars_cache = zeros(nvars, nvars)
    r∂x_cache = zeros(nvars, nx)
    r∂u_cache = zeros(nvars, nu)
    model = ContactParticle19(nx, nu, nw, nv, nc, nvars,
        mass, gravity, side, timestep, friction_coefficient,
        r_cache,
        r∂vars_cache,
        r∂x_cache,
        r∂u_cache,
        x -> x,
        x -> x,
        x -> x,
        x -> x,
        )
end

function unpack_vars(model::ContactParticle19, vars)
    v = vars[1:model.nv]
    Γ = vars[model.nv .+ (1:6*model.nc)]
    S = vars[model.nv + 6*model.nc .+ (1:6*model.nc)]
    return v, Γ, S
end

function pack_vars(model::ContactParticle19, v, Γ, S)
    return [v; Γ; S]
end

function residual(model::ContactParticle19, x, u, vars, κ)
    v, Γ, S = unpack_vars(model, vars)
    residual(model, x, u, v, Γ, S, κ)
end

function residual(model::ContactParticle19, x, u, v, Γ, S, κ)
    timestep = model.timestep
    mass = model.mass
    gravity = model.gravity
    cf = model.friction_coefficient
    p2 = x[1:3]
    v15 = x[4:6]
    v25 = v[1:3]
    p1 = p2 - timestep * v15
    p3 = p2 + timestep * v25

    vtan = v25[2:3]

    sγ = S[1:1]
    sψ = S[2:2]
    sβ = S[3:6]

    γ = Γ[1:1]
    ψ = Γ[2:2]
    β = Γ[3:6]

    N = [0 0 1]
    D = [1 0 0;
         0 1 0]
    P = [+D;
         -D]

    res = [
        mass * (p3 - 2p2 + p1)/timestep - timestep * mass * [0,0, gravity] - N' * γ - P' * β - u * timestep;
        sγ - (p3[3:3] .- model.side/2);
        sψ - (cf * γ - [sum(β)]);
        sβ - (P * v25 + ψ[1]*ones(4));
        Γ .* S .- κ[1];
        ]
    return res
end

function add_symbolics!(model::ContactParticle19)
    nx = model.nx
    nu = model.nu
    nvars = model.nvars
    @variables x_[1:nx]
    @variables u_[1:nu]
    @variables vars_[1:nvars]
    @variables vars_[1:nvars]
    @variables κ_[1:1]

    # generate expressions
    r = residual(model, x_, u_, vars_, κ_)
    r∂vars = Symbolics.jacobian(r, vars_)
    r∂x = Symbolics.jacobian(r, x_)
    r∂u = Symbolics.jacobian(r, u_)

    #build in place functions
    code_r = build_function(r, x_, u_, vars_, κ_)[2]
    code_r∂vars = build_function(r∂vars, x_, u_, vars_, κ_)[2]
    code_r∂x = build_function(r∂x, x_, u_, vars_, κ_)[2]
    code_r∂u = build_function(r∂u, x_, u_, vars_, κ_)[2]

    # eval code into functions
    fct_r = eval(code_r)
    fct_r∂vars = eval(code_r∂vars)
    fct_r∂x = eval(code_r∂x)
    fct_r∂u = eval(code_r∂u)
    model.fct_r = fct_r
    model.fct_r∂vars = fct_r∂vars
    model.fct_r∂x = fct_r∂x
    model.fct_r∂u  = fct_r∂u
    return nothing
end

function newton_solve(model::ContactParticle19, x, u; btol=1e-5, rtol=1e-2btol, verbose=false)
    κ0 = [btol]
    κ = ones(1)
    vars = ones(model.nvars)
    for i = 1:10
        rvio, bvio = violation!(model, x, u, vars, κ0)
        if (rvio <= rtol && bvio <= btol)
            verbose && println("rvio$(scn(rvio, exp_digits=2)) bvio$(scn(bvio, exp_digits=2)) ")
            break
        end
        verbose && println("i j  rvio     bvio     α        Δ")
        for j = 1:20

            model.fct_r(model.r_cache, x, u, vars, κ)
            model.fct_r∂vars(model.r∂vars_cache, x, u, vars, κ)
            rvio, bvio = violation(model)
            (rvio <= rtol && bvio <= btol) && break
            Δ = - model.r∂vars_cache \ model.r_cache
            α = conesearch(model, vars, Δ)
            α = linesearch(model, x, u, vars, κ, Δ, rvio, bvio, α)
            vars += α * Δ
            verbose && println(
                "$i " *
                "$j " *
                "$(scn(rvio, exp_digits=2)) " *
                "$(scn(bvio, exp_digits=2)) " *
                "$(scn(α, exp_digits=2)) " *
                "$(scn(norm(Δ, Inf), exp_digits=2)) "
            )
        end
        κ = max.(κ0, κ/10)
    end

    # Solution and gradients
    v = vars[1:model.nv]
    model.fct_r∂vars(model.r∂vars_cache, x, u, vars, κ0)
    model.fct_r∂x(model.r∂x_cache, x, u, vars, κ0)
    model.fct_r∂u(model.r∂u_cache, x, u, vars, κ0)
    v∂x = - model.r∂vars_cache \ model.r∂x_cache
    v∂u = - model.r∂vars_cache \ model.r∂u_cache
    return v, v∂x[1:3,:], v∂u[1:3,:]
end

function conesearch(model, vars::Vector{T}, Δ::Vector{T}) where T
    α = 1.0
    for k = 1:300
        all((0.95vars[model.nv+1:end] .+ α .* Δ[model.nv+1:end]) .> 0.0) && break
        α *= 0.5
        if k == 300
            @warn "conesearch"
        end
    end
    return α
end

function linesearch(model::ContactParticle19, x, u, vars, κ, Δ, rvio, bvio, α)
    for i = 1:10
        rvio_cand, bvio_cand = violation!(model, x, u, vars + α * Δ, κ)
        (rvio_cand <= rvio || bvio_cand < bvio) && break
        α *= 0.5
    end
    return α
end

function violation!(model::ContactParticle19, x, u, vars, κ)
    model.fct_r(model.r_cache, x, u, vars, κ)
    rvio, bvio = violation(model)
    return rvio, bvio
end

function violation(model::ContactParticle19)
    nv = model.nv
    rvio = norm(model.r_cache[1:nv], Inf)
    bvio = norm(model.r_cache[nv+1:end], Inf)
    return rvio, bvio
end

function dynamics(model::ContactParticle19, y, x, u; btol=1e-3, rtol=1e-2btol)
    v, v∂x, v∂u = newton_solve(model, x, u, btol=btol, rtol=rtol, verbose=false)
    y[1:3] .= x[1:3] + model.timestep * v
    y[4:6] .= v
    return nothing
end

function dynamics_jacobian_state(model::ContactParticle19, dx, x, u; btol=1e-3, rtol=1e-2btol)
    v, v∂x, v∂u = newton_solve(model, x, u, btol=btol, rtol=rtol, verbose=false)
    dx[1:3,:] .= model.timestep * v∂x
    dx[1:3,1:3] .+= I(3)
    dx[4:6,:] .= v∂x
    return v∂x
end

function dynamics_jacobian_input(model::ContactParticle19, du, x, u; btol=1e-3, rtol=1e-2btol)
    v, v∂x, v∂u = newton_solve(model, x, u, btol=btol, rtol=rtol, verbose=false)
    du[1:3,:] .= model.timestep * v∂u
    du[4:6,:] .= v∂u
    return nothing
end



# model = ContactParticle19()
# add_symbolics!(model)


#
# x0 = rand(model0.nx)
# u0 = rand(model0.nu)
# κ0 = [1e-5]
# vars0 = rand(model0.nvars)
# newton_solve(model0, x0, u0; btol=1e-5, rtol=1e-7)
# # @benchmark newton_solve(model0, x0, u0; btol=1e-5, rtol=1e-7)
#
# y0 = zeros(model0.nx)
# dx0 = zeros(model0.nx, model0.nx)
# du0 = zeros(model0.nx, model0.nu)
# @benchmark dynamics(model0, y0, x0, u0; btol=1e-5, rtol=1e-7)
# y0
# @benchmark dynamics_jacobian_state(model0, dx0, x0, u0; btol=1e-5, rtol=1e-7)
# dx0
# @benchmark dynamics_jacobian_input(model0, du0, x0, u0; btol=1e-5, rtol=1e-7)
# du0
#
