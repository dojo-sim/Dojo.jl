# task objective push an object pack to the origin of the frame

# we want to find a policy that accomplishes the task successfully for various
# initializations and that is robust to disturbances

# we rely on an implicit policy
# u_star = argmin_u F(x,u;θ)
#           u ∈ U

# we make several choices
    # the task objective: simple cost function
    # the policy's optimization problem: simple QP with data parameterized by x and θ

# Training
    # we optimize the task objective over several rollouts wrt θ
    # dOBJ/dθ = dOBJ/dX * dX/dU * dU/dθ


################################################################################
# Smoothly Differentiable QP
################################################################################
using Symbolics
using LinearAlgebra
using BenchmarkTools


function pack_vars(u, s, γ)
    return [u; s; γ]
end

function pack_data(Qd, q, Av, b, c)
    return [Qd; q; Av; b; c]
end

function unpack_vars(vars::AbstractVector, nx::Int, nu::Int, nc::Int)
    ns = nc + 2nu
    off = 0
    u = vars[off .+ SUnitRange(1,nu)]; off += nu
    s = vars[off .+ SUnitRange(1,ns)]; off += ns
    γ = vars[off .+ SUnitRange(1,ns)]; off += ns
    return u, s, γ
end

function unpack_data(data::AbstractVector, nx::Int, nu::Int, nc::Int)
    off = 0
    Qd = data[off .+ SUnitRange(1,nu)]; off += nu
    q  = data[off .+ SUnitRange(1,nu)]; off += nu
    Av = data[off .+ SUnitRange(1,nc*nu)]; off += nc*nu
    b  = data[off .+ SUnitRange(1,nc)]; off += nc
    c  = data[off .+ SUnitRange(1,nu)]; off += nu
    return Qd, q, Av, b, c
end

function residual_slow(u, s, γ, ρ, Qd, q, Av, b, c)
    Q = Diagonal(Qd)
    A = reshape(Av, (nc,nu))
    D = [A; I(nu); -I(nu)]
    d = [b; c; c]
    r = [Q*u + q + D'*γ;
         D*u + d - s;
         γ .* s .- ρ]
    return r
end

function generate_residual(nx::Int, nu::Int, nc::Int)
    ns = nc+2nu
    Symbolics.@variables u[1:nu], s[1:ns], γ[1:ns], ρ[1:1]
    Symbolics.@variables Qd[1:nu], q[1:nu], Av[1:nc*nu], b[1:nc], c[1:nu]
    Q = Diagonal(Qd)
    A = reshape(Av, (nc,nu))
    D = [A; I(nu); -I(nu)]
    d = [b; c; c]

    vars = [u; s; γ]
    data = [Qd; q; Av; b; c]

    r = residual_slow(u, s, γ, ρ, Qd, q, Av, b, c)
    ∇vars = Symbolics.jacobian(r, vars)
    ∇data = Symbolics.jacobian(r, data)

    residual_code = build_function(r, vars, ρ, data)[2]
    ∇vars_residual_code = build_function(∇vars, vars, ρ, data)[2]
    ∇data_residual_code = build_function(∇data, vars, ρ, data)[2]
    return residual_code, ∇vars_residual_code, ∇data_residual_code
end

abstract type QPOptions
end

@with_kw struct QPOptions19{T} <: QPOptions
    outer_iterations::Int=10
    inner_iterations::Int=15
    residual_tol::T=1e-5
    iterations_linesearch::Int=4
    decay_relaxation::T=0.2
    decay_linesearch::T=0.5
    relaxation_tol::T=3e-4
    verbose::Bool=false
end

function solve(vars::AbstractVector{T}, data::AbstractVector{T};
        options::QPOptions=QPOptions19(), differentiate::Bool=false) where T
    ∇vars = zeros(nr, nr)
    Δvars = zeros(nr)
    r = zeros(nr)
    ρ = [0.1]
    cnt = 0
    for i = 1:options.outer_iterations
        (norm(r, Inf) < options.residual_tol) && (ρ[1] < options.relaxation_tol) && break
        for j = 1:options.inner_iterations
            cnt += 1
            residual_fast(r, vars, ρ, data)
            (norm(r, Inf) < options.residual_tol) && break
            ∇vars_residual_fast(∇vars, vars, ρ, data)
            Δvars = - ∇vars \ r
            α = linesearch(r, vars, Δvars, data, ρ, options)
            vars += α*Δvars
            options.verbose && println("I ", cnt,
                "  ρ", scn(ρ[1], digits=0),
                "  r", scn(norm(r, Inf), digits=0),
                "  ∇r", scn(norm(∇vars, 2), digits=0),
                )
        end
        ρ .*= options.decay_relaxation
    end
    options.verbose && println("I ", cnt,
        "  ρ", scn(ρ[1], digits=0),
        "  r", scn(norm(r, Inf), digits=0),
        "  ∇r", scn(norm(∇vars, 2), digits=0),
        )

    if differentiate
        ∇data = zeros(nr, nd)
        ∇data_residual_fast(∇data, vars, ρ, data)
        grad = -∇vars \ ∇data
        return vars, grad
    end
    return vars
end

function linesearch(r, vars, Δvars, data, ρ, options::QPOptions)
    α = 1.0
    r_candidate = zeros(nr)
    for i = 1:options.iterations_linesearch
        residual_fast(r_candidate, vars+α*Δvars, ρ, data)
        (norm(r_candidate, Inf) < norm(r, Inf)) && break
        α *= options.decay_linesearch
    end
    return α
end

nx = 1
nu = 1
nc = 2
ns = nc + 2nu
nr = nu + 2ns
nd = nu + nu + nu*nc + nc + nu
nθ = nd * (nx+1)

Qd = rand(nu)
q = rand(nu)
Av = rand(nc*nu)
b = rand(nc)
c = ones(nu)

A = reshape(Av, (nc,nu))
Q = Diagonal(Qd)
D = [A; I(nu); -I(nu)]
d = [b; c; c]


residual_code, ∇vars_residual_code, ∇data_residual_code = generate_residual(nx, nu, nc)
residual_fast = eval(residual_code)
∇vars_residual_fast = eval(∇vars_residual_code)
∇data_residual_fast = eval(∇data_residual_code)
r = zeros(nr)
∇vars = zeros(nr,nr)
∇data = zeros(nr,nd)
u = srand(nu)
s = srand(ns)
γ = srand(ns)
ρ = srand(1)
data = pack_data(Qd, q, Av, b, c)
vars = pack_vars(u, s, γ)
residual_fast(r, vars, ρ, data)
∇vars_residual_fast(∇vars, vars, ρ, data)
∇data_residual_fast(∇data, vars, ρ, data)
# @benchmark residual_fast($r, $vars, $ρ, $data)
# @benchmark ∇vars_residual_fast($∇vars, $vars, $ρ, $data)
# @benchmark ∇data_residual_fast($∇data, $vars, $ρ, $data)



u0 = 0.1*ones(nu)
s0 = 0.1*ones(ns)
γ0 = 0.1*ones(ns)
vars0 = [u0; s0; γ0]
vars0, grad0 = solve(vars0, data, differentiate=true, options=QPOptions19(verbose=true))
u1, s1, γ1 = unpack_vars(vars0, nx, nu, nc)

norm(D*u1 + d - s1, Inf)
norm(s1 .* γ1, Inf)

################################################################################
# Dataset
################################################################################
N = 20
X0 = [[i/(N-1)] for i=0:N-1]
U0 = [Float64.(x .> 0.5) for x in X0]
training_data0 = [X0, U0]

abstract type TrainingOptions
end

@with_kw struct TrainingOptions19{T} <: TrainingOptions
    iterations::Int=600
    cost_tol::T=1e-6
    gradient_tol::T=1e-5
    iterations_linesearch::Int=4
    decay_relaxation::T=0.1
    decay_linesearch::T=0.5
    relaxation_tol::T=3e-4
    hessian_reg::T=1e-5
    verbose::Bool=true
end

function train(θ, training_data; QP_options::QPOptions=QPOptions19(),
        training_options::TrainingOptions=TrainingOptions19())
    reg = training_options.hessian_reg
    for i = 1:training_options.iterations
        c, g, H = training_cost(θ, training_data, differentiate=true, QP_options=QP_options)
        (norm(c, Inf) < training_options.cost_tol) && (norm(g, Inf) < training_options.gradient_tol) && break
        Δθ = - (H + reg*I(nθ)) \ g
        α = linesearch(c, θ, Δθ, training_data, training_options)
        θ += α * Δθ
        training_options.verbose && println("I ", i,
            "  c", scn(c, digits=2),
            "  g", scn(norm(g, Inf), digits=0),
            "  H", scn(norm(H, Inf), digits=0),
            "  α", scn(norm(α, Inf), digits=0),
            "  Δθ", scn(norm(Δθ, Inf), digits=0),
            )
    end
    return θ
end

function linesearch(c, θ, Δθ, training_data, options::TrainingOptions)
    α = 1.0
    c_candidate = 0.0
    for i = 1:options.iterations_linesearch
        c_candidate = training_cost(θ+α*Δθ, training_data, differentiate=false)
        (norm(c_candidate, Inf) < norm(c, Inf)) && break
        α *= options.decay_linesearch
    end
    return α
end

function training_cost(θ, training_data; differentiate::Bool=false,
        QP_options::QPOptions=QPOptions19())
    X, U = training_data
    N = length(X)
    c = 0.0
    g = zeros(nθ)
    H = zeros(nθ, nθ)
    for i = 1:N
        if differentiate
            QPdata = reshape(θ, nd, nx+1) * [X[i]; 1.0] # add bias term to x
            vars, ∇data_full = solve(0.1*ones(nr), QPdata, differentiate=true, options=QP_options)
            u, s, γ = unpack_vars(vars, nx, nu, nc)
            ∇data = ∇data_full[1:nu, :] # extract gradient of u wrt to data
            # ∇θ = ∇data * ForwardDiff.jacobian(θ -> sigmoid.(reshape(θ, nd, nx+1)*[X[i]; 1.0]), θ)
            ∇θ = ∇data * ForwardDiff.jacobian(θ -> reshape(θ, nd, nx+1)*[X[i]; 1.0], θ)
            g += ∇θ' * (u - U[i])
            H += ∇θ' * ∇θ
        else
            QPdata = reshape(θ, nd, nx+1) * [X[i]; 1.0] # add bias term to x
            vars = solve(0.1*ones(nr), QPdata, differentiate=false, options=QP_options)
            u, s, γ = unpack_vars(vars, nx, nu, nc)
            # verbose && (@show scn.(u))
        end
        c += 0.5 * transpose(u - U[i]) * (u - U[i])
    end
    differentiate && return c/N, g/N, H/N
    return c/N
end

function sigmoid(x)
    return 1 / (1 + exp(-x))
end

function policy(x, θ; QP_options::QPOptions=QPOptions19())
    QPdata = reshape(θ, nd, nx+1) * [x; 1.0] # add bias term to x
    vars = solve(0.1*ones(nr), QPdata, differentiate=false, options=QP_options)
    u, s, γ = unpack_vars(vars, nx, nu, nc)
    return u
end


# TODO optimize this initial guess
Qd0 = 0.1ones(nu)
q0 = zeros(nu)
Av0 = 0.1*ones(nc*nu)
b0 = 0.1ones(nc)
c0 = 0.1ones(nu)
data0 = pack_data(Qd0, q0, Av0, b0, c0)
θ0 = reshape([zeros(nd, nx) data0], nθ)

training_cost(θ0, training_data0, differentiate=false)
training_cost(θ0, training_data0, differentiate=true)

QP_options0 = QPOptions19(relaxation_tol=3e-4)
training_options0 = TrainingOptions19(hessian_reg=1e-4)
θ1 = train(θ0, training_data0, QP_options=QP_options0, training_options=training_options0)

Xp = Vector(-1:0.005:2)
Up = [policy([x], θ1, QP_options=QP_options0)[1] for x in Xp]
scatter(Xp, Up)

# verify solver
# plot results for slver

# train using positive and negative examples
