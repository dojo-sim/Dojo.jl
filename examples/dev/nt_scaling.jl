using LinearAlgebra
using Random

################################################################################
# No scaling
################################################################################
function pack_qp_vars(x, λ, γ, s)
    return [x; λ; γ; s]
end

function unpack_qp_vars(vars)
    off = 0
    x = vars[off .+ (1:nx)]; off += nx
    λ = vars[off .+ (1:nc)]; off += nc
    γ = vars[off .+ (1:ns)]; off += ns
    s = vars[off .+ (1:ns)]; off += ns
    return x, λ, γ, s
end

function qp_residual(vars, data, μ)
    P,C,G,h,A,b = data
    x, λ, γ, s = unpack_qp_vars(vars)
    r = [P*x + c + A'*λ + G'*γ;
         A*x - b;
         G*x + s - h + μ*ones(ns) ./γ - s;]
    return r
end

function qp_update(vars, α, Δ, μ)
    x, λ, γ, s = unpack_qp_vars(vars)
    v = deepcopy(vars)
    v[1:nx+nc+ns] += α * Δ
    Δγ = Δ[nx+nc+1:end]
    v[end-ns+1:end] += -s + μ*ones(ns) ./ γ - s .* α .* Δγ ./ γ
    return v
end

function qp_line_search(r, vars, Δ, data, μ)
    α = 1.0
    for i = 1:20
        vars_candidate = qp_update(vars, α, Δ, μ)
        _, _, γcand, scand = unpack_qp_vars(vars_candidate)
        all([γcand; scand] .> 0) && break
        α *= 0.5
    end
    for i = 1:20
        vars_candidate = qp_update(vars, α, Δ, μ)
        r_candidate = qp_residual(vars_candidate, data, μ)
        (norm(r_candidate, Inf) <= norm(r, Inf)) && break
        α *= 0.5
    end
    return α
end

function qp_solve!(vars, data)
    P,C,G,h,A,b = data
    μ = 1e-1
    iter = 0

    for i = 1:10
        for j = 1:20
            iter += 1
            x, λ, γ, s = unpack_qp_vars(vars)
            r = qp_residual(vars, data, μ)
            (norm(r, Inf) < 1e-6) && break
            R = [P A' G';
                 A zeros(nc,nc) zeros(nc,ns);
                 G zeros(ns,nc) -Diagonal(s ./ γ);]
            Δ = - R \ r
            α = qp_line_search(r, vars, Δ, data, μ)
            vars = qp_update(vars, α, Δ, μ)
            r0 = qp_residual(vars, data, μ)
            println("iter ", iter, " i ", i, "  j ", j,
                "  r", scn(norm(r0,Inf), exp_digits=2),
                "  α", scn(α),
                "  Δ", scn(norm(Δ,Inf), exp_digits=2),)
        end
        μ *= 0.1
    end
    return vars
end

################################################################################
# Nesterov Todd scaling
################################################################################
function qps_residual(vars, data, μ)
    P,C,G,h,A,b = data
    x, λ, γ, s = unpack_qp_vars(vars)
    W = Diagonal(sqrt.(s ./ γ))
    Wi = Diagonal(sqrt.(γ ./ s))
    γb = W*γ
    r = [P*x + c + A'*λ + G'*Wi*γb;
         A*x - b;
         G*x - h + μ*W*ones(ns) ./ γb;]
    return r
end

function qps_update(vars, α, Δ, μ)
    v = deepcopy(vars)
    x, λ, γ, s = unpack_qp_vars(vars)
    W = Diagonal(sqrt.(s ./ γ))
    Wi = Diagonal(sqrt.(γ ./ s))

    Δx = Δ[1:nx]
    Δλ = Δ[nx .+ (1:nc)]
    Δγb = Δ[nx + nc .+ (1:ns)]
    Δγ = Wi * Δγb

    v[1:nx] += α * Δx
    v[nx .+ (1:nc)] += α * Δλ
    v[nx + nc .+ (1:ns)] += α * Δγ
    v[end-ns+1:end] += -s + μ*ones(ns) ./ γ - s .* α .* Δγ ./ γ
    return v
end

function qps_line_search(r, vars, Δ, data, μ)
    α = 1.0
    for i = 1:20
        vars_candidate = qps_update(vars, α, Δ, μ)
        _, _, γcand, scand = unpack_qp_vars(vars_candidate)
        all([γcand; scand] .> 0) && break
        α *= 0.5
    end
    for i = 1:20
        vars_candidate = qps_update(vars, α, Δ, μ)
        r_candidate = qps_residual(vars_candidate, data, μ)
        (norm(r_candidate, Inf) <= norm(r, Inf)) && break
        α *= 0.5
    end
    return α
end

function qps_solve!(vars, data)
    P,C,G,h,A,b = data
    μ = 1e-1
    iter = 0

    for i = 1:10
        for j = 1:20
            iter += 1
            x, λ, γ, s = unpack_qp_vars(vars)
            W = Diagonal(sqrt.(s ./ γ))
            Wi = Diagonal(sqrt.(γ ./ s))
            r = qps_residual(vars, data, μ)
            (norm(r, Inf) < 1e-6) && break
            R = [P A' G'*Wi;
                 A zeros(nc,nc) zeros(nc,ns);
                 G zeros(ns,nc) -W;]
            Δ = - R \ r
            α = qps_line_search(r, vars, Δ, data, μ)
            vars = qps_update(vars, α, Δ, μ)
            r0 = qps_residual(vars, data, μ)
            println("iter ", iter, " i ", i, "  j ", j,
                "  r", scn(norm(r0,Inf), exp_digits=2),
                "  α", scn(α),
                "  Δ", scn(norm(Δ,Inf), exp_digits=2),)
        end
        μ *= 0.1
    end
    return vars
end


nx = 4
nc = 2
ns = 2

P = Diagonal(exp.(10randn(nx)))
c = exp.(5randn(nx))
G = exp.(5randn(ns,nx))
h = exp.(5randn(ns))
A = exp.(5randn(nc,nx))
b = exp.(5randn(nc))

data = (P,c,G,h,A,b)

x0 = rand(nx)
λ0 = rand(nc)
γ0 = rand(ns)
s0 = rand(ns)

vars0 = pack_qp_vars(x0, λ0, γ0, s0)

vars_sol0 = qp_solve!(vars0, data)
xsol0, λsol0, γsol0, ssol0 = unpack_qp_vars(vars_sol0)
γsol0 .* ssol0

vars_sol1 = qps_solve!(vars0, data)
xsol1, λsol1, γsol1, ssol1 = unpack_qp_vars(vars_sol1)
γsol1 .* ssol1
