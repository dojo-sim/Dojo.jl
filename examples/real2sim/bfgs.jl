# Taken from Nocedal and Wright, Algorithm 6.1
using LinearAlgebra

function bfgs_solve(f, g, x0; H0=I(length(x0)), ϵtol=1e-4, iter=100)
    x = copy(x0)
    H = H0
    for k = 1:iter
        (norm(g(x), Inf) < ϵtol) && break
        p = - H * g(x)
        α = bfgs_linesearch(f, g, x, p)
        gprev = g(x)
        x = x + α*p
        s = α*p
        y = g(x) - gprev
        ρ = 1 / (y'*s)
        H = (I - ρ*s*y')*H*(I - ρ*s*y') + ρ*s*s'
        println("k:", k, "   x:", scn.(x), "H:", scn.(diag(H)))
    end
    return x
end

function bfgs_linesearch(f, g, x, p; iter=10)
    α = 1.0
    fprev = f(x)
    for k = 1:iter
        (f(x + α*p) <= fprev) && break
        α /= 2
    end
    return α
end

function strong_wolfe_linesearch(f, g, x, p;
    α0 = 0, αmax = 2.0, α1 = 1.0, iter = 20,
    c1 = 1e-4, c2 = 0.9,
    )
    ϕ(α) = f(x + α*p)
    pn = p ./ norm(p)
    ψ(α) = g(x + α*p)' * pn
    er = ϕ(0)
    gr = ψ(0)

    for i = 1:iter
        e0 = ϕ(α0)
        e1 = ϕ(α1)
        if (e1 > er + c1 * α1 * gr) || (e1 > e0 && i > 1)
            αstar = zoom(α0, α1)
            return αstar
        end
        g1 = ψ(α1)
        if abs(g1) <= -c2 * gr
            αstar = α1
            return αstar
        end
        if g1 >= 0
            αstar = zoom(α1, α0)
            return αstar
        end
        α0 = α1
        α1 = min(1.2*α1, αmax)
    end
    return nothing
end

function zoom(α)

nx = 5
Q = Diagonal(rand(nx))
q = - rand(nx)
c = 1.0
f(x) = 0.5 * x'*Q*x + q'*x + c
g(x) = Q*x + q
h(x) = Q

x0 = rand(nx)

xsol = bfgs_solve(f, g, x0)
