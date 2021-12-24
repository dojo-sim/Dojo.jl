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


nx = 5
Q = Diagonal(rand(nx))
q = - rand(nx)
c = 1.0
f(x) = 0.5 * x'*Q*x + q'*x + c
g(x) = Q*x + q
h(x) = Q

x0 = rand(nx)

xsol = bfgs_solve(f, g, x0)
