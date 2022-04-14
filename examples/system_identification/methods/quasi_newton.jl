# Taken from Nocedal and Wright, Algorithm 6.1
using LinearAlgebra

function quasi_newton_solve(f, fgH, x0; ftol=-Inf, gtol=1e-4, iter=100,
        lower=-Inf, upper=Inf, reg = 1e-3, Δrot::Int=0, n_sample0 = 50, Δn_sample=2, n_sample_max=500)
    x = copy(x0)
    X = [copy(x0)]
    rot = 0.0
    n_sample = n_sample0
    for k = 1:iter
        rot += Δrot
        n_sample = min(n_sample + Δn_sample, n_sample_max)
        fe, ge, He = fgH(x, rot=rot, n_sample=n_sample)
        He += reg * I
        ((norm(ge, Inf) < gtol) || (fe < ftol)) && break
        p = - He \ ge
        α = clamped_linesearch(f, x, p, fe, rot=rot, n_sample=n_sample, lower=lower, upper=upper)
        x = clamp.(x + α*p, lower, upper)
        push!(X, copy(x))
        println("k:", k,
            "   f:", scn(fe, digits=3),
            "   ∇:", scn.(norm(ge, Inf)),
            "   α:", scn.(α),
            # "   g:", scn.(-g(x, rot=rot)),
            # "   p:", scn.(p),
            "   x:", scn.(x),
            # "   H:", scn.(H),
            )
    end
    return x, X
end

function clamped_linesearch(f, x, p, fprev; rot=0.0, n_sample=50, iter=4,
        lower=-Inf, upper=Inf)
    α = 10.0
    for k = 1:iter
        xc = clamp.(x + α*p, lower, upper)
        (f(xc,rot=rot,n_sample=n_sample) <= fprev) && break
        α /= 3
        (k == iter) && (α = 0.01 / norm(p,Inf); @show α)
    end
    return α
end


# nx = 2
# Q = Diagonal(rand(nx))
# q = - rand(nx)
# c = 1.0
# f(x;rot=0) = 0.5 * x'*Q*x + q'*x + c
# g(x;rot=0) = Q*x + q
# H(x;rot=0) = Q
#
# x0 = rand(nx)
#
# xsol = quasi_newton_solve(f, g, H, x0)
