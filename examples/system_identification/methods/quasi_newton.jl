# Taken from Nocedal and Wright, Algorithm 6.1
using LinearAlgebra

function quasi_newton_solve(f, fgH, x0; ftol=-Inf, gtol=1e-4, iter=100, α0=1.0,
        lower=-Inf, upper=Inf, reg=1e-3, reg_min=1e-9, reg_max=1e6,
        Δrot::Int=0, n_sample0 = 50, Δn_sample=2, n_sample_max=500)

    x = copy(x0)
    X = [copy(x0)]
    rot = 0.0
    n_sample = n_sample0
    ls_failure = false

    for k = 1:iter
        rot += Δrot
        n_sample = min(n_sample + Δn_sample, n_sample_max)
        fe, ge, He = fgH(x, rot=rot, n_sample=n_sample)
        He += reg * I
        if ls_failure
            reg = clamp(reg * 2, reg_min, reg_max)
        else
            reg = clamp(reg / 1.5, reg_min, reg_max)
        end
        @show reg
        ((norm(ge, Inf) < gtol) || (fe < ftol)) && break
        p = - He \ ge
        α, ls_failure = clamped_linesearch(f, x, p, fe, α0=α0, rot=rot, n_sample=n_sample, lower=lower, upper=upper)
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

function clamped_linesearch(f, x, p, fprev; α0=1.0, rot=0.0, n_sample=50, iter=4,
        lower=-Inf, upper=Inf)
    α = α0
    ls_failure = false
    for k = 1:iter
        xc = clamp.(x + α*p, lower, upper)
        (f(xc,rot=rot,n_sample=n_sample) <= fprev) && break
        α /= 3
        (k == iter) && (α = 0.001 / norm(p,Inf); @show α; ls_failure = true)
    end
    return α, ls_failure
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
