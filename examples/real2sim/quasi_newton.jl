# Taken from Nocedal and Wright, Algorithm 6.1
using LinearAlgebra

function bfgs_solve(f, g, x0; H0=Matrix(Diagonal(ones(length(x0)))),
    ftol = -Inf, gtol=1e-4, iter=100, lower=-Inf, upper=Inf)
    x = copy(x0)
    H = H0
    Ix = Matrix(Diagonal(ones(length(x0))))
    rot = 0.0
    for k = 1:iter
        # rot += 1
        fe = f(x, rot=rot)
        ge = g(x, rot=rot)
        ((norm(ge, Inf) < gtol) || (fe < ftol)) && break
        p = - H * ge
        α = bfgs_linesearch(f, g, x, p, rot=rot, lower=lower, upper=upper)
        # α = strong_wolfe_linesearch(f, g, x, p, rot=rot, lower=lower, upper=upper)
        gprev = ge
        x = clamp.(x + α*p, lower, upper)
        s = α*p
        sn = s ./ norm(s)
        y = g(x, rot=rot) - gprev
        y = y - min(y'*sn, -0.1)*sn
        ρ = 1 / (y'*s) #TODO fix this
        # if y'*s > 1e-10
            # @show "here"
        H = (Ix - ρ*s*y')*H*(Ix - ρ*s*y') + ρ*s*s'
        @show isposdef(H)
        # end
        # reg = 1e-6
        # while !isposdef(H)
        #     @show reg
        #     @show isposdef(H)
        #     reg *= 2
        #     (reg > 1) && break
        #     H = (1-reg) * H + reg*Ix
        # end

        println("k:", k, "    f:", scn(f(x,rot=rot), digits=3), "   ∇:", scn.(norm(g(x,rot=rot), Inf)),
            "    g:", scn.(-g(x, rot=rot)),
            "    p:", scn.(p),
            "    α:", scn.(α),
            "    x:", scn.(x),
            # "    H:", scn.(H),
            )
    end
    return x
end

function bfgs_linesearch(f, g, x, p; rot=0.0, iter=10,
        lower=-Inf, upper=Inf)
    α = 1.0
    fprev = f(x,rot=rot)
    for k = 1:iter
        xc = clamp.(x + α*p, lower, upper)
        (f(xc,rot=rot) <= fprev) && break
        α /= 2
    end
    return α
end

function strong_wolfe_linesearch(f, g, x, p; rot=0.0,
    α0=0.1, αmax=2.0, α1=1.0, iter=20,
    c1=1e-4, c2=0.9, lower=-Inf, upper=Inf,
    )
    ϕ(α) = f(clamp.(x + α*p, lower, upper), rot=rot)
    pn = p ./ norm(p)
    ψ(α) = g(clamp.(x + α*p, lower, upper), rot=rot)' * pn
    er = ϕ(0)
    gr = ψ(0)

    for i = 1:iter
        e0 = ϕ(α0)
        e1 = ϕ(α1)
        if (e1 > er + c1 * α1 * gr) || (e1 > e0 && i > 1)
            αstar = zoom(α0, α1, x, p, f, g, rot=rot, iter=20, c1=c1, c2=c2, lower=lower, upper=upper)
            return αstar
        end
        g1 = ψ(α1)
        if abs(g1) <= -c2 * gr
            αstar = α1
            return αstar
        end
        if g1 >= 0
            αstar = zoom(α1, α0, x, p, f, g, rot=rot, iter=20, c1=c1, c2=c2, lower=lower, upper=upper)
            return αstar
        end
        α0 = α1
        α1 = min(1.2*α1, αmax)
    end
    return nothing
end

function zoom(αlo, αhi, x, p, f, g; rot=0.0,
    iter=20, c1=1e-4, c2=0.9,
    lower=-Inf, upper=Inf,
    )

    ϕ(α) = f(clamp.(x + α*p, lower, upper))
    pn = p ./ norm(p)
    ψ(α) = g(clamp.(x + α*p, lower, upper))' * pn
    er = ϕ(0)
    gr = ψ(0)

    for k = 1:iter
        αj = (αhi + αlo)/2
        ej = ϕ(αj)
        elo = ϕ(αlo)
        if (ej > er + c1*αj*gr) || (ej >= elo)
            αhi = αj
        else
            gj = ψ(αj)
            if abs(gj) <= -c2*gr
                αstar = αj
                return αstar
            elseif gj*(αhi - αlo) >= 0
                αhi = αlo
            end
            αlo = αj
        end
    end
    return (αhi + αlo)/2
end


nx = 5
Q = Diagonal(rand(nx))
q = - rand(nx)
c = 1.0
f(x;rot=0) = 0.5 * x'*Q*x + q'*x + c
g(x;rot=0) = Q*x + q
h(x) = Q

x0 = rand(nx)

xsol = bfgs_solve(f, g, x0)
