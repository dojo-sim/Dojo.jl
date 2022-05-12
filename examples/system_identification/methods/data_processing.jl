
function newton_solve(f, g, H, x0; ftol=-Inf, gtol=1e-4, iter=100,
        lower=-Inf, upper=Inf, reg = 1e-3)
    x = copy(x0)
    X = [copy(x0)]
    for k = 1:iter
		fe = f(x)
		ge = g(x)
		He = H(x)
        He += reg * norm(He, Inf) * I
        ((norm(ge, Inf) < gtol) || (fe < ftol)) && break
        p = - He \ ge
        α = clamped_linesearch(f, x, p, lower=lower, upper=upper)
        x = clamp.(x + α*p, lower, upper)
        push!(X, copy(x))
        println("k:", k,
            "   f:", Dojo.scn(fe, digits=3),
            "   ∇:", Dojo.scn.(norm(ge, Inf)),
            "   α:", Dojo.scn.(α),
            # "   g:", scn.(-g(x, rot=rot)),
            # "   p:", scn.(p),
            # "   x:", scn.(x),
            # "   H:", scn.(H),
            )
    end
    return x, X
end

function clamped_linesearch(f, x, p; iter=10,
        lower=-Inf, upper=Inf)
    α = 1.0
    fprev = f(x)
    for k = 1:iter
        xc = clamp.(x + α*p, lower, upper)
        (f(xc) <= fprev) && break
        α /= 2
    end
    return α
end

################################################################################
# Generate & Save Dataset
################################################################################
generate_hardware_dataset(N=50, sleep_ratio=0.01)

################################################################################
# Load Dataset
################################################################################
params0, trajs0, pairs0 = open_dataset(:hardwarebox; N=50)

pairs0
params0
data0 = params0[:data]


function quad_loss(x, timestep, x0, v0)
	loss = 0.0
	for k = 1:7
		t = (k-1) * timestep
		Δx = x[k][3] - (x0 + v0*t - t^2/2 * 9.81)
		loss += Δx^2
	end
	return loss
end

function pairs_loss(pairs, θ)
	loss = 0.0
	for i = 1:50
		ps = pairs[(i-1)*6 .+ (1:6)]
		x = [[ps[1][1][1:3]]; [p[2][1:3] for p in ps]]
		loss += quad_loss(x, θ[1], θ[1 + (i-1)*2 .+ (1:2)]...)
	end
	return loss / 50
end

function pairs_grad(pairs, θ)
	grad = FiniteDiff.finite_difference_gradient(θ -> pairs_loss(pairs, θ), θ)
	return grad
end

function pairs_hess(pairs, θ)
	hess = FiniteDiff.finite_difference_hessian(θ -> pairs_loss(pairs, θ), θ)
	return hess
end


θ0 = ones(101)
pairs_loss(pairs0, θ0)
pairs_grad(pairs0, θ0)
plot(Gray.(abs.(pairs_hess(pairs0, θ0))))

f00(θ) = pairs_loss(pairs0, θ)
g00(θ) = pairs_grad(pairs0, θ)
H00(θ) = pairs_hess(pairs0, θ)


θsol, θSOL = newton_solve(f00, g00, H00, θ0, gtol=1e-6, reg = 0.)
θsol
