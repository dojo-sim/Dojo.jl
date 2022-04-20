function ip_solver(d, ϕ, θ, nd::Int, nϕ::Int; rtol=1e-6, μ0=1e-3, v0=ones(nd+2nϕ), show_plot=false)
    plt = plot(xlims=[1e-12,1e2], ylims=[1e-12,1e2], legend=false, aspectratio=1.0,
        axis=:log)

    # residual
    function r(v, μ)
        x, s, γ = unpack(v)
        return [d(x, γ, θ); s - ϕ(x); s .* γ .- μ]
    end
    function c(v, μ)
        x, s, γ = unpack(v)
        return norm(d(x, γ, θ), Inf), norm([s - ϕ(x); s .* γ .- μ], Inf)
    end
    function ∇r(v, μ)
        FiniteDiff.finite_difference_jacobian(v -> r(v, μ), v)
    end

    # Intialisation
    x, s, γ = unpack(v0)
    v = [x; s; γ]
    μ = 1.0 # useless value
    iter = 0
    for j = 1:40
        plotting(v, μ, plt, show_plot)
        (norm(r(v, μ0), Inf) < rtol) && break
        iter += 1
        g = r(v, μ0)
        H = ∇r(v, 0.0)
        Δaff = - H \ g
        αaff = conesearch(v, Δaff)
        ν, σ = centering(v, nϕ, αaff, Δaff)
        μ = max(ν*σ, μ0)
        g = r(v, μ)
        Δ = - H \ g
        α = conesearch(v, Δ)
        α = linesearch(v, Δ, c(v,μ0), c, μ, α)
        v += α * Δ
        Δx, Δs, Δγ = unpack(Δ)
        x, s, γ = unpack(v)
        println("j $j  α$(scn(α))  r$(scn(norm(g[1:nd+nϕ],Inf)))  sγ$(scn(norm(s'*γ,Inf)))  Δs$(scn(norm(Δs,Inf)))  Δγ$(scn(norm(Δγ,Inf)))  μ$(scn(μ))")
    end
    println("iter $(iter)  g $(scn(norm(r(v,μ0),Inf)))")
    x, s, γ = unpack(v)
    return x, s, γ, iter
end

function centering(v, nϕ, αaff, Δaff)
    deg = nϕ
    x, s, γ = unpack(v)
    ν = 1/deg * s' * γ
    xaff, saff, γaff = unpack(v + αaff * Δaff)
    νaff = 1/deg * saff' * γaff
    σ = min(1, max(0, νaff/ν))^3
    return ν, σ
end

function plotting(v, μ, plt, show_plot)
    x, s, γ = unpack(v)
    scatter!(plt, s, γ, color=RGB(0.0,0.0,0.0), α=0.2, markersize=8.0)
    annotate!(plt, s[1], γ[1], string(Int(round(log(10, μ)))), color=:blue)
    show_plot && display(plt)
end

function linesearch(v, Δ, vio0, c, μ, α)
    for i = 1:20
        vio = c(v + α * Δ, μ)
        any(vio .< vio0) && break
        α *= 0.5
    end
    return α
end

function conesearch(v, Δ)
    α = 1.0
    for i = 1:500
        x, s, γ = unpack(v + α * Δ)
        all([s; γ] .> 0.0) && break
        α *= 0.9
    end
    return α
end


nd = 1
nϕ = 1

h = 0.02
m = 1.0
g = [-10.0]
x0 = [0.0]
x1 = [0.0]
u = [20.0]
θ = [x0; x1; u] # force

function d(x, γ, θ)
    x0 = θ[1:1]
    x1 = θ[2:2]
    u = θ[3:3]
    return m*(x - 2x1 + x0)/h - h*m*g - γ - u*h
end

function ϕ(x)
    return x
end

function unpack(v)
    x = v[1:nd]
    s = v[nd .+ (1:nϕ)]
    γ = v[nd + nϕ .+ (1:nϕ)]
    return x, s, γ
end

xi = 2.2e+0
μi = 1e-5
v0 = [xi*ones(nd); 1e-0*ones(nϕ); 1e-0*ones(nϕ)]
ip_solver(d, ϕ, θ, nd, nϕ; rtol=1e-8, μ0=μi, v0=v0, show_plot=true)

v0 = [xi*ones(nd); 1e-0*ones(nϕ); μi*ones(nϕ)]
ip_solver(d, ϕ, θ, nd, nϕ; rtol=1e-8, μ0=μi, v0=v0, show_plot=true)

v0 = [xi*ones(nd); μi*ones(nϕ); 1e-0*ones(nϕ)]
ip_solver(d, ϕ, θ, nd, nϕ; rtol=1e-8, μ0=μi, v0=v0, show_plot=true)

v0 = [xi*ones(nd); μi*ones(nϕ); μi*ones(nϕ)]
ip_solver(d, ϕ, θ, nd, nϕ; rtol=1e-8, μ0=μi, v0=v0, show_plot=true)

v0 = [xi*ones(nd); sqrt(μi)*ones(nϕ); sqrt(μi)*ones(nϕ)]
ip_solver(d, ϕ, θ, nd, nϕ; rtol=1e-8, μ0=μi, v0=v0, show_plot=true)

v0 = [xi*ones(nd); 1/sqrt(μi)*ones(nϕ); 1/sqrt(μi)*ones(nϕ)]
ip_solver(d, ϕ, θ, nd, nϕ; rtol=1e-8, μ0=μi, v0=v0, show_plot=true)


function simulate!(x0, x1, U; rtol=1e-8, μ0=1e-5, warmstart=true, show_plot=false)
    θ = [x0; x1; U[1]]
    v = μ0 * ones(nd + 2nϕ)
    vs = []
    iters = []
    for i = 1:length(U)-1
        !warmstart && (v = μ0 * ones(nd + 2nϕ))
        x, s, γ, iter = ip_solver(d, ϕ, θ, nd, nϕ; rtol=rtol, μ0=μ0, v0=v, show_plot=show_plot)
        v = [x; s; γ]
        x0 = x1
        x1 = x
        θ = [x0; x1; U[i+1]]
        push!(vs, v)
        push!(iters, iter)
    end
    iters
    return vs, iters
end

N = 100
x0 = [1.0]
x1 = [1.0]
U = [100*(rand(1).-0.5) for i=1:N]
ts = [i*h for i = 0:N-2]
vs_w, iters_w = simulate!(x0, x1, U, rtol=1e-8, μ0=1e-8)
vs, iters = simulate!(x0, x1, U, rtol=1e-8, μ0=1e-5, warmstart=false)

sum(iters) / sum(iters_w)

plt = plot(layout=(3,1), ylims=[0, Inf])
plot!(plt[1], ts, [v[1] for v in vs], color=:black, linewidth=3.0, label="x", ylabel="particle altitude", xlabel="time")
plot!(plt[2], ts, [v[3] for v in vs], color=:black, linewidth=3.0, label="γ", ylabel="particle altitude", xlabel="time")
plot!(plt[3], ts, iters, color=:red, linewidth=3.0, label="baseline", ylabel="iterations", xlabel="time")
plot!(plt[3], ts, iters_w, color=:blue, linewidth=3.0, label="warmstart")




us = 0.0:0.15:20.0
θs = [[u] for u in us]

plt = plot()
for μ0 in [1e-6, 1e-5, 1e-4, 1e-3, 1e-2]#, 1e-1, 1e-0]
    # xs = [ip_solver(d, ϕ, θ, nd, nϕ, rtol=1e-12, μ0=μ0)[1][1] for θ in θs]
    plot!(plt, us, xs)
end
display(plt)


mech = get_sphere(gravity = 0.00, contact_type=:linear)
initialize!(mech, :sphere, x=[0,0,1.0])
simulate!(mech, 1.0)
M = full_matrix(mech.system)
M[1:6,1:6]
M[1:6,7:12]
M[1:6,13:18]'
M[7:12,1:6]


using QDLDL
using SparseArrays
n = 5
N = 1.0*rand(n,n)
A1 = rand(n,n)
A1 = A1'A1
A2 = rand(n,n)
A2 = A2'A2
A = sparse([A1 N'
            N  A2])
qdldl(A)
eigvals(Matrix(A))

M[13:18,1:6]
M[7:12,7:12]
M[7:12,13:18]
M[13:18,7:12]
M[13:18,13:18]
