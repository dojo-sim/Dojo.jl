function ip_solver(d, ϕ, θ, nd::Int, nϕ::Int; rtol=1e-6, μ0=1e-3, v0=ones(nd+2nϕ), show_plot=false)
    plt = plot(xlims=[1e-8,10.0], ylims=[1e-8,10.0], legend=false, aspectratio=1.0, axis=:log)

    # residual
    function r(v, μ)
        x, s, γ = unpack(v)
        return [d(x, γ, θ); s - ϕ(x); s .* γ .- μ]
    end
    function ∇r(v, μ)
        FiniteDiff.finite_difference_jacobian(v -> r(v, μ), v)
    end

    # Intialisation
    x0, s0, γ0 = unpack(v0)
    μ = max(μ0, s0' * γ0)
    v = v0
    iter = 0
    for i = 1:10
        (norm(r(v, μ0), Inf) < rtol) && break
        for j = 1:10
            iter += 1
            x, s, γ = unpack(v)
            scatter!(plt, s, γ, color=RGB(0.0,0.0,0.0), α=0.2, markersize=8.0)
            annotate!(plt, s[1], γ[1], string(Int(round(log(10, μ)))), color=:blue)
            show_plot && display(plt)
            g = r(v, μ)
            (norm(g, Inf) < rtol) && break
            H = ∇r(v, μ)
            Δ = - H \ g
            α = linesearch(v, Δ, g, r, μ)
            v += α * Δ
            println("i $i   j $j  α $(scn(α)) g $(scn(norm(g,Inf))) Δ $(scn(norm(Δ,Inf)))")
        end
        μ = max(μ0, μ/10.0)
    end
    println("iter $(iter)  g $(scn(norm(r(v,μ0),Inf)))")
    return unpack(v)
end


function linesearch(v, Δ, g0, r, μ)
    α = 1.0
    for i = 1:20
        x, s, γ = unpack(v + α * Δ)
        all([s; γ] .> 0.0) && break
        α *= 0.5
    end
    for i = 1:20
        g = r(v + α * Δ, μ)
        (norm(g, Inf) <= norm(g0, Inf)) && break
        α *= 0.5
    end
    return α
end
plt = plot(1:0.01:2, 1:0.01:2)
annotate!(plt, 1, 1, "my text", :red)

nd = 1
nϕ = 1

h = 0.1
m = 1.0
g = [-10.0]
x0 = [0.0]
x1 = [0.0]
θ = [4.0] # force

function d(x, γ, θ)
    return m*(x - 2x1 + x0)/h - h*m*g - γ - θ*h
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

v0 = [zeros(nd); 1e-0*ones(nϕ); 1e-0*ones(nϕ)]
ip_solver(d, ϕ, θ, nd, nϕ; rtol=1e-10, μ0=1e-8, v0=ones(nd+2nϕ), show_plot=true)
v0 = [zeros(nd); 1e-0*ones(nϕ); 1e-8*ones(nϕ)]
ip_solver(d, ϕ, θ, nd, nϕ; rtol=1e-10, μ0=1e-8, v0=v0, show_plot=true)
v0 = [zeros(nd); 1e-8*ones(nϕ); 1e-0*ones(nϕ)]
ip_solver(d, ϕ, θ, nd, nϕ; rtol=1e-10, μ0=1e-8, v0=v0, show_plot=true)
v0 = [zeros(nd); 1e-8*ones(nϕ); 1e-8*ones(nϕ)]
ip_solver(d, ϕ, θ, nd, nϕ; rtol=1e-10, μ0=1e-8, v0=v0, show_plot=true)














us = 0.0:0.15:20.0
θs = [[u] for u in us]

plt = plot()
for μ0 in [1e-6, 1e-5, 1e-4, 1e-3, 1e-2]#, 1e-1, 1e-0]
    # xs = [ip_solver(d, ϕ, θ, nd, nϕ, rtol=1e-12, μ0=μ0)[1][1] for θ in θs]
    plot!(plt, us, xs)
end
display(plt)
