using Pkg
Pkg.activate(@__DIR__)

using FiniteDiff
using Plots
using LinearAlgebra
using Dojo
using ForwardDiff
using LDLFactorizations
using QDLDL
using BenchmarkTools
using AMD
using MeshCat

vis = Visualizer()
open(vis)

mech = get_sphere(contact_type=:linear)
initialize!(mech, :sphere)
simulate!(mech, 1.0)
M = full_matrix(mech.system)
M[1:6,1:6]
M[1:6, 13:18]
M[13:18, 1:6]


function visualize_particle(x, vis; h=0.01, name=:robot, animation=MeshCat.Animation(Int(ceil(1/h))))
    N = length(x)
    setobject!(vis[name], HyperRectangle(-Vec(0.25,0.25,0.0), Vec(0.5,0.5,0.5)))
    for i = 1:N
        atframe(animation, i) do
            settransform!(vis[name], Translation(x[i]...))
        end
    end
    setanimation!(vis, animation)
    return vis, animation
end

function ip_solver(r, ∇r, c, θ, nd::Int, nc::Int; rtol=1e-6, μ0=1e-3*ones(nc),
        v0=ones(nd+2nc), show_plot=false, I=40)

    plt = plot(xlims=[1e-12,1e2], ylims=[1e-12,1e2], legend=false,
        aspectratio=1.0, axis=:log)

    # Intialisation
    v = v0
    μ = 1.0*ones(nc) # useless value
    iter = 0
    println("\n")
    println("j α       r        sγ       Δs       Δγ       μ")
    for j = 1:I
        plotting(v, μ, plt, show_plot)
        (norm(r(v, μ0, θ), Inf) < rtol) && break
        iter += 1
        g = r(v, μ0, θ)
        H = ∇r(v, 0.0, θ)
        Δaff = get_direction(H, g, v)
        αaff = conesearch(v, Δaff)
        ν, σ = centering(v, nc, αaff, Δaff)
        μ = max.(ν .* σ, μ0)

        # Predictor Corrector
        g = r(v, μ, θ)
        Δ = get_direction(H, g, v)
        α = conesearch(v, Δ)
        α = linesearch(v, Δ, c(v, μ0, θ), c, μ, α, θ)
        v += α .* Δ
        Δz, ΔΓ, ΔS = unpack(Δ)

        # Naive
        # α = linesearch(v, Δaff, c(v, μ0, θ), c, μ, αaff, θ)
        # v += α * Δaff
        # Δz, ΔΓ, ΔS = unpack(Δaff)

        z, Γ, S = unpack(v)
        println("$j" *
            # "$(scn.(α[end-5:end])) " *
            "$(scn(mean(α[end-5:end]))) " *
            "$(scn(norm(g[1:nd+nc],Inf), exp_digits=2)) " *
            "$(scn(norm(S'*Γ,Inf), exp_digits=2)) " *
            "$(scn(norm(ΔS,Inf), exp_digits=2)) " *
            "$(scn(norm(ΔΓ,Inf), exp_digits=2)) " *
            "$(scn(mean(μ)))")
            # "$(scn.(μ))")
    end
    x, S, Γ = unpack(v)
    β = Γ[3:6]
    @show scn.(β)
    println("iter $(iter)  g $(scn(norm(r(v, μ0, θ),Inf)))")
    return v, iter
end

function get_direction(H, g, v)
    x, Γ, S = unpack(v)
    gx, gΓ, gS = unpack(g)
    Hc = compact_hessian(H, v)
    gc = [gx; gΓ - E*gS ./ Γ]
    Δc = - Hc \ gc
    ΔΓ = Δc[nd .+ (1:nc)]
    Δ = [Δc; - (ΔΓ .* S + gS) ./ Γ]
    # Δ = - H \ g
    return Δ
end

function compact_hessian(H, v)
    x, Γ, S = unpack(v)
    Hc = H[1:nd+nc, 1:nd+nc]
    Hc[nd .+ (1:nc), nd .+ (1:nc)] += E*Diagonal(-S ./ Γ)
    return Hc
end

function centering(v, nc, αaff, Δaff)
    # deg = nc
    z, Γ, S = unpack(v)
    # ν = 1/deg * S' * Γ
    ν = 1/1 * S .* Γ
    xaff, Γaff, Saff = unpack(v + αaff .* Δaff)
    # νaff = 1/deg * Saff' * Γaff
    νaff = 1/1 * Saff .* Γaff
    σ = min.(1, max.(0, νaff ./ ν)).^3
    return ν, σ
end

function plotting(v, μ, plt, show_plot)
    z, Γ, S = unpack(v)
    for i = 1:nc
        scatter!(plt, S[i:i], Γ[i:i], color=colors[i], marker=0.2, markersize=12.0,
            markershape=marker_shapes[i])
        annotate!(plt, S[i], Γ[i], string(i), color=:white, seriesalpha=0.4)
    end
    show_plot && display(plt)
end

function linesearch(v, Δ, vio0, c, μ, α, θ)
    for i = 1:20
        vio = c(v + α .* Δ, μ, θ)
        any(vio .< vio0) && break
        α *= 0.5
    end
    return α
end

function conesearch(v, Δ)
    α = ones(nc)
    for j = 1:500
        z, Γ, S = unpack(0.99*v + [ones(nd); α; α] .* Δ)
        for i = 1:nc
            all([S[i]; Γ[i]] .> 1e-10) && break
            α[i] *= 0.9
        end
    end
    return [ones(nd); α; α]
end

function unpack(v)
    off = 0
    z = v[1:nd]; off += nd
    Γ = v[off .+ (1:nc)]; off += nc
    S = v[off .+ (1:nc)]; off += nc
    return z, Γ, S
end

function ϕ(x)
    return x[3:3]
end
N = ForwardDiff.jacobian(ϕ, ones(3))

function velocity_mapping(x, θ)
    x1 = θ[nd .+ (1:nd)]
    v = (x - x1)[1:2]/h
    return [v; -v]
end

function slacks(x, Γ, S, θ)
    # β = Γ[1:4]
    # ψ = Γ[5:5]
    # γ = Γ[6:6]
    # η = S[1:4]
    # ζ = S[5:5]
    # s = S[6:6]
    # [
    #  + η - (velocity_mapping(x, θ) + ψ[1]*ones(4))
    #  - ζ + (cf*γ-[sum(β)]);
    #  + s - ϕ(x);
    #  ]
    # β = Γ[1:4]
    # γ = Γ[5:5]
    # ψ = Γ[6:6]
    # η = S[1:4]
    # s = S[5:5]
    # ζ = S[6:6]
    # [
    #  + η - (velocity_mapping(x, θ) + ψ[1]*ones(4))
    #  + s - ϕ(x);
    #  - ζ + (cf*γ-[sum(β)]);
    #  ]
    γ = Γ[1:1]
    ψ = Γ[2:2]
    β = Γ[3:6]
    s = S[1:1]
    ζ = S[2:2]
    η = S[3:6]
    [
    + s - ϕ(x);
    - ζ + (cf*γ-[sum(β)]);
    + η - (velocity_mapping(x, θ) + ψ[1]*ones(4))
    ]
end

function d(x, γ, β, θ)
    x0 = θ[1:nd]
    x1 = θ[nd .+ (1:nd)]
    u = θ[2nd .+ (1:nd)]
    return m*(x - 2x1 + x0)/h - h*m*g - N'*γ - M'*β - u*h
end

# residual
function r(v, μ, θ)
    z, Γ, S = unpack(v)
    x = z
    # β = Γ[1:4]
    # ψ = Γ[5:5]
    # γ = Γ[6:6]
    # η = S[1:4]
    # ζ = S[5:5]
    # s = S[6:6]
    # β = Γ[1:4]
    # γ = Γ[5:5]
    # ψ = Γ[6:6]
    # η = S[1:4]
    # s = S[5:5]
    # ζ = S[6:6]
    γ = Γ[1:1]
    ψ = Γ[2:2]
    β = Γ[3:6]
    s = S[1:1]
    ζ = S[2:2]
    η = S[3:6]
    return [d(x, γ, β, θ); slacks(x, Γ, S, θ); S .* Γ .- μ]
end

function c(v, μ, θ)
    z, Γ, S = unpack(v)
    x = z
    # β = Γ[1:4]
    # ψ = Γ[5:5]
    # γ = Γ[6:6]
    # η = S[1:4]
    # ζ = S[5:5]
    # s = S[6:6]
    # β = Γ[1:4]
    # γ = Γ[5:5]
    # ψ = Γ[6:6]
    # η = S[1:4]
    # s = S[5:5]
    # ζ = S[6:6]
    γ = Γ[1:1]
    ψ = Γ[2:2]
    β = Γ[3:6]
    s = S[1:1]
    ζ = S[2:2]
    η = S[3:6]
    return norm(d(x, γ, β, θ), Inf), norm(slacks(x, Γ, S, θ), Inf), norm(S .* Γ .- μ, Inf)
end

function ∇r(v, μ, θ)
    # FiniteDiff.finite_difference_jacobian(v -> r(v, μ, θ), v)
    ForwardDiff.jacobian(v -> r(v, μ, θ), v)
end


nd = 3
ncon = 1
nc = 6ncon

cf = 0.3
h = 0.01
M = ForwardDiff.jacobian(x -> velocity_mapping(x, zeros(3nd)), ones(3))
E = ForwardDiff.jacobian(S -> slacks(ones(nd), ones(nc), S, ones(3nd)), ones(nc))
m = 1.0
g = [0.0, 0.0, -10.0]
x0 = [0.0, 0.0, 0.0]
x1 = [0.0, 0.0, 0.0]
u = 0*[100.0, 100.0, 100.0]
θ = [x0; x1; u] # force
colors = [:yellow, :pink, :red, :green, :orange, :magenta]
marker_shapes = [:circle, :rect, :diamond, :hexagon, :utriangle, :dtriangle]

zi = 0.2e+0
μi = 1e-8
rtol0 = 1e-10
S0 = 1e-0*[μi,μi,μi,μi,μi,1]
# S0 = 1e-0*[μi,1,μi,μi,μi,μi]
# S0 = [1,μi,μi,μi,μi,μi]
Γ0 = 1e-0*[1,1,1,1,1,μi]
# Γ0 = 1e-0*[1,μi,1,1,1,1]
# Γ0 = [μi,1,1,1,1,1]
v0 = [zi*ones(nd); 1e-0*Γ0; 1e-0*S0]
ip_solver(r, ∇r, c, θ, nd, nc; rtol=rtol0, μ0=μi*ones(nc), v0=v0, show_plot=false, I=1)

v0 = [zi*ones(nd); 1e-0*S0; 1e-0*Γ0]
ip_solver(r, ∇r, c, θ, nd, nc; rtol=rtol0, μ0=μi*ones(nc), v0=v0, show_plot=true)

v0 = [zi*ones(nd); μi*ones(nc); 1e-0*ones(nc)]
ip_solver(r, ∇r, c, θ, nd, nc; rtol=rtol0, μ0=μi, v0=v0, show_plot=true)

v0 = [zi*ones(nd); μi*ones(nc); μi*ones(nc)]
ip_solver(r, ∇r, c, θ, nd, nc; rtol=rtol0, μ0=μi, v0=v0, show_plot=true)

v0 = [zi*ones(nd); sqrt(μi)*ones(nc); sqrt(μi)*ones(nc)]
ip_solver(r, ∇r, c, θ, nd, nc; rtol=rtol0, μ0=μi, v0=v0, show_plot=true)

v0 = [zi*ones(nd); 1/sqrt(μi)*ones(nc); 1/sqrt(μi)*ones(nc)]
ip_solver(r, ∇r, c, θ, nd, nc; rtol=rtol0, μ0=μi, v0=v0, show_plot=true)


function simulate!(x0, x1, U; rtol=1e-8, μ0=1e-5, warmstart=true, show_plot=false)
    θ = [x0; x1; U[1]]
    v = μ0 * ones(nd + 2nc)
    vs = []
    iters = []
    for i = 1:length(U)-1
        !warmstart && (v = μ0 * ones(nd + 2nc))
        v, iter = ip_solver(r, ∇r, c, θ, nd, nc; rtol=rtol, μ0=μ0, v0=v, show_plot=show_plot)
        x, S, Γ = unpack(v)
        x0 = x1
        x1 = x
        θ = [x0; x1; U[i+1]]
        push!(vs, v)
        push!(iters, iter)
    end
    iters
    return vs, iters
end

cf = 1.10
H = 222
x0 = [0.0, 0.0, 0.1]
x1 = [0.0, 0.01, 0.1]
U = [1*1e-9*[1e3, 1e3, 1e2] .* (rand(3).-0.5) for i=1:H]
ts = [i*h for i = 0:H-2]
# rtol0 = 1e-10
# μ0sim = 1e-8
rtol0 = 1e-6
μ0sim = 1e-4
# rtol0 = 1e-8
# μ0sim = 1e-10
# warm starting works well with rtol >= μ0
# warm starting doesn't work with rtol < μ0
vs_w, iters_w = simulate!(x0, x1, U[1:H], rtol=rtol0, μ0=μ0sim, show_plot=false)
# vs, iters = simulate!(x0, x1, U[1:H], rtol=rtol0, μ0=μ0sim, warmstart=false, show_plot=false)

traj = [v[1:3] for v in vs]
traj_w = [v[1:3] for v in vs_w]
h = 0.01
vis, anim = visualize_particle(traj_w, vis, h=h, name=:robot_w)
# vis, anim = visualize_particle(traj, vis, h=h, name=:robot)


sum(iters[2:end]) / sum(iters_w[2:end])

plt = plot(layout=(3,1), ylims=[0, Inf])
plot!(plt[1], ts, [v[1] for v in vs], color=:black,
    linewidth=3.0, label="x", ylabel="particle altitude", xlabel="time")
plot!(plt[1], ts, [v[2] for v in vs], color=:black,
    linewidth=3.0, label="x", ylabel="particle altitude", xlabel="time")
plot!(plt[1], ts, [v[3] for v in vs], color=:black,
    linewidth=3.0, label="x", ylabel="particle altitude", xlabel="time")
plot!(plt[2], ts, [v[10] for v in vs], color=:black,
    linewidth=3.0, label="γ", ylabel="particle altitude", xlabel="time")
plot!(plt[2], ts, [v[12] for v in vs], color=:black,
    linewidth=3.0, label="γ", ylabel="particle altitude", xlabel="time")
plot!(plt[2], ts, [v[13] for v in vs], color=:black,
    linewidth=3.0, label="γ", ylabel="particle altitude", xlabel="time")
plot!(plt[2], ts, [v[14] for v in vs], color=:black,
    linewidth=3.0, label="γ", ylabel="particle altitude", xlabel="time")
plot!(plt[2], ts, [v[15] for v in vs], color=:black,
    linewidth=3.0, label="γ", ylabel="particle altitude", xlabel="time")
plot!(plt[3], ts, iters, color=:red, linewidth=3.0,
    label="baseline", ylabel="iterations", xlabel="time")
plot!(plt[3], ts, iters_w, color=:blue, linewidth=3.0,
    label="warmstart")


θ_fail = [1.0749936604216475, 0.345451001543372, 1.7595285488739987e-8,
    1.0730488202711659, 0.34350615200672757, 1.1864041044755351e-8,
    -39.31172751966181, 0.7359591976611735, 32.736306617565916];
v_fail = [1.0730488202711659, 0.34350615200672757, 1.1864041044755351e-8,
    1.1864041044566393e-8, 2.08743427470252e-6, 2.462924699301821e-7,
    5.589014923184363e-11, 6.336856929609757e-8, 3.0960514914330893e-7,
    0.7667077899163286, 1.548305196669243e-7, 1.04979274779736e-6,
    1.1982359485922102e-8, 0.22199137332345495, 0.008017814574647229];
# v_fail[end-11:end] = max.(v_fail[end-11:end], 1e-8)
rtol0 = 1e-10
μ0fail = 1e-8
rtol0 = 1e-6
μ0fail = 1e-4
ip_solver(r, ∇r, c, θ_fail, nd, nc; rtol=rtol0, μ0=μ0fail, v0=v_fail, show_plot=true, I=100)
ip_solver(r, ∇r, c, θ_fail, nd, nc; rtol=rtol0, μ0=μ0fail, show_plot=true)

plt = plot(xlims=[1e-12,1e2], ylims=[1e-12,1e2], legend=false, aspectratio=1.0, axis=:log)
plotting(v_fail, μ0fail, plt, true)








v0 = 0.1*ones(15)
H0 = ∇r(v0, μ0sim, [x0; x1; U[1]])

impact = [1,2,3,4,10]
H0[impact, impact]
Hc0 = compact_hessian(H0, v0)
Hc0 - Hc0'

plot(Gray.(1e9abs.(H0)))
plot(Gray.(1e9abs.(H0+H0')))
plot(Gray.(1e9abs.(H0-H0')))
H0[1:3,1:3]
H0[4:9,4:9]
H0[4:9,end-5:end]
H0[end-5:end,4:9]
H0[end-5:end,end-5:end]
H0[4:9,1:3]
H0[1:3,4:9]'

Hc0 = compact_hessian(H0, rand(15))
plot(Gray.(1e9abs.(Hc0)))
plot(Gray.(1e9abs.(Hc0 + Hc0')))
plot(Gray.(1e9abs.(Hc0 - Hc0')))
Hc0[1:3,1:3]
Hc0[1:3,4:9]
Hc0[4:9,1:3]
Hc0[4:9,4:9]
Hc0 - Hc0'
function permutation_matrix(p)
    n = length(p)
    P = zeros(n,n)
    for i = 1:n
        P[p[i],i] = 1
    end
    return P
end

function permute_matrix(A, p)
    P = permutation_matrix(p)
    return P' * A * P
end

P = permutation_matrix([1,2,3,6,7,8,9,4,5])
Hp0 = P' * Hc0 * P
plot(Gray.(1e9abs.(Hp0 - Hp0')))




# L = 500
# A = sprand(L,L,5/L) + I
# A = A * A'
# # A = sparse(Diagonal([1,2,3,4,-0.0]))
# b = sprand(L,0.5)
# LDLT = ldl(A)  # LDLᵀ factorization of A
#
# @benchmark x = LDLT \ b   # solves Ax = b
#
# F = qdldl(A)
# x = solve(F, b)
# @benchmark solve!(F, b)
#


A = sprand(100, 100, 1/200)
A = A * A'
plot(Gray.(100*abs.(Matrix(A))))
p_amd = amd(A)
p_symamd = symamd(A)
p_colamd = colamd(A)

Ap = permute_matrix(A, p_symamd)
Ap = permute_matrix(A, p_amd)
plot(Gray.(100*abs.(Matrix(A))))
plot(Gray.(100*abs.(Matrix(Ap))))

p_amd = amd(sparse(Hc0))
p_symamd = symamd(sparse(Hc0))
p_colamd = colamd(sparse(Hc0))
Hc0_a = permute_matrix(Hc0, p_amd)
Hc0_s = permute_matrix(Hc0, p_symamd)
Hc0_c = permute_matrix(Hc0, p_colamd)
plot(Gray.(100*abs.(Matrix(Hc0))))
plot(Gray.(100*abs.(Matrix(Hc0_a))))
plot(Gray.(100*abs.(Matrix(Hc0_s))))
plot(Gray.(100*abs.(Matrix(Hc0_c))))