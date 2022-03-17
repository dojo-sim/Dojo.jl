using LinearAlgebra
using Plots
using FiniteDiff


function pack_vars(p3, γ)
    return [p3; γ]
end

function unpack_vars(vars)
    p3 = vars[1:3]
    γ = vars[4:4]
    return p3, γ
end

function pack_data(p1, p2, ρ)
    return [p1; p2; ρ]
end

function unpack_data(data)
    p1 = data[1:3]
    p2 = data[4:6]
    ρ = data[7]
    return p1, p2, ρ
end

function residual(vars, data)
    p1, p2, ρ = unpack_data(data)
    p3, γ = unpack_vars(vars)
    rdyn = m*(p3 - 2p2 + p1)/h - h*m*[0,0,g] + [0 0 1]'*γ
    rimp = γ - 1/ρ * min.(p3[3:3], 0)
    return [rdyn; rimp]
end

function jacobian_residual(vars, data)
    FiniteDiff.finite_difference_jacobian(vars -> residual(vars,data), vars)
end

function linesearch(vars, data, Δ, r)
    α = 1.0
    for i = 1:10
        vars_cand = vars + α*Δ
        r_cand = residual(vars_cand, data)
        (norm(r_cand, Inf) <= norm(r, Inf)) && break
        α = 1.0
    end
    return α
end

function soft_step(vars, data)
    for k = 1:2
        for i = 1:10
            r = residual(vars, data)
            (norm(r, Inf) < 1e-6) && break
            J = jacobian_residual(vars, data)
            @show J
            Δ = - J \ r
            α = linesearch(vars, data, Δ, r)
            vars += α * Δ
        end
        p1, p2, ρ = unpack_data(data)
        ρ /= 10.0
        data = pack_data(p1, p2, ρ)
    end
    return unpack_vars(vars)
end

function soft_simulate(vars, data, N)
    p1, p2, ρ = unpack_data(data)
    P = [p1, p2]
    Γ = []
    for t = 1:N
        p3, γ = soft_step(vars, data)
        p1, p2, ρ = unpack_data(data)
        data = pack_data(p2, p3, 1.0)
        vars = pack_vars(p3, zeros(nγ))
        push!(P, p3)
        push!(Γ, γ)
    end
    return P, Γ
end


# problem parameters
m = 1.0
g = -9.81
h = 0.01

p10 = [0,0,1.0]
p20 = [0,0,1.0]
ρ0 = 1.0
data0 = pack_data(p10, p20, ρ0)

p30 = [0,0,1.0]
γ0 = [0.0]
nγ = 1
vars0 = pack_vars(p30, γ0)

residual(vars0, data0)
jacobian_residual(vars0, data0)
soft_step(vars0, data0)
P0, Γ0 = soft_simulate(vars0, data0, 100)

plot([p[3] for p in P0], linewidth=3.0, xlabel="time step", ylabel="altitude")
P3 = [p[3] for p in P0]
plot(log.(10, abs.(P3)))



X1 = -5:0.05:0
plot(X1, 0.5*(X1 .- 1).^2)
X2 = 0:0.05:10
plot(X2, -log.(X2 .+ 1))
plot([X1;X2], [0.5*((X1 .- 1).^2 .- 1); -log.(X2 .+ 1)])


function soft_potential(y, ρ)
    if y > 0
        p = -log(y+ρ) + log(ρ)
    else
        p = 1/(2ρ)*(y-1)^2 - 1/(2ρ)
    end
    return p
end

function gradient_soft_potential(y, ρ)
    if y > 0
        g = -1/(y+ρ)
    else
        g = 1/ρ * (y-1)
    end
    return g
end



plt = plot()
Y = 0:0.01:10
for (i,ρ) in enumerate([1, 0.1, 0.01])#, 1e-2, 1e-3, 1e-4, 1e-5]
    # plot!(plt, Y, log.(10, -soft_potential.(Y, ρ)), linewidth=2i)
    plot!(plt, Y, soft_potential.(Y, ρ), linewidth=2i)
end
display(plt)
