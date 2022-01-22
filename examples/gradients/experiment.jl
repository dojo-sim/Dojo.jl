using Dojo
using Plots

include("methods.jl")

# cone
contact_type = :linear_contact
# contact_type = :contact

# scale system for nice plots
mech = getmechanism(:box2d, Δt=0.1, g=-1.0, cf=1.0, contact=true, contact_type=contact_type);

# inputs
Fs = 0.01:0.01:0.2
Fsref = 0.001:0.001:0.2

# options
undercut = 1.0

mode = :friction
# mode = :impact
if mode == :friction && contact_type == :contact
    filename = "nonlinear_friction"
elseif mode == :friction && contact_type == :linear_contact
    filename = "linear_friction"
else
    filename = "impact"
end


plt = plot(layout=(3,1), size=(500,800), legend=:topleft)

# hard contact
for btol = 1e-10
    res = [box2d_dojo(mech, F, btol=btol, mode=mode) for F = Fsref]
    write_csv(["F", "x", "grad"], [[Fsref[i]; res[i]...] for i=1:length(Fsref)], "data/$(filename)_dojo_"*scn(btol)[2:end]*".csv")
    x = [r[1] for r in res]
    ∇x = [r[2] for r in res]
    plot!(plt[1,1], Fsref,  x, xlabel="F", ylabel="x", linewidth=2.0,
        label="true dynamics", color=:black)
    plot!(plt[2,1], Fsref, ∇x, xlabel="F", ylabel="∇x", linewidth=2.0,
        label="true dynamics", color=:black)#ylims=(-0.1, 1.5))
    plot!(plt[3,1], Fsref, ∇x, xlabel="F", ylabel="∇x", linewidth=2.0,
        label="true dynamics", color=:black)#, ylims=(-0.1, 1.5))
    display(plt)
end

# smooth contact
for btol in [1e-4, 1e-5, 1e-6, 1e-7, 1e-8] .* undercut
    res = [box2d_dojo(mech, F, btol=btol, mode=mode, undercut=undercut) for F = Fsref]
    write_csv(["F", "x", "grad"], [[Fsref[i]; res[i]...] for i=1:length(Fsref)], "data/$(filename)_dojo_"*scn(btol)[2:end]*".csv")
    x = [r[1] for r in res]
    ∇x = [r[2] for r in res]
    plot!(plt[1,1], Fsref,  x, xlabel="F", ylabel="x", linewidth=3.0, linestyle=:dot,
        label="Dojo btol="*scn(btol))
    plot!(plt[2,1], Fsref, ∇x, xlabel="F", ylabel="∇x", linewidth=3.0, linestyle=:dot,
        label="Dojo btol="*scn(btol)[2:end])
    display(plt)
end

# gradient bundles
for Σ in [1e-2, 1e-3, 1e-4]
    res = [box2d_gradientbundle(mech, F, N=500, Σ=Σ*I, mode=mode) for F = Fs]
    write_csv(["F", "x", "gb0", "gb1"], [[Fs[i]; res[i]...] for i=1:length(Fs)], "data/$(filename)_gb_"*scn(Σ)[2:end]*".csv")
    x = [r[1] for r in res]
    ∇x0 = [r[2] for r in res]
    ∇x1 = [r[3] for r in res]
    # plot!(plt[1,1], Fs,  x, xlabel="F", ylabel="x", linewidth=3.0, linestyle=:dot)
    plot!(plt[3,1], Fs, ∇x0, xlabel="F", ylabel="∇x", linewidth=6.0,
        linestyle=:dot,  label="0th order G.B. Σ="*scn(Σ)[2:end])
    plot!(plt[3,1], Fs, ∇x1, xlabel="F", ylabel="∇x", linewidth=3.0,
        linestyle=:dash, label="1st order G.B. Σ="*scn(Σ)[2:end])
    display(plt)
end

################

################
δu1 = [rand(3) for k=1:10]
B1 = rand(6,3)
Bf1 = vec(B1)
δx1 = [B1*δu for δu in δu1]
δx1 = [rand(6) for k=1:10]

leastsquares(δx1, δu1)





using LinearAlgebra
using Symbolics
################################################################################
# Impact simplified
################################################################################
m = 1.0
Δt = 0.1
gv = 10.0
ρ = 1e-4

function residual(z3, γ, u, ρ)
    z0 = ρ / (m*gv)
    r = [
        m * (z0 - z3)/Δt + Δt*(γ + u - m*gv),
        z3 * γ - ρ,
        ]
    return r
end

function analytical_root(u, ρ)
    z0 = ρ / (m*gv)
    a = - m / Δt
    b = Δt * (u - m*gv) + m * z0 / Δt
    c = Δt * ρ
    Δ = b^2 - 4a*c
    z3 = (-b - sqrt(Δ)) / (2a)
    γ = ρ / z3
    return z3, γ
end

@variables z3_, γ_, u_, ρ_
r_ = residual(z3_, γ_, u_, ρ_)
z3s_, γs_ = analytical_root(u_, ρ_)
∇x = Symbolics.jacobian(r_, [z3_, γ_])
∇u = Symbolics.jacobian(r_, [u_])
dz3du = (- inv(∇x) * ∇u)[1]
dz3du_fct = eval(build_function(dz3du, z3_, γ_))
gradmap = dz3du_fct(analytical_root(u_, ρ_)...)
gradmap_fct = eval(build_function(gradmap, u_, ρ_))

plot(gradmap_fct.(U, 1e-6))

using Plots; pyplot()
y=range(0.001,stop=20.0,length=200)
x=range(-15,stop=3,length=100)
f(x,y) = gradmap_fct(y, exp(log(10) * x))
plot(x,y,f,st=:surface,camera=(0,30))

z30 = z0
γ0 = m*gv
u0 = 0.0
ρ0 = ρ
residual(z30, γ0, u0, ρ0)
z30s, γ0s = analytical_root(u0, ρ0)
residual(z30s, γ0s, u0, ρ0)

U = [u for u = 0:0.01:20]
Z3 = [analytical_root(u, 1e-6)[1] for u in U]
dZ3dU = [dz3du_fct(z3,1e-6/z3) for z3 in Z3]
plot(U, Z3)
plot(U, dZ3dU)



################################################################################
# Linear Friction simplified
################################################################################
m = 1.0
Δt = 0.1
μγ = 10.0
ρψ = 1e-4
ρβ = 1e-4*ones(2)

function residual(x, β, sψ, η, ψ, u, ρψ, ρβ)
    r = [
        m*x/Δt + Δt*(u+β[1]-β[2]);
        sψ - (μγ - β[1] - β[2]);
        η - [x, -x] ./Δt - ψ*ones(2);
        sψ * ψ - ρψ;
        η .* β .- ρβ;
        ]
    return r
end

function optim_root(u, ρψ, ρβ, tol=1e-8)
    z = ones(7)
    for k = 1:40
        r = residual(unpack_z(z)..., u, ρψ, ρβ)
        (norm(r, Inf) < tol) && break
        Δz = - jacz_fct(z, u, ρψ, ρβ) \ r
        x, β, sψ, η, ψ = unpack_z(z)
        Δx, Δβ, Δsψ, Δη, Δψ = unpack_z(Δz)
        αψ = ort_step_length([sψ; ψ], [Δsψ; Δψ], τ=0.90)
        αβ = ort_step_length([η; β], [Δη; Δβ], τ=0.90)
        z += min(αψ, αβ) * Δz
        (k == 40) && (@show "failure")
    end
    return unpack_z(z)
end

function unpack_z(z)
    x = z[1]
    β = z[2:3]
    sψ = z[4]
    η = z[5:6]
    ψ = z[7]
    return x, β, sψ, η, ψ
end
function pack_z(x, β, sψ, η, ψ)
    return [x; β; sψ; η; ψ]
end

u0 = 1.0
ρψ0 = 1e-4
ρβ0 = 1e-6
optim_root(u0, ρψ0, ρβ0)

@variables x_, β_[1:2], sψ_, η_[1:2], ψ_, u_, ρψ_, ρβ_[1:2]
z_ = [x_; β_; sψ_; η_; ψ_]
r_ = residual(x_, β_, sψ_, η_, ψ_, u_, ρψ_, ρβ_)
jacz_ = Symbolics.jacobian(r_, [x_; β_; sψ_; η_; ψ_])
jacθ_ = Symbolics.jacobian(r_, [u_; ρψ_; ρβ_])

jacz_fct = eval(build_function(jacz_, z_, u_, ρψ_, ρβ_)[1])

dxdu_fct = eval(build_function((jacz_ \ jacθ_)[1,1], z_, u_, ρψ_, ρβ_))
function dxdu_map(u, ρψ, ρβ)
    z = pack_z(optim_root(u, ρψ, ρβ)...)
    return dxdu_fct(z, u, ρψ, ρβ)
end

U = [u for u = 0:0.01:20]
plot(U, [dxdu_map(u, 1e-4, 1e-4*ones(2)) for u in U])
dxdu_map(U[1], 1e-6, 1e-6*ones(2))


using Plots; pyplot()
y=range(0.001,stop=20.0,length=200)
x=range(-15,stop=-1,length=100)
f(x,y) = gradmap_fct(y, exp(log(10) * x))
plot(x,y,f,st=:surface,camera=(90,30))

z30 = z0
γ0 = m*gv
u0 = 0.0
ρ0 = ρ
residual(z30, γ0, u0, ρ0)
z30s, γ0s = analytical_root(u0, ρ0)
residual(z30s, γ0s, u0, ρ0)

Z3 = [analytical_root(u, 1e-6)[1] for u in U]
dZ3dU = [dz3du_fct(z3,1e-6/z3) for z3 in Z3]
plot(U, Z3)
plot(U, dZ3dU)
UnderlyingUnderlying
