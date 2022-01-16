# Utils
function module_dir()
    return joinpath(@__DIR__, "..", "..")
end

# Activate package
using Pkg
Pkg.activate(module_dir())

using MeshCat

# Open visualizer
vis = Visualizer()
open(vis)

# Include new files
include(joinpath(module_dir(), "examples", "loader.jl"))


include("methods.jl")

conetype = :linear
# conetype = :soc
mech = getmechanism(:box2d, Δt=0.05, g=-10.0, cf=1.0, contact=true, conetype=conetype);
btol = 1e-4
Fs = 0:0.5:20
Fsref = 0:0.05:20
plt = plot(layout=(3,1), size=(500,800), legend=:topleft)
# mode = :friction
mode = :impact
undercut = (mode == :impact) ? 1.01 : 1.50
if mode == :friction && conetype == :soc
    filename = "nonlinear_friction"
elseif mode == :friction && conetype == :linear
    filename = "linear_friction"
else
    filename = "impact"
end

for btol = 1e-10
    res = [box2d_dojo(mech, F, btol=btol, mode=mode) for F = Fsref]
    x = [r[1] for r in res]
    ∇x = [r[2] for r in res]
    plot!(plt[1,1], Fsref,  x, xlabel="F", ylabel="x", linewidth=2.0,
        label="true dynamics")
    plot!(plt[2,1], Fsref, ∇x, xlabel="F", ylabel="∇x", linewidth=2.0,
        label="true dynamics", ylims=(-0.1, 1.5))
    plot!(plt[3,1], Fsref, ∇x, xlabel="F", ylabel="∇x", linewidth=2.0,
        label="true dynamics", ylims=(-0.1, 1.5))
    display(plt)
end

# for btol in exp.(log(10)* (-10:1:-3))
for btol in [1e-3, 1e-4, 1e-5, 1e-6, 1e-8]
# for btol = 3e-4
    res = [box2d_dojo(mech, F, btol=btol, mode=mode, undercut=undercut) for F = Fsref]
    write_csv(["F", "x", "grad"], [[Fsref[i]; res[i]...] for i=1:length(Fsref)], "data/$(filename)_dojo_"*scn(btol)[2:end]*".csv")
    x = [r[1] for r in res]
    ∇x = [r[2] for r in res]
    # plot!(plt[1,1], Fsref,  x, xlabel="F", ylabel="x", linewidth=3.0, linestyle=:dot,
        # label="Dojo btol="*scn(btol))
    plot!(plt[2,1], Fsref, ∇x, xlabel="F", ylabel="∇x", linewidth=3.0, linestyle=:dot,
        label="Dojo btol="*scn(btol)[2:end])
    display(plt)
end

for Σ in [1e-1, 1e-2, 1e-3]
# for Σ in [4e-2]
    res = [box2d_gradientbundle(mech, F, N=100, Σ=Σ*I, mode=mode) for F = Fs]
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
display(plt)


δu1 = [rand(3) for k=1:10]
B1 = rand(6,3)
Bf1 = vec(B1)
δx1 = [B1*δu for δu in δu1]
δx1 = [rand(6) for k=1:10]

leastsquares(δx1, δu1)
