using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

# ## Setup
using Dojo 
using Plots
include("methods/gradient_bundles.jl")
include("methods/least_squares.jl")
include("methods/utilities.jl")

# ## cone
contact_type = :linear
## contact_type = :nonlinear

# ## scale system for nice plots
mech = get_mechanism(:block2d, 
    timestep=0.1, 
    g=-1.0, 
    friction_coefficient=1.0, 
    contact=true, 
    contact_type=contact_type);

# ## inputs
Fs = 0.01:0.01:0.2
Fsref = 0.001:0.001:0.2

# ## options
undercut = 1.0

mode = :friction
## mode = :impact
if mode == :friction && contact_type == :nonlinear
    filename = "nonlinear_friction"
elseif mode == :friction && contact_type == :linear
    filename = "linear_friction"
else
    filename = "impact"
end

plt = plot(layout=(3,1), size=(500,800), legend=:topleft)

# ## hard contact
for btol = 1e-10
    res = [block2d_dojo(mech, F, btol=btol, mode=mode) for F = Fsref]
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

# ## smooth contact
for btol in [1e-4, 1e-5, 1e-6, 1e-7, 1e-8] .* undercut
    res = [block2d_dojo(mech, F, btol=btol, mode=mode, undercut=undercut) for F = Fsref]
    write_csv(["F", "x", "grad"], [[Fsref[i]; res[i]...] for i=1:length(Fsref)], "data/$(filename)_dojo_"*scn(btol)[2:end]*".csv")
    x = [r[1] for r in res]
    ∇x = [r[2] for r in res]
    plot!(plt[1,1], Fsref,  x, xlabel="F", ylabel="x", linewidth=3.0, linestyle=:dot,
        label="Dojo btol="*scn(btol))
    plot!(plt[2,1], Fsref, ∇x, xlabel="F", ylabel="∇x", linewidth=3.0, linestyle=:dot,
        label="Dojo btol="*scn(btol)[2:end])
    display(plt)
end

# ## gradient bundles
for Σ in [1e-2, 1e-3, 1e-4]
    res = [block2d_gradientbundle(mech, F, N=500, Σ=Σ*I, mode=mode) for F = Fs]
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

