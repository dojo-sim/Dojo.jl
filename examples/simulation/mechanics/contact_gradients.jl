# ### Setup
# PKG_SETUP_2
using Dojo 
using DojoEnvironments
using Plots
using LinearAlgebra

# ### Define gradient functions
function get_gradients(mechanism::Mechanism, F, mode; 
    rtol=1e-10, btol=1e-10, 
    undercut=1.0, no_progress_undercut=1.0)

    (mode == :friction) && (idx = 1)
    (mode == :impact) && (idx = 2)
    opts = SolverOptions(;
        rtol, btol, undercut, no_progress_undercut)

    initialize!(mechanism, :block2d, 
        position=[0, 0], velocity=[0, 0], 
        orientation=0, angular_velocity=0)

    x0 = get_minimal_state(mechanism)
    u0 = zeros(3)
    u0[idx] += F 

    _, ∇u = get_minimal_gradients!(mechanism, x0, u0; opts)

    Dojo.step_minimal_coordinates!(mechanism, x0, u0; opts)

    z1 = Dojo.get_next_state(mechanism)
    x1 = maximal_to_minimal(mechanism, z1)

    ∂x∂F = ∇u[idx,idx]
    Δx = x1[idx] - x0[idx]

    return Δx, ∂x∂F
end

# ### Gravity = Friction = 1 for nice plots
mech = get_mechanism(:block2d; 
    timestep=0.1, gravity=-1.0, 
    friction_coefficient=1.0, contact=true)

# ### Force range
F_ref = 0.5:0.01:1.5
tolerances = [1e-4, 1e-6, 1e-8, 1e-10]

# ### Gradient plots
gradient_plot = plot(layout=(2,2), size=(800,600), legend=:topleft);

# ### Impact plots
for btol in tolerances
    gradients = [get_gradients(mech, F, :impact; btol) for F = F_ref]
    x = [r[1] for r in gradients]
    ∇x = [r[2] for r in gradients] 
    plot!(gradient_plot[1,1], F_ref,  x, xlabel="Fz", ylabel="z", linewidth=3.0, linestyle=:dot,
        label="Dojo btol ="*Dojo.scn(btol))
    plot!(gradient_plot[2,1], F_ref, ∇x, xlabel="Fz", ylabel="∇z", linewidth=3.0, linestyle=:dot,
        label="Dojo btol ="*Dojo.scn(btol)[2:end])
end

# ### Friction plots
for btol in tolerances
    gradients = [get_gradients(mech, F, :friction; btol) for F = F_ref]
    x = [r[1] for r in gradients]
    ∇x = [r[2] for r in gradients]
    plot!(gradient_plot[1,2], F_ref,  x, xlabel="Fx", ylabel="x", linewidth=3.0, linestyle=:dot,
        label="Dojo btol ="*Dojo.scn(btol))
    plot!(gradient_plot[2,2], F_ref, ∇x, xlabel="Fx", ylabel="∇x", linewidth=3.0, linestyle=:dot,
        label="Dojo btol ="*Dojo.scn(btol)[2:end])
end

# ### Show plots
display(gradient_plot) # left: Impact, right: Friction