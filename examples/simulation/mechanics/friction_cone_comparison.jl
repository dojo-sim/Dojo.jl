# using Pkg
# Pkg.activate(joinpath(@__DIR__, ".."))
# Pkg.instantiate()

# ## Setup
using Dojo
using DojoEnvironments

# ## Parameters
friction_coefficient = 0.25
x0 = [-1.5, -0.50, 0.25]
q0 = one(Quaternion)
v0 = [4, 0.80, 0.0]
ω0 = zeros(3)

# ## Create blocks for comparison
mech_linear = get_mechanism(:block; friction_coefficient, contact_type=:linear, color=orange);
mech_nonlinear = get_mechanism(:block; friction_coefficient, contact_type=:nonlinear, color=cyan);

# ## Simulate block with linear friction 
initialize!(mech_linear, :block; position=x0, velocity=v0, orientation=q0, angular_velocity=ω0)
storage_linear = simulate!(mech_linear, 4.0; record=true)

# ## Simulate block with nonlinear friction 
initialize!(mech_nonlinear, :block; position=x0, velocity=v0, orientation=q0, angular_velocity=ω0)
storage_nonlinear = simulate!(mech_nonlinear, 4.0; record=true)

# ## Visualize
vis = Visualizer()
set_camera!(vis, cam_pos=[-0.01, 0.0, 100.0], zoom=30)
set_light!(vis, ambient=0.65, direction="Positive")

vis, animation = visualize(mech_linear, storage_linear; vis, name=:linear, return_animation=true);
vis = visualize(mech_nonlinear, storage_nonlinear; vis, name=:nonlinear, animation)
