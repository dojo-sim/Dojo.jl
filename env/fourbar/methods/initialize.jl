function get_fourbar(; timestep::T=0.01, gravity=[0.0; 0.0; -9.81],
    model=:fourbar, spring=0.0, damper=1.0) where T

    path = joinpath(@__DIR__, "../deps/$(String(model)).urdf")
    mech = Mechanism(path, false, T, gravity=gravity, timestep=timestep, spring=spring, damper=damper)
    return mech
end

function initialize_fourbar!(mechanism::Mechanism; θ=0.0, ω1=0.0, ω2=0.0) where T
    zero_velocity!(mechanism)
    set_minimal_coordinates_velocities!(mechanism, get_joint_constraint(mechanism, :jointb1); xmin=[ -θ, ω2])
    set_minimal_coordinates_velocities!(mechanism, get_joint_constraint(mechanism, :joint12); xmin=[+2θ, 0])
    set_minimal_coordinates_velocities!(mechanism, get_joint_constraint(mechanism, :jointb3); xmin=[ +θ, ω1])
    set_minimal_coordinates_velocities!(mechanism, get_joint_constraint(mechanism, :joint34); xmin=[-2θ, 0])
end


# vis = Visualizer()
# open(vis)
