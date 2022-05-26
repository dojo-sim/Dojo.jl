function get_fourbar(; 
    timestep=0.01, 
    gravity=[0.0; 0.0; -9.81],
    model=:fourbar, 
    spring=0.0, 
    damper=0.0,
    parse_damper=true,
    T=Float64)

    path = joinpath(@__DIR__, "../deps/$(String(model)).urdf")
    mech = Mechanism(path; floating=false, T,
        gravity, 
        timestep, 
        parse_damper)

    # Adding springs and dampers
    set_springs!(mech.joints, spring)
    set_dampers!(mech.joints, damper)

    return mech
end

function initialize_fourbar!(mechanism::Mechanism; 
    angle=0.0, 
    angular_velocity=szeros(2)) where T

    zero_velocity!(mechanism)
    set_minimal_coordinates_velocities!(mechanism, get_joint(mechanism, :jointb1); 
        xmin=[ -angle, angular_velocity[1]])
    set_minimal_coordinates_velocities!(mechanism, get_joint(mechanism, :joint12); 
        xmin=[+2angle, 0])
    set_minimal_coordinates_velocities!(mechanism, get_joint(mechanism, :jointb3); 
        xmin=[ +angle, angular_velocity[2]])
    set_minimal_coordinates_velocities!(mechanism, get_joint(mechanism, :joint34); 
        xmin=[-2angle, 0])
end