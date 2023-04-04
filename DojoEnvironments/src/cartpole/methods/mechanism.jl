function get_cartpole(; 
    timestep=0.1, 
    input_scaling=timestep, 
    gravity=-9.81, 
    slider_mass=1.0,
    pendulum_mass=1.0,
    pendulum_length=1.0,
    color=RGBA(0.7, 0.7, 0.7, 1.0)
    springs=0.0, 
    dampers=0.0,
    limits=false,
    joint_limits=Dict(),
    T=Float64)

    # mechanism
    origin = Origin{Float64}()
    slider = Capsule(1.5 * 0.075, 1.0, slider_mass, 
        orientation_offset=RotX(0.5 * π), color)
    pendulum = Capsule(radius, pendulum_length, pendulum_mass, color)
    bodies = [slider, pendulum]
    
    joint_origin_slider = JointConstraint(Prismatic(origin, slider, Y_AXIS))
    joint_slider_pendulum = JointConstraint(Revolute(slider, pendulum, X_AXIS; 
        child_vertex=-0.5*pendulum_length*Z_AXIS))
    joints = [joint_origin_slider, joint_slider_pendulum]

    mech = Mechanism(origin, bodies, joints;
        gravity, timestep, input_scaling)

    # springs and dampers
    set_springs!(mech.joints, springs)
    set_dampers!(mech.joints, dampers)

    # joint limits    
    if limits
        joints = set_limits(mech, joint_limits)

        mech = Mechanism(Origin{T}(), mech.bodies, joints;
            gravity, timestep, input_scaling)
    end

    # construction finished
    return mech
end

function initialize_cartpole!(mech::Mechanism{T,Nn,Ne,Nb}; 
    mode=:up, 
    pendulum_length=1.0) where {T,Nn,Ne,Nb}

# origin to slider
set_maximal_configurations!(mech.origin, mech.bodies[1])
set_maximal_velocities!(mech.bodies[1], 
    v=[0.0; 0.0; 0.0],
    ω=zeros(3))

# slider to pendulum
if mode == :down
    set_maximal_configurations!(mech.bodies[1], mech.bodies[2], 
        Δx=[0.0; 0.0; -0.5 * pendulum_length], 
        Δq=RotX(π))
    set_maximal_velocities!(mech.bodies[2], 
        v=zeros(3), 
        ω=zeros(3))
elseif mode == :up
    set_maximal_configurations!(mech.bodies[1], mech.bodies[2], 
        Δx=[0.0; 0.0; 0.5 * pendulum_length], 
        Δq=RotX(0))
    set_maximal_velocities!(mech.bodies[2], 
        v=zeros(3), 
        ω=zeros(3))
end
end

function mujoco_inertia!(mech)
mech.bodies[1].m = 1.0063
mech.bodies[1].J = Diagonal([0.106974, 0.106974, 0.00636812])

mech.bodies[2].m = 0.4321
mech.bodies[2].J = Diagonal([0.0422274, 0.0422274, 0.0012155])
end