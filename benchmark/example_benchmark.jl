using Dojo

###############################
### Example benchmark for Atlas
###############################
mech = get_mechanism(:atlas, 
    timestep=0.1, 
    gravity=-9.81, 
    friction_coefficient=0.5, 
    damper=25.0, 
    spring=1.0, 
    contact_feet=true, 
    contact_body=true, 
    model_type=:simple)

Dojo.initialize!(mech, :atlas_stance,
    body_position=[0.0, 0.0, 1.0], 
    body_orientation=[0.0, 0.2, 0.1])

SUITE["atlas"] = @benchmarkable simulate!($mech, 2.5, opts=SolverOptions(rtol=1.0e-6, btol=1.0e-6)) samples=2



###############################
### Example benchmark for Pendulum
###############################
mechanism = get_mechanism(:pendulum,
    timestep=0.01,
    gravity=-9.81,
    damper=5.0,
    spring=0.0)

function controller!(mechanism, t)
    ## Target state
    x_goal = [1.0 * π; 0.0]

    ## Current state
    x = get_minimal_state(mechanism)

    ## Gains
    K = [5.0 0.5] * 0.1

    # Control inputs
    u = -K * (x - x_goal)
    set_input!(mechanism, u)
end

initialize!(mechanism, :pendulum,
    angle=0.0 * π,
    angular_velocity=0.0);

SUITE["pendulum"] = @benchmarkable simulate!($mechanism, 2.0, $controller!) samples=2