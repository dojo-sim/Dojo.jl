# ### Setup
# PKG_SETUP
using Dojo
using DojoEnvironments
using LinearAlgebra

# ### Environment
quadrotor_env = get_environment(:quadrotor_waypoint; horizon=200)

# ### Controllers
trans_z_mode = normalize([1;1;1;1])

function velocity_controller!(environment, v_des)
    v_is = get_state(environment)[9]
    error_z = v_des - v_is
    thrust = (10*error_z - 1*v_is + 5.1)*trans_z_mode # P, D, feedforward

    rpm = thrust*20
    set_input!(environment, rpm)
end

function position_controller!(environment, z_des)
    z_is = get_state(environment)[3]
    v_des = z_des - z_is
    velocity_controller!(environment, v_des)
end

function controller!(environment, k)
    position_controller!(environment, 0.3)
end

# ### Simulate
simulate!(quadrotor_env, controller!; record=true)

# ### Visualize
vis = visualize(quadrotor_env)
render(vis)
