################################################################################
# Step
################################################################################
function dynamics(y, env::Environment, x, u, w)
	step(env, x, u)[1]
    y .= env.state
end

function dynamics_jacobian_state(dx, env::Environment, x, u, w)
	step(env, x, u, diff=true)
    dx .= env.dynamics_jacobian_state
end

function dynamics_jacobian_input(du, env::Environment, x, u, w)
	# step(env, x, u, diff=true) # this is run in dynamics_jacobian_state
	du .= env.dynamics_jacobian_input
end