"""
	dynamics(y, env, x, u, w)

	evaluates an environment's dynamics (in-place)

	y: state at next time step 
	env: Environment 
	x: state at current time step 
	u: input 
	w: system parameters
"""
function dynamics(y, env::Environment, x, u, w)
	step(env, x, u)[1]
    y .= env.state
end

"""
	dynamics_jacobian_state(dx, env, x, u, w)

	evaluates an environment's dynamics Jacobian wrt state (in-place)

	dx: state Jacobian at next time step 
	env: Environment 
	x: state at current time step 
	u: input 
	w: system parameters
"""
function dynamics_jacobian_state(dx, env::Environment, x, u, w)
	step(env, x, u, diff=true)
    dx .= env.dynamics_jacobian_state
end

"""
	dynamics_jacobian_input(du, env, x, u, w)

	evaluates an environment's dynamics Jacobian wrt input (in-place)

	du: input Jacobian at next time step 
	env: Environment 
	x: state at current time step 
	u: input 
	w: system parameters
"""
function dynamics_jacobian_input(du, env::Environment, x, u, w)
	# step(env, x, u, diff=true) # this is run in dynamics_jacobian_state
	du .= env.dynamics_jacobian_input
end