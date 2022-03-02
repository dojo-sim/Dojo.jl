set_camera!(env.vis, 
	zoom=5.0, 
	cam_pos=[0.0, -5.0, 0.0])

z = [minimal_to_maximal(env.mechanism, x) for x in x_sol]
t = 1 #10, 20, 30, 41
set_robot(env.vis, env.mechanism, z[t])