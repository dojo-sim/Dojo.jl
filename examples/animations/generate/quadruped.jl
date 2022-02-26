x_shift = deepcopy(x_sol)
for x in x_shift
    x[3] += 0.01
end
z = [minimal_to_maximal(env.mechanism, x) for x in x_shift]

t = 1 #10, 20, 30, 41
set_robot(env.vis, env.mechanism, z[41])

# ## visualize
x_view = [[x_shift[1] for t = 1:15]..., x_shift..., [x_shift[end] for t = 1:15]...]
visualize(env, x_view)

set_camera!(vis, cam_pos=[0,-50,0], zoom=30)
set_floor!(env.vis, z=0.01)
