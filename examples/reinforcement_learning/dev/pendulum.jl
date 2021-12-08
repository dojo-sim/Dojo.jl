################################################################################
# Script
################################################################################
vis = Visualizer()
env = make("Pendulum", 
    mode=:max, 
    # mode=:min,
    vis=vis)

open(vis)

for i = 1:1
    observation = reset(env)
    for t = 1:200
        render(env)
        sleep(0.001)
        u = sample(env.aspace)
        # u = [4.0]
        observation, reward, done, info = step(env, u)
        if done
            println("Episode finished after $(t+1) timesteps")
            break
        end
    end
end
close(env)
