@setup_workload begin
    mechanisms = [
        :ant,
        :atlas,
        :block,
        :block2d,
        :cartpole,
        :dzhanibekov,
        :exoskeleton,
        :fourbar,
        :halfcheetah,
        :hopper,
        :humanoid,
        :npendulum,
        :nslider,
        :panda,
        :pendulum,
        :quadrotor,
        :quadruped,
        :raiberthopper,
        :slider,
        :snake,
        :sphere,
        :tippetop,
        :twister,
        :uuv,
        :walker,
        :youbot,
    ]

    environments = [
        :ant_ars,
        :cartpole_dqn,
        :pendulum,
        :quadruped_waypoint,
        :quadruped_sampling,
        :quadrotor_waypoint,
        :uuv_waypoint,
        :youbot_waypoint,
    ]

    @compile_workload begin
        # Simulate all mechanisms
        for name in mechanisms 
            mech = get_mechanism(name) 
            initialize!(mech, name)
            simulate!(mech, mech.timestep * 2)
        end

        # Simulate all environments
        for name in environments 
            env = get_environment(name; horizon=2)
            x = get_state(env)
            u = DojoEnvironments.input_map(env, nothing)
            set_input!(env, nothing)
            step!(env, x)
            step!(env, x, u)
            simulate!(env)
        end
    end
end