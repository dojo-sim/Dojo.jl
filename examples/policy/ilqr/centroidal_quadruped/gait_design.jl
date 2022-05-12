
function trotting_configuration(model::CentroidalQuadruped114, Tm::Int; timestep=0.01, velocity=0.15,
        foot_height=0.08, mode::Symbol=:left)

    nq = model.nq
    nx = 2nq
    body_height = RD.centroidal_quadruped.body_height
    foot_x = RD.centroidal_quadruped.foot_x
    foot_y = RD.centroidal_quadruped.foot_y
    stride = velocity * timestep * (Tm - 1)


    orientations = [zeros(3) for i=1:Tm]
    body_positions = [[0;0;body_height] for i=1:Tm]

    flight_positions = [zeros(3) for i=1:Tm]
    for (i,θ) in enumerate(range(0, stop=π, length=Tm))
        x = -stride/2 * cos(θ)
        z = foot_height * sin(θ)
        flight_positions[i] = [x; 0; z]
    end
    stance_positions = [[s;0;0] for s in range(stride/2, stop=-stride/2, length=Tm)]

    foot_1_positions = [[+foot_x; +foot_y; 0.0] for i=1:Tm]
    foot_2_positions = [[+foot_x; -foot_y; 0.0] for i=1:Tm]
    foot_3_positions = [[-foot_x; +foot_y; 0.0] for i=1:Tm]
    foot_4_positions = [[-foot_x; -foot_y; 0.0] for i=1:Tm]

    if mode == :left
        foot_1_positions .+= flight_positions
        foot_2_positions .+= stance_positions
        foot_3_positions .+= stance_positions
        foot_4_positions .+= flight_positions
    else
        foot_1_positions .+= stance_positions
        foot_2_positions .+= flight_positions
        foot_3_positions .+= flight_positions
        foot_4_positions .+= stance_positions
    end

    configurations = [[
        body_positions[i];
        orientations[i];
        foot_1_positions[i];
        foot_2_positions[i];
        foot_3_positions[i];
        foot_4_positions[i];
        ] for i=1:Tm]
    return configurations
end

function trotting_gait(model::CentroidalQuadruped114, Tm::Int; timestep=0.01, velocity=0.15,
        foot_height=0.08)

    T = 2Tm - 1
    nq = model.nq

    configurations_left = trotting_configuration(centroidal_quadruped, Tm,
        timestep=timestep, velocity=velocity, foot_height=foot_height, mode=:left)
    configurations_right = trotting_configuration(centroidal_quadruped, Tm,
        timestep=timestep, velocity=velocity, foot_height=foot_height, mode=:right)

    configurations = [configurations_left[1:end-1]; configurations_right]
    # add forward velocity
    for (i,q) in enumerate(configurations)
        q[1] += (i-1) * velocity * timestep
        q[7] += (i-1) * velocity * timestep
        q[10] += (i-1) * velocity * timestep
        q[13] += (i-1) * velocity * timestep
        q[16] += (i-1) * velocity * timestep
    end
    # compute finite difference velocity
    states = [[q; zeros(nq)] for q in configurations]
    for i = 2:T
        states[i][nq .+ (1:nq)] = (configurations[i] - configurations[i-1]) / timestep
    end
    states[1][nq .+ (1:nq)] = states[end][nq .+ (1:nq)]
    return states
end
