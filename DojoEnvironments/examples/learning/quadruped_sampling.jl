using Dojo
using DojoEnvironments
using Random
using LinearAlgebra

rng = MersenneTwister(1)
N = 2500
M = 50
paramcontainer = [[0.1; 0; 1; 0; -1.5]]
paramstorage = [[0.1; 0; 1; 0; -1.5]]
bias = zeros(5)
distance = 0.0
explore_factor = 0.1
distancestorage = zeros(M)

env = get_environment(:quadruped_sampling; horizon=N, timestep=0.002, limits=false, gravity=-9.81, contact_body=false)


legmovement(k,a,b,c,offset) = a*cos(k*b*0.01*2*pi+offset)+c
Kp = [100;80;60]
Kd = [5;4;3]

function controller!(x, k)
    angle21 = legmovement(k,paramcontainer[1][2],paramcontainer[1][1],paramcontainer[1][3],0)
    angle22 = legmovement(k,paramcontainer[1][2],paramcontainer[1][1],paramcontainer[1][3],pi)
    angle31 = legmovement(k,paramcontainer[1][4],paramcontainer[1][1],paramcontainer[1][5],-pi/2)
    angle32 = legmovement(k,paramcontainer[1][4],paramcontainer[1][1],paramcontainer[1][5],pi/2)

    u = zeros(12)

    for i=1:4
         θ1 = x[12+(i-1)*6+1]
        dθ1 = x[12+(i-1)*6+2]
         θ2 = x[12+(i-1)*6+3]
        dθ2 = x[12+(i-1)*6+4]
         θ3 = x[12+(i-1)*6+5]
        dθ3 = x[12+(i-1)*6+6]

        if i == 1 || i == 4
            u[(i-1)*3+1] = Kp[1]*(0-θ1) + Kd[1]*(0-dθ1)
            u[(i-1)*3+2] = Kp[2]*(angle21-θ2) + Kd[2]*(0-dθ2)
            u[(i-1)*3+3] = Kp[3]*(angle31-θ3) + Kd[3]*(0-dθ3)
        else
            u[(i-1)*3+1] = Kp[1]*(0-θ1) + Kd[1]*(0-dθ1)
            u[(i-1)*3+2] = Kp[2]*(angle22-θ2) + Kd[2]*(0-dθ2)
            u[(i-1)*3+3] = Kp[3]*(angle32-θ3) + Kd[3]*(0-dθ3)
        end
    end

    return u
end

function reset_state!(env)
    initialize!(env.mechanism, :quadruped; body_position=[0;0;-0.43], hip_angle=0, thigh_angle=paramcontainer[1][3], calf_angle=paramcontainer[1][5]) 

    calf_state = get_body(env.mechanism, :FR_calf).state
    position = get_sdf(get_contact(env.mechanism, :FR_calf_contact), Dojo.current_position(calf_state), Dojo.current_orientation(calf_state))

    initialize!(env.mechanism, :quadruped; body_position=[0;0;-position-0.43], hip_angle=0, thigh_angle=paramcontainer[1][3], calf_angle=paramcontainer[1][5]) 
end

function rollout(env; record=false)
    for k=1:N
        x = DojoEnvironments.get_state(env)
        if x[3] < 0 || !all(isfinite.(x))
            display("  failed")
            break
        end
        u = controller!(x, k)
        DojoEnvironments.step!(env, x, u; k, record)
    end
end

for i=1:M
    display("run: $i")

    if bias == zeros(5)
        paramcontainer[1] += randn!(rng, zeros(5))*explore_factor
    else
        paramcontainer[1] += randn!(rng, zeros(5))*0.002 + normalize(bias)*0.01
    end

    reset_state!(env)
    x0 = DojoEnvironments.get_state(env)[1]

    try
        rollout(env)
    catch
        display("  errored")
        paramcontainer[1] = paramstorage[end]
        bias = zeros(5)
        explore_factor *= 0.9
    end

    new_state = DojoEnvironments.get_state(env)
    distancenew = new_state[1] - x0

    if !all(isfinite.(new_state)) || distancenew <= distance
        display("  unsuccessful")
        !all(isfinite.(new_state)) && display("  nans")
        paramcontainer[1] = paramstorage[end]
        bias = zeros(5)
        explore_factor *= 0.9
    else
        display("  successful")
        distance = distancenew
        push!(paramstorage,paramcontainer[1])
        bias = paramstorage[end]-paramstorage[end-1]
        explore_factor = 0.1
    end

    display("  distance: $distancenew")
    distancestorage[i] = distance
end

reset_state!(env)
rollout(env; record=true)

visualize(env.mechanism, env.storage)