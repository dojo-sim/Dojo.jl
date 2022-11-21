using Dojo
using DojoEnvironments
using IterativeLQR
using LinearAlgebra

env = get_environment(:panda, 
    representation=:minimal, 
    timestep=0.002,
    gravity=-0*9.81);

# env = DojoEnvironments.panda(; 
#     representation=:minimal, 
#     timestep=0.002,
#     gravity=-0*9.81);

joint1 = get_joint(env.mechanism,:joint1)
joint2 = get_joint(env.mechanism,:joint2)
joint3 = get_joint(env.mechanism,:joint3)
joint4 = get_joint(env.mechanism,:joint4)
joint5 = get_joint(env.mechanism,:joint5)
joint6 = get_joint(env.mechanism,:joint6)
joint7 = get_joint(env.mechanism,:joint7)
joint8 = get_joint(env.mechanism,:jointf1)
joint9 = get_joint(env.mechanism,:jointf2)

joints = [joint1;joint2;joint3;joint4;joint5;joint6;joint7;joint8;joint9]

function read_pos_and_vel()
    f1 = open("src/panda/methods/pos.txt", "r")
    f2 = open("src/panda/methods/vel.txt", "r")
    while !eof(f1)
        for i = 1:9
            s = readline(f1)
            if i != 1
                append!(qpos, -parse(Float64, s))
            else
                append!(qpos, parse(Float64, s))
            end
        end
    end
    close(f1)
    while !eof(f2)
        for i = 1:9
            s = readline(f2)
            if i != 1
                append!(qvel, -parse(Float64, s))
            else
                append!(qvel, parse(Float64, s))
            end
        end
    end
    close(f2)
end

function read_ctrl()
    ctrls = []
    f = open("src/panda/methods/ctrl.txt", "r")
    while !eof(f)
        s = readline(f)
        append!(ctrls, parse(Float64, s))
    end
    close(f)
end

function fill_storage()
    for i = first:last
        for j = 1:7
            set_minimal_coordinates!(env.mechanism, joints[j], [qpos[i*9+j]])
        end
        for j = 1:10
            storage.x[j][i-first+1] = env.mechanism.bodies[j].state.x1
            storage.q[j][i-first+1] = env.mechanism.bodies[j].state.q1
        end
    end
end

# After smooth
# step = []
# for i = 1:9
#     append!(step, (qpos[last*9+i]-qpos[first*9+i])/(last-first))
# end
# storage_smooth = Storage(500, 11)
# for i = first:last
#     set_minimal_coordinates!(env.mechanism, joint1, [qpos[first*9+1]+step[1]*(i-first)])
#     set_minimal_coordinates!(env.mechanism, joint2, [qpos[first*9+2]+step[2]*(i-first)])
#     set_minimal_coordinates!(env.mechanism, joint3, [qpos[first*9+3]+step[3]*(i-first)])
#     set_minimal_coordinates!(env.mechanism, joint4, [qpos[first*9+4]+step[4]*(i-first)])
#     set_minimal_coordinates!(env.mechanism, joint5, [qpos[first*9+5]+step[5]*(i-first)])
#     set_minimal_coordinates!(env.mechanism, joint6, [qpos[first*9+6]+step[6]*(i-first)])
#     set_minimal_coordinates!(env.mechanism, joint7, [qpos[first*9+7]+step[7]*(i-first)])
#     # set_minimal_coordinates!(env.mechanism, joint8, [qpos[first*9+9]+step[9]*(i-first)])
#     # set_minimal_coordinates!(env.mechanism, joint9, [qpos[first*9+8]+step[8]*(i-first)])
#     for j = 1:11
#         storage_smooth.x[j][i-first+1] = env.mechanism.bodies[j].state.x1
#         storage_smooth.q[j][i-first+1] = env.mechanism.bodies[j].state.q1
#     end
# end

function initialize()
    for j = 1:7
        set_minimal_coordinates!(env.mechanism, joints[j], [qpos[first*9+j]])
    end
end

function pid!(joint, pos, vel, err_sum)
    kp = 100
    kd = 10
    ki = 10
    p = (pos-minimal_coordinates(env.mechanism, joint)[1])*kp
    d = (vel-minimal_velocities(env.mechanism, joint)[1])*kd
    err_sum += pos-minimal_coordinates(env.mechanism, joint)[1]
    i = err_sum*ki
    ctrl = p + d + i
    return ctrl, err_sum
end

function controller!(mechanism, t)
    err_sum = zeros(9)
    for j = 1:9
        if j == 8 || j == 9
            write(f, string(0.0))
            write(f, '\n')
        else
            ctrl, err_sum[j] = pid!(joints[j], qpos[(t+first-1)*9+j], qvel[(t+first-1)*9+j], err_sum[j])
            set_input!(joints[j], [ctrl])
            ctrls[(t-1)*9+j] = ctrl
            write(f, string(ctrl))
            write(f, '\n')
        end
    end
end

# Read qpos and qvel
qpos = zeros(0)
qvel = zeros(0)
read_pos_and_vel()
first = 20
last = Int(length(qpos)/9-1)
steps = last-first+1

# Visualize the trajectory from mujoco
storage = Storage(steps, 11)
fill_storage()
visualize(env.mechanism, storage);

# Do PID control and write control inputs into file
ctrls = zeros(steps*9)
f = open("ctrl.txt", "w")
initialize()
storage_pid = simulate!(env.mechanism, steps*0.002, controller!, record = true)
visualize(env.mechanism, storage_pid);
close(f)

# Visualize the smooth trajectory
# visualize(env.mechanism, storage_smooth);

# Read contol inputs
read_ctrl()

# Do the optimization
n = env.num_states
m = env.num_inputs
T = steps

dyn = IterativeLQR.Dynamics(
    (y, x, u, w) -> dynamics(y, env, x, u, w), 
    (dx, x, u, w) -> dynamics_jacobian_state(dx, env, x, u, w),
    (du, x, u, w) -> dynamics_jacobian_input(du, env, x, u, w),
    n, n, m)
model = [dyn for t = 1:T-1]

# Initial state
x1 = zeros(18)
for i = 1:7
    x1[(i-1)*2+1] = qpos[first*9+i]
end

# Terminal state
xT = zeros(18)
for i = 1:7
    xT[(i-1)*2+1] = qpos[last*9+i]
end

# Ctrl input from pid controller
u_guess = [ctrls[(t-1)*9+1:(t-1)*9+9] for t = 1:T-1]

# Trajectory of initial guess
traj = IterativeLQR.rollout(model, x1, u_guess)
storage = generate_storage(env.mechanism, [env.representation == :minimal ? minimal_to_maximal(env.mechanism, x) : x for x in traj])
visualize(env.mechanism, storage, vis=env.vis, build=true)

# Objective function
t = [1 0.001 1 0.001 1 0.001 1 0.001 1 0.001 1 0.001 1 0.001 0.001 0.001 0.001 0.001]; # Currently dont care about fingers

function ot(x, u)
    return transpose(x-xT) * Diagonal(1.0 * t) * (x-xT) + transpose(u) * Diagonal(1.0 * ones(m)) * u
end

function oT(x, u)
    return transpose(x-xT) * Diagonal(1.0 * t) * (x-xT)
end
# ot = (x, u, w) -> transpose(x-xT) * Diagonal(1.0 * t) * (x-xT) + transpose(u) * Diagonal(1.0 * ones(m)) * u
# oT = (x, u, w) -> transpose(x-xT) * Diagonal(1.0 * t) * (x-xT)
ct = IterativeLQR.Cost(ot, n, m)
cT = IterativeLQR.Cost(oT, n, 0)
obj = [[ct for t = 1:T-1]..., cT]

# constraints
function initial_constraint(x, u)
    cons = x - x1
    for i = 1:9
        cons[(i-1)*2+2] = 0
    end
    return cons
end

function terminal_constraint(x, u)
    cons = x - xT
    for i = 1:9
        cons[(i-1)*2+2] = 0
    end
    return cons
end
con1 = IterativeLQR.Constraint(initial_constraint, n, 0)
cont = IterativeLQR.Constraint()
conT = IterativeLQR.Constraint(terminal_constraint, n, 0)
cons = [con1, [cont for t = 2: T-1]..., conT]

#Solver
s = IterativeLQR.Solver(model, obj, cons,
    options = IterativeLQR.Options(
        max_iterations = 30, 
        verbose = true))
IterativeLQR.initialize_controls!(s, u_guess)
IterativeLQR.initialize_states!(s, x1)

#Solve
# @time IterativeLQR.solve!(s)
IterativeLQR.solve!(s)

#Solution
traj_sol, u_sol = IterativeLQR.get_trajectory(s)
@show IterativeLQR.eval_obj(s.m_data.obj.costs, s.m_data.x, s.m_data.u, s.m_data.w)
@show s.s_data.iter[1]
@show norm(terminal_constraint(s.m_data.x[T], zeros(0)), Inf)

#Visualize
storage = generate_storage(env.mechanism, [env.representation == :minimal ? minimal_to_maximal(env.mechanism, x) : x for x in traj_sol])
visualize(env.mechanism, storage, vis=env.vis, build=true)