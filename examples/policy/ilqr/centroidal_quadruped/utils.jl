function reset!(dynamics_model::Simulator; r_tol=1e-5, κ_tol = 1e-3, undercut = 5.0)
    dynamics_model.ip.opts.r_tol = r_tol
    dynamics_model.ip.opts.κ_tol = κ_tol
    dynamics_model.ip.opts.undercut = undercut
    return nothing
end

function trajectory_generator(solver::Solver, dynamics_model::Simulator, gait::Vector;
        u_guess=[[4,0,16, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0.0] for i=1:length(gait)-1],
        initial_disturbance = zeros(36),
        )

    # ## initialization
    ############################################################################
    parameters = deepcopy(gait)
    T = length(parameters)
    x1 = deepcopy(parameters[1])
    xT = deepcopy(parameters[end])

    # ## initialize
    dyn = solver.problem.model.dynamics
    x_guess = iLQR.rollout(dyn, x1 + initial_disturbance, u_guess, parameters)

    iLQR.initialize_controls!(solver, u_guess)
    iLQR.initialize_states!(solver, x_guess)
    reset!(dynamics_model)
    local_continuation_callback!(solver::Solver) = continuation_callback!(solver, dynamics_model)

    # ## solve
    @time iLQR.constrained_ilqr_solve!(solver, augmented_lagrangian_callback! = local_continuation_callback!)

    # ## solution
    x_sol, u_sol = iLQR.get_trajectory(solver)
    K_sol = solver.policy.K
    return x_sol, u_sol, K_sol, solver.data.max_violation[1]
end

function solver_generator(dynamics_model::Simulator, gait::Vector;
        u_hover=[4,0,16, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0.0],
        )

    # ## initialization
    ############################################################################
    nq = dynamics_model.model.nq
    nx = 2nq
    nu = dynamics_model.model.nu
    nw = dynamics_model.model.nw
    nu_infeasible = 6
    parameters = deepcopy(gait)
    T = length(parameters)
    x1 = deepcopy(parameters[1])
    xT = deepcopy(parameters[end])

    # ## dynamics
    ############################################################################
    function f1(y, x, u, w)
        @views u_ctrl = u[1:nu]
        @views x_di = x[1:nx]
        RD.dynamics(dynamics_model, view(y, 1:nx), x_di, u_ctrl, w)
        return nothing
    end

    function f1x(dx, x, u, w)
        @views u_ctrl = u[1:nu]
        @views x_di = x[1:nx]
        dx .= 0.0
        RD.dynamics_jacobian_state(dynamics_model, view(dx, 1:nx, 1:nx), x_di, u_ctrl, w)
        return nothing
    end

    function f1u(du, x, u, w)
        @views u_ctrl = u[1:nu]
        @views x_di = x[1:nx]
        du .= 0.0
        RD.dynamics_jacobian_input(dynamics_model, view(du, 1:nx, 1:nu), x_di, u_ctrl, w)
        return nothing
    end

    function ft(y, x, u, w)
        @views u_ctrl = u[1:nu]
        @views x_di = x[1:nx]
        RD.dynamics(dynamics_model, view(y, 1:nx), x_di, u_ctrl, w)
        return nothing
    end

    function ftx(dx, x, u, w)
        @views u_ctrl = u[1:nu]
        @views x_di = x[1:nx]
        dx .= 0.0
        RD.dynamics_jacobian_state(dynamics_model, view(dx, 1:nx, 1:nx), x_di, u_ctrl, w)
        return nothing
    end

    function ftu(du, x, u, w)
        @views u_ctrl = u[1:nu]
        @views x_di = x[1:nx]
        RD.dynamics_jacobian_input(dynamics_model, view(du, 1:nx, 1:nu), x_di, u_ctrl, w)
        return nothing
    end

    # user-provided dynamics and gradients
    dyn1 = iLQR.Dynamics(f1, f1x, f1u, nx, nx, nu)
    dynt = iLQR.Dynamics(ft, ftx, ftu, nx, nx, nu)

    dyn = [dyn1, [dynt for t = 2:T-1]...]

    # ## objective
    ############################################################################
    function ot(x, u, w)
        J = 0.0
        qbody = [1e-0, 1e-0, 1e+2]
        qfoot = [1e-0, 1e-0, 1e+2]
        q = 1e-0 * [1e-0*qbody; 1e+0*ones(3); [qfoot; qfoot; qfoot; qfoot]]
        v = 3e+1 * [1e-1*ones(3); 1e-0*ones(3); 1e-1*ones(12)]
        r = 1e-1 * [ones(6); [1e-1,1,1e-3]; [1e-1,1,1e-3]; [1e-1,1,1e-3]; [1e-1,1,1e-3]]
        ex = x[1:nx] - w
        eu = u[1:nu] - u_hover
        J += 0.5 * transpose(ex) * Diagonal([q; v]) * ex
        J += 0.5 * transpose(eu) * Diagonal(r) * eu
        return J
    end

    function oT(x, u, w)
        J = 0.0
        qbody = [1e-0, 1e-0, 1e+1]
        qfoot = [1e-0, 1e-0, 1e+2]
        q = 1e-0 * [1e-0*qbody; 1e-1*ones(3); [qfoot; qfoot; qfoot; qfoot]]
        v = 3e+1 * [1e-1*ones(3); 1e-0*ones(3); 1e-1*ones(12)]
        ex = x[1:nx] - w
        J += 0.5 * transpose(ex) * Diagonal([q; v]) * ex
        return J
    end

    ct = iLQR.Cost(ot, nx, nu, num_parameter=nx)
    cT = iLQR.Cost(oT, nx, 0, num_parameter=nx)
    obj = [[ct for t = 1:(T - 1)]..., cT]


    # ## constraints
    ############################################################################
    ul = -1.0 * [1e-3*ones(nu_infeasible); 1e3ones(nu-nu_infeasible)]
    uu = +1.0 * [1e-3*ones(nu_infeasible); 1e3ones(nu-nu_infeasible)]

    function con1(x, u, w)
        [
            1e-2 * (ul - u[1:nu]);
            1e-2 * (u[1:nu] - uu);
        ]
    end

    function cont(x, u, w)
        [
            1e-2 * (ul - u[1:nu]);
            1e-2 * (u[1:nu] - uu);
        ]
    end

    function goal(x, u, w)
        [
            1e-2 * (x[1:nq+3] - xT[1:nq+3]);
        ]
    end

    con_policy1 = iLQR.Constraint(con1, nx, nu, num_parameter=nx, indices_inequality=collect(1:2nu))
    con_policyt = iLQR.Constraint(cont, nx, nu, num_parameter=nx, indices_inequality=collect(1:2nu))
    con_policyT = iLQR.Constraint(goal, nx, 0)

    cons = [con_policy1, [con_policyt for t = 2:T-1]..., con_policyT]


    # ## problem
    ############################################################################
    opts = iLQR.Options(line_search=:armijo,
        max_iterations=75,
        max_dual_updates=30,
        objective_tolerance=1e-3,
        lagrangian_gradient_tolerance=1e-3,
        constraint_tolerance=1e-4,
        initial_constraint_penalty=1e-3,
        scaling_penalty=3.0,
        max_penalty=1e4,
        verbose=true)

    p = iLQR.Solver(dyn, obj, cons, options=opts, parameters=parameters)
end

function initial_disturbance_sampler(; configuration_amplitude=0.03, velocity_amplitude=0.01)
    nq = 18
    configuration = configuration_amplitude * (rand(nq) .- 0.5) .* [1,1,1, 1,1,1, 1,1,1, 1,1,1, 1,1,1, 1,1,1]
    velocity = velocity_amplitude * (rand(nq) .- 0.5) .* [1,1,1, 1,1,1, 1,1,1, 1,1,1, 1,1,1, 1,1,1]
    initial_disturbance = [configuration; velocity]
    return project_state(initial_disturbance)
end

function project_state(x)
    nq = 18
    configuration = x[1:nq]
    velocity = x[nq .+ (1:nq)]
    configuration[[3,9,12,15,18]] = max.(0, configuration[[3,9,12,15,18]])
    xp = [configuration; velocity]
    return xp
end


function trajectory_sampler(solver::Solver, dynamics_model::Simulator, gait::Vector, sample_number::Int;
        u_guess=[[4,0,16, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0.0] for i=1:length(gait)-1],
        configuration_amplitude=0.03,
        velocity_amplitude=0.01,
        )

    for i = 1:sample_number
        initial_disturbance = initial_disturbance_sampler(;
            configuration_amplitude=configuration_amplitude,
            velocity_amplitude=velocity_amplitude)

        x_sol, u_sol, K_sol, max_violation = trajectory_generator(solver, dynamics_model, gait,
            u_guess=u_guess,
            initial_disturbance=initial_disturbance,
            )

        if max_violation < 1e-3
            JLD2.jldsave(joinpath(@__DIR__, "dataset/centroidal_quadruped_$i.jld2"),
                x_sol=x_sol, u_sol=u_sol, K_sol=K_sol)
        end
    end

end
