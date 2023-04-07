################################################################################
# Block (2D)
################################################################################
struct Block2D end

function block2d(;
    representation=:minimal,
    timestep=0.05,
    gravity=-9.81,
    friction_coefficient=0.8,
    info=nothing,
    seed=1,
    contact=true,
    vis=Visualizer(),
    name=:robot,
    infeasible_control=false,
    opts_step=SolverOptions(rtol=3.0e-4, btol=3.0e-4, undercut=1.5),
    opts_grad=SolverOptions(rtol=3.0e-4, btol=3.0e-4, undercut=1.5),
    T=Float64)

    mechanism = get_mechanism(:block2d;
        timestep,
        gravity,
        friction_coefficient)

    initialize!(mechanism, :block2d)

    if representation == :minimal
        nx = minimal_dimension(mechanism)
    elseif representation == :maximal
        nx = maximal_dimension(mechanism)
    end
    nu_inf = input_dimension(mechanism)
    nu = infeasible_control ? nu_inf : nu_inf - 2 # remove last 2 controls
    no = nx

    aspace = BoxSpace(nu,
        low=(-1.0e8 * ones(nu)),
        high=(1.0e8 * ones(nu)))
    ospace = BoxSpace(no,
        low=(-Inf * ones(no)),
        high=(Inf * ones(no)))

    rng = MersenneTwister(seed)
    z = get_maximal_state(mechanism)
    x = representation == :minimal ? maximal_to_minimal(mechanism, z) : z
    fx = zeros(nx, nx)
    fu = zeros(nx, nu)

    u_prev = zeros(nu)
    control_mask = infeasible_control ? I(nu) : [I(nu) zeros(nu, 2)]
    control_scaling = Diagonal(ones(nu))

    build_robot(mechanism, vis=vis, name=name)

    TYPES = [Block2D, T, typeof(mechanism), typeof(aspace), typeof(ospace), typeof(info)]
    Environment{TYPES...}(mechanism, representation, aspace, ospace,
        x, fx, fu,
        u_prev, control_mask' * control_scaling,
        nx, nu, no,
        info,
        [rng], vis,
        opts_step, opts_grad)
end

# function initialize_block2d!(mechanism::Mechanism{T};
#     position=[0, 1],
#     velocity=[0, 0],
#     orientation=0,
#     angular_velocity=0) where T

#     if length(mechanism.contacts) > 0
#         model = mechanism.contacts[1].model
#         side = model.collision.contact_origin[2]
#         offset = model.collision.contact_radius
#         z = side + offset
#     else
#         z = 0.0
#     end

#     body = mechanism.bodies[1]

#     set_maximal_configurations!(body,
#         x=[0; position] + [0, 0.0 , z],
#         q=RotX(orientation))
#     set_maximal_velocities!(body,
#         v=[0; velocity],
#         ω=[angular_velocity, 0, 0])
# end
