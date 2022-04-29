using DirectTrajectoryOptimization
const DTO = DirectTrajectoryOptimization

function centroidal_quadruped_dyn(model, env, h, y, x, u, w)

    # dimensions
    nq = model.nq
    nu = model.nu

    # configurations

    q1⁻ = x[1:nq]
    q2⁻ = x[nq .+ (1:nq)]
    q2⁺ = y[1:nq]
    q3⁺ = y[nq .+ (1:nq)]

    # control
    u_control = u[1:nu]
    γ = u[nu .+ (1:4)]
    β = u[nu + 4 .+ (1:16)]

    E = [1.0 0.0 -1.0 0.0;
         0.0 1.0 0.0 -1.0] # friction mapping
    J = J_func(model, env, q2⁺)
    λ = transpose(J[1:12, :]) * [
                        [E * β[0  .+ (1:4)]; γ[1]];
                        [E * β[4  .+ (1:4)]; γ[2]];
                        [E * β[8  .+ (1:4)]; γ[3]];
                        [E * β[12 .+ (1:4)]; γ[4]]
                       ]
    [
     q2⁺ - q2⁻;
     dynamics(model, h, q1⁻, q2⁺, u_control, zeros(model.nw), λ, q3⁺)
    ]
end

function centroidal_quadruped_dyn1(model, env, h, y, x, u, w)
    nx = 2 * model.nq
    [
     centroidal_quadruped_dyn(model, env, h, y, x, u, w);
     y[nx .+ (1:5)] - [u[model.nu .+ (1:4)]; u[end]];
     y[nx + 5 .+ (1:nx)] - x[1:nx];
     y[nx + 5 + nx .+ (1:model.nu)] - u[1:model.nu];
    ]
end

function centroidal_quadruped_dynt(model, env, h, y, x, u, w)
    nx = 2 * model.nq
    [
     centroidal_quadruped_dyn(model, env, h, y, x, u, w);
     y[nx .+ (1:5)] - [u[model.nu .+ (1:4)]; u[end]];
     y[nx + 5 .+ (1:nx)] - x[nx + 5 .+ (1:nx)];
     y[nx + 5 + nx .+ (1:model.nu)] - u[1:model.nu];
    ]
end

function contact_constraints_inequality_1(model, env, h, x, u, w)
    nq = model.nq
    nu = model.nu
    nx = 2nq

    q2 = x[1:nq]
    q3 = x[nq .+ (1:nq)]

    u_control = u[1:nu]
    γ = u[nu .+ (1:4)]
    β = u[nu + 4 .+ (1:16)]
    ψ = u[nu + 4 + 16 .+ (1:4)]
    η = u[nu + 4 + 16 + 4 .+ (1:16)]
    sα = u[nu + 4 + 16 + 4 + 16 .+ (1:1)]

    ϕ = ϕ_func(model, env, q3)[1:4]

    μ = model.μ_world
    fc = μ .* γ[1:4] - [sum(β[0 .+ (1:4)]); sum(β[4 .+ (1:4)]); sum(β[8 .+ (1:4)]); sum(β[12 .+ (1:4)]);]

    [
     -ϕ;
     -fc;
     β .* η .- sα;
     ψ .* fc  .- sα;
    ]
end

function contact_constraints_inequality_t(model, env, h, x, u, w)
    nq = model.nq
    nu = model.nu
    nx = 2nq

    q2 = x[1:nq]
    q3 = x[nq .+ (1:nq)]

    u_control = u[1:nu]
    γ = u[nu .+ (1:4)]
    β = u[nu + 4 .+ (1:16)]
    ψ = u[nu + 4 + 16 .+ (1:4)]
    η = u[nu + 4 + 16 + 4 .+ (1:16)]
    sα = u[nu + 4 + 16 + 4 + 16 .+ (1:1)]

    ϕ = ϕ_func(model, env, q3)[1:4]
    γ⁻ = x[nx .+ (1:4)]
    sα⁻ = x[nx + 4 .+ (1:1)]

    μ = model.μ_world
    fc = μ .* γ[1:4] - [sum(β[0 .+ (1:4)]); sum(β[4 .+ (1:4)]); sum(β[8 .+ (1:4)]); sum(β[12 .+ (1:4)]);]

    [
     -ϕ;
     -fc;
     γ⁻ .* ϕ .- sα⁻;
     β .* η .- sα;
     ψ .* fc  .- sα;
    ]
end


function contact_constraints_inequality_T(model, env, h, x, u, w)
    nq = model.nq
    nx = 2nq

    q2 = x[1:nq]
    q3 = x[nq .+ (1:nq)]

    ϕ = ϕ_func(model, env, q3)[1:4]
    γ⁻ = x[nx .+ (1:4)]
    sα⁻ = x[nx + 4 .+ (1:1)]

    [
     -ϕ;
     γ⁻ .* ϕ .- sα⁻;
    ]
end

function contact_constraints_equality(model, env, h, x, u, w)
    nq = model.nq
    nu = model.nu

    q2 = x[1:nq]
    q3 = x[nq .+ (1:nq)]

    γ = u[nu .+ (1:4)]
    β = u[nu + 4 .+ (1:16)]
    ψ = u[nu + 4 + 16 .+ (1:4)]
    η = u[nu + 4 + 16 + 4 .+ (1:16)]
    sα = u[nu + 4 + 16 + 4 + 16 .+ (1:1)]

    E = [1.0 0.0 -1.0 0.0;
         0.0 1.0 0.0 -1.0]
    v = (q3 - q2) ./ h[1]
    vT = vcat([E' * v[(i-1) * 3 .+ (1:2)] for i = 1:4]...)

    ψ_stack = vcat([ψi * ones(4) for ψi in ψ]...)

    [
     η - vT - ψ_stack;
    ]
end
