function num_var(model::CentroidalQuadruped114)
    nq = model.nq
    nc = model.nc
    nb = 2

    nq + nc + nc + nc + nc * nb + nc + nc * nb
end

function indices_z(model::CentroidalQuadruped114)
    nq = model.nq
    nc = model.nc
    nb = 2
    q = collect(1:nq)
    γ = collect(nq .+ (1:nc))
    sγ = collect(nq + nc .+ (1:nc))
    ψ = collect(nq + nc + nc .+ (1:nc))
    b = collect(nq + nc + nc + nc .+ (1:(nc * nb)))
    sψ = collect(nq + nc + nc + nc + nc * nb .+ (1:nc))
    sb = collect(nq + nc + nc + nc + nc * nb + nc .+ (1:(nc * nb)))
    IndicesZ(q, γ, sγ, ψ, b, sψ, sb)
end

function RD.indices_optimization(model::CentroidalQuadruped114)
    nq = model.nq
    nc = model.nc
    nb = 2
    nz = num_var(model)

    IndicesOptimization(
        nz,
        nz,
        [collect(nq .+ (1:nc)), collect(nq + nc .+ (1:nc))],
        [collect(nq .+ (1:nc)), collect(nq + nc .+ (1:nc))],
        [
            [collect([nq + 2nc + 1, (nq + 2nc + nc + 0 .+ (1:2))...]), collect([nq + 2nc + nc + nc * nb + 1, (nq + 2nc + nc + nc * nb + nc + 0 .+ (1:2))...])],
            [collect([nq + 2nc + 2, (nq + 2nc + nc + 2 .+ (1:2))...]), collect([nq + 2nc + nc + nc * nb + 2, (nq + 2nc + nc + nc * nb + nc + 2 .+ (1:2))...])],
            [collect([nq + 2nc + 3, (nq + 2nc + nc + 4 .+ (1:2))...]), collect([nq + 2nc + nc + nc * nb + 3, (nq + 2nc + nc + nc * nb + nc + 4 .+ (1:2))...])],
            [collect([nq + 2nc + 4, (nq + 2nc + nc + 6 .+ (1:2))...]), collect([nq + 2nc + nc + nc * nb + 4, (nq + 2nc + nc + nc * nb + nc + 6 .+ (1:2))...])],
            # [collect([nq + 2nc + 5, (nq + 2nc + nc + 8 .+ (1:2))...]), collect([nq + 2nc + nc + nc * nb + 5, (nq + 2nc + nc + nc * nb + nc + 8 .+ (1:2))...])],
        ],
        [
            [collect([nq + 2nc + 1, (nq + 2nc + nc + 0 .+ (1:2))...]), collect([nq + 2nc + nc + nc * nb + 1, (nq + 2nc + nc + nc * nb + nc + 0 .+ (1:2))...])],
            [collect([nq + 2nc + 2, (nq + 2nc + nc + 2 .+ (1:2))...]), collect([nq + 2nc + nc + nc * nb + 2, (nq + 2nc + nc + nc * nb + nc + 2 .+ (1:2))...])],
            [collect([nq + 2nc + 3, (nq + 2nc + nc + 4 .+ (1:2))...]), collect([nq + 2nc + nc + nc * nb + 3, (nq + 2nc + nc + nc * nb + nc + 4 .+ (1:2))...])],
            [collect([nq + 2nc + 4, (nq + 2nc + nc + 6 .+ (1:2))...]), collect([nq + 2nc + nc + nc * nb + 4, (nq + 2nc + nc + nc * nb + nc + 6 .+ (1:2))...])],
            # [collect([nq + 2nc + 5, (nq + 2nc + nc + 8 .+ (1:2))...]), collect([nq + 2nc + nc + nc * nb + 5, (nq + 2nc + nc + nc * nb + nc + 8 .+ (1:2))...])],
        ],
        collect(1:(nq + nc + nc * nb + nc)),
        collect(nq + nc + nc * nb + nc .+ (1:nc)),
        collect(nq + nc + nc * nb + nc + nc .+ (1:(nc * 3))),
        [collect(nq + nc + nc * nb + nc + nc + (i - 1) * 3 .+ (1:3)) for i = 1:nc],
        collect(nq + nc + nc * nb + nc .+ (1:(nc + nc * 3))))
end

Trajectory(model::CentroidalQuadruped114, T) = Trajectory(model, T, nc=model.nc, nb=(2 * model.nc))
GradientTrajectory(model::CentroidalQuadruped114, T) = GradientTrajectory(model, T, nc=model.nc, nb=(2 * model.nc))

function residual(model::CentroidalQuadruped114, mass_matrix, dynamics_bias, kinematics, kinematics_jacobians, z, θ, μ)

    # dimensions
    nq = model.nq
    nu = model.nu
    nw = model.nw
    nc = model.nc
    nb = 2
    nf = length(friction_coefficients(model))

    # unpack data
    q0 = θ[1:nq]
    q1 = θ[nq .+ (1:nq)]
    u1 = θ[2nq .+ (1:nu)]
    w1 = θ[2nq + nu .+ (1:nw)]
    friction = θ[2nq + nu + nw .+ (1:nf)]
    h = θ[2nq + nu + nw + nf .+ (1:1)]

    # unpack variables
    q2 = z[1:nq]
    γ1 = z[nq .+ (1:nc)]
    sγ1 = z[nq + nc .+ (1:nc)]
    ψ1 = z[nq + nc + nc .+ (1:nc)]
    b1 = z[nq + nc + nc + nc .+ (1:(nc * nb))]
    sψ1 = z[nq + nc + nc + nc + (nc * nb) .+ (1:nc)]
    sb1 = z[nq + nc + nc + nc + (nc * nb) + nc .+ (1:(nc * nb))]

    # rotation matrix at surface
    rotation = [I(3) for i = 1:nc] #TODO non-flat surfaces

    # tangential velocity
    v1 = (q2 - q1) / h[1]
    vT = [(rotation[i] * kinematics_jacobians[i](q2) * v1)[1:2] for i = 1:nc]

    # contact forces
    J = contact_jacobian(model, q1)
    λ1 = transpose(J) * vcat([transpose(rotation[i]) * [b1[(i - 1) * nb .+ (1:nb)]; γ1[i]] for i = 1:nc]...)

    # residual
    [
    dynamics(model, mass_matrix, dynamics_bias, h, q0, q1, u1, w1, λ1, q2);
    signed_distance(model, q2) - sγ1;
    vcat([vT[i] - sb1[(i - 1) * nb .+ (1:nb)] for i = 1:nc]...);
    vcat([friction[i] * γ1[i] - ψ1[i] for i = 1:nc]...);
    γ1 .* sγ1 .- μ[1];
    vcat([cone_product([ψ1[i]; b1[(i - 1) * nb .+ (1:nb)]], [sψ1[i]; sb1[(i - 1) * nb .+ (1:nb)]]) - [μ[1]; 0.0; 0.0] for i = 1:nc]...);
    ]
end

function codegen_dynamics(model::CentroidalQuadruped114)
    mm = q -> mass_matrix(model, q)
    db = (q, q̇) -> dynamics_bias(model, q, q̇)
    return mm, db
end

function dynamics(model::CentroidalQuadruped114, mass_matrix, dynamics_bias, h, q0, q1, u1, w1, λ1, q2)
    # evalutate at midpoint
    qm1 = 0.5 * (q0 + q1)
    vm1 = (q1 - q0) / h[1]
    qm2 = 0.5 * (q1 + q2)
    vm2 = (q2 - q1) / h[1]

    D1L1, D2L1 = lagrangian_derivatives(mass_matrix, dynamics_bias, qm1, vm1)
    D1L2, D2L2 = lagrangian_derivatives(mass_matrix, dynamics_bias, qm2, vm2)

    d = 0.5 * h[1] * D1L1 + D2L1 + 0.5 * h[1] * D1L2 - D2L2 # variational integrator (midpoint)
    d .+= transpose(input_jacobian(model, qm2)) * u1             # control inputs
    d .+= λ1                                                # contact impulses
    @show "here"
    # d[(floating_base_dim(model) + 1):end] .-= model.friction_joint .* [vm2[(floating_base_dim(model) + 1):end] - [vm2[1:3]] # joint friction
    d[1:3] .-= model.friction_joint .* (vm2[1:3] - vm2[7:9]) # joint friction
    d[1:3] .-= model.friction_joint .* (vm2[1:3] - vm2[10:12]) # joint friction
    d[1:3] .-= model.friction_joint .* (vm2[1:3] - vm2[13:15]) # joint friction
    d[1:3] .-= model.friction_joint .* (vm2[1:3] - vm2[16:18]) # joint friction

    d[7:9]   .-= model.friction_joint .* (vm2[7:9] - vm2[1:3]) # joint friction
    d[10:12] .-= model.friction_joint .* (vm2[10:12] - vm2[1:3]) # joint friction
    d[13:15] .-= model.friction_joint .* (vm2[13:15] - vm2[1:3]) # joint friction
    d[16:18] .-= model.friction_joint .* (vm2[16:18] - vm2[1:3]) # joint friction
    return d
end
