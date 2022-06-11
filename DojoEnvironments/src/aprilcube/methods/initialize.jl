function get_aprilcube(;
    timestep=0.01,
    gravity=[0.0; 0.0; -9.81],
    side=1.0,
    friction_coefficient=0.8,
    color=RGBA(0.3, 0.3, 0.3, 1.0),
    options=ColliderOptions(),
    T=Float64)

    body = Box(side, side, side, 1.0, color=color, name=:aprilcube)
    origin = Origin{T}(name=:origin)
    bodies = [body]

    joints = [JointConstraint(Floating(origin, bodies[1]), name=:floating_joint)]
    mechanism = Mechanism(origin, bodies, joints,
        timestep=timestep,
        gravity=gravity)

    contact = [0.0, 0.0, 0.0]
    normal = [0.0, 0.0, 1.0]
    center_of_mass = szeros(3)
    collider = SoftCollider(body.mass, body.inertia, center_of_mass,
        x -> box_densities(x, b=side/2*ones(3)),
        x -> box_gradient(x, b=side/2*ones(3)),
        N=5000, particle_noise=0.03)
    contacts = [soft_contact_constraint(get_body(mechanism, :aprilcube), normal, collider,
        friction_coefficient=friction_coefficient,
        collider_origin=-collider.center_of_mass,
        name=:contact_1)]
    set_minimal_coordinates!(mechanism, get_joint(mechanism, :floating_joint), [0.0; 0.0; 0.5; zeros(3)])
    mechanism = Mechanism(origin, bodies, joints, contacts,
        gravity=gravity,
        timestep=timestep)
    mechanism.contacts[1].model.collision.collider.options = options
    return mechanism
end

function initialize_aprilcube!(mechanism::Mechanism{T};
    position=zeros(3),
    orientation=one(Quaternion),
    velocity=zeros(3),
    angular_velocity=zeros(3)) where T

    r = 0.50
    joint = get_joint(mechanism, :floating_joint)
    zero_velocity!(mechanism)
    set_minimal_coordinates!(mechanism, joint, [position + [0.0, 0.0, r] rotation_vector(orientation)])
    set_minimal_velocities!(mechanism, joint, [velocity; angular_velocity])
end

################################################################################
# Utils
################################################################################
function box_sdf(p::AbstractVector; b::AbstractVector=ones(3))
    q = abs.(p) - b
    return norm(max.(q, 0)) + min(maximum(q), 0);
end

function density_map(x)
    d = 2.0 - 120x/(0.03*4)
    return clamp(d, 0.0, 120.0)
end

function box_densities(ps::AbstractMatrix; b::AbstractVector=ones(3))
    N = size(ps)[1]
    d = zeros(N)
    for i = 1:N
        d[i] = density_map(box_sdf(ps[i,:], b=b))
    end
    return d
end

function box_gradient(p; b::AbstractVector=ones(3), δ=0.01)
    FiniteDiff.finite_difference_jacobian(p -> box_sdf(p, b=b), p, absstep=δ, relstep=0.01)
end

# plot(-0.1:0.001:0.1, density_map.(-0.1:0.001:0.1))

# collider = SoftCollider(1.0, rand(3,3), rand(3),
#    x -> box_densities(x, b=0.5/2*ones(3)),
#    x -> box_gradient(x, b=0.5/2*ones(3)),
#    N=5000, particle_noise=0.03)
# histogram(collider.densities)


# xrange = -0.5:0.01:0.5
# xrange = -0.0:0.01:0.0
# yrange = -0.5:0.01:0.5
# N = length(xrange)
# z = 0.1
# particles = slice_particles(xrange, yrange, z)
# nerf_object = get_nerf_object()
# density = density_query(nerf_object, particles)
# plot(density[290:360])
# density = reshape(density, (N,N))
#
# using Plots
# plt = heatmap(xrange, yrange, density,
#     c=cgrad([:blue, :white,:red, :yellow]),
#     xlabel="x", ylabel="y",
#     title="bunny")
#
#
#
# p = [0.9,0,0.0]
# b = [1,1,1.0]
# sdf_box(p, b)
# # density 0->100 in 6cm
