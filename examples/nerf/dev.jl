include(joinpath(module_dir(), "examples/nerf/utils.jl"))

abstract type NerfContact{T,N} <: Contact{T,N}
end
data_dim(model::NerfContact) = 6 # [offset, p]

abstract type NerfObject{T}
end

mutable struct NerfObject133{T} <: NerfObject{T}
    render_kwargs_test::Dict
    iso::T
    scale::T
end

function sampler(p::AbstractVector, ϵ)
    ps = ϵ * [
        +0 +0 +0;
        +0 +0 +1;
        +0 +1 +0;
        +1 +0 +0;
        +0 +0 -1;
        +0 -1 +0;
        -1 +0 +0;
        +1 +1 +1;
        +1 +1 -1;
        +1 -1 +1;
        +1 -1 -1;
        -1 +1 +1;
        -1 +1 -1;
        -1 -1 +1;
        -1 -1 -1;]
    for i in 1:size(ps)[1]
        ps[i,:] .+= p
    end
    return ps
end

function nerf_density(nerf::NerfObject, p::AbstractVector)
    # p in Nerf frame
    # p = Matrix(p')
    ps = sampler(p, FDEPS)
    ps = convert(Matrix{Float32}, ps)
    d = density_query(nerf.render_kwargs_test, ps) ./ nerf.scale
    d = mean(d)
    return convert(Float64, d)
end

function nerf_density_gradient(nerf::NerfObject, p::AbstractVector)
    # p and ∇ in Nerf frame
    p = convert(Vector{Float32}, p)
    p = Matrix(p')
    ∇ = density_gradient_query(nerf.render_kwargs_test, p) ./ nerf.scale
    return convert(Vector{Float64}, ∇[1,:])
end

mutable struct NerfContact133{T,N} <: NerfContact{T,N}
    ainv3::Adjoint{T,SVector{3,T}} # inverse matrix
    contact_point::SVector{3,T}
    offset::SVector{3,T} # offset of the contact point expressed in the world frame
    nerf::NerfObject{T}

    function NerfContact133(body::Body{T}, normal::AbstractVector, nerf::NerfObject;
            p = szeros(T, 3), offset::AbstractVector = szeros(T, 3)) where T
        V1, V2, V3 = orthogonalcols(normal) # gives two plane vectors and the original normal axis
        A = [V1 V2 V3]
        Ainv = inv(A)
        ainv3 = Ainv[3,SA[1; 2; 3]]'
        new{Float64,2}(ainv3, p, offset, nerf)
    end
end


function constraint(mechanism, contact::ContactConstraint{T,N,Nc,Cs}) where {T,N,Nc,Cs<:NerfContact{T,N}}
    model = contact.model
    body = get_body(mechanism, contact.parent_id)
    x3, q3 = next_configuration(body.state, mechanism.timestep)
    p = inv(q3) * (model.contact_point - x3) # position of the origin point in the frame of centered on P expressed in the body's frame
    c = model.nerf.iso - nerf_density(model.nerf, p)  # > 0
    SVector{1,T}(c - contact.impulses_dual[2][1])
end

function constraint_jacobian_velocity(model::NerfContact, x3::AbstractVector, q3::UnitQuaternion,
    x2::AbstractVector, v25::AbstractVector, q2::UnitQuaternion, ϕ25::AbstractVector, λ, timestep)
    # V = model.ainv3 * timestep
    # Ω = model.ainv3 * ∂vrotate∂q(model.contact_point, q3) * rotational_integrator_jacobian_velocity(q2, ϕ25, timestep)
    # return [V Ω]
    p = inv(q3) * (model.contact_point - x3)
    ∂nerf∂p = FiniteDiff.finite_difference_jacobian(p -> nerf_density(model.nerf,p), p, absstep=FDEPS, relstep=FDEPS)
    # ∂nerf∂p ./= norm(∂nerf∂p) + 1e-2
    X = -∂nerf∂p * -rotation_matrix(inv(q3))
    Q = -∂nerf∂p * ∂qrotation_matrix_inv(q3, model.contact_point - x3)
    return [X Q] * integrator_jacobian_velocity(q2, ϕ25, timestep)
    # inv(q3) * (model.contact_point - x3)
    # FiniteDiff.finite_difference_jacobian(
    #     vϕ -> -nerf_density(model.nerf,
    #         inv(next_orientation(q2, vϕ[SUnitRange(4,6)], timestep)) *
    #         (model.contact_point - next_position(x2, vϕ[SUnitRange(1,3)], timestep))),
    #     [v25; ϕ25], absstep=FDEPS, relstep=FDEPS)
end

function constraint_jacobian_configuration(model::NerfContact, x3::AbstractVector, q3::UnitQuaternion,
    x2::AbstractVector, v25::AbstractVector, q2::UnitQuaternion, ϕ25::AbstractVector, λ, timestep)
    # X = model.ainv3
    # Q = model.ainv3 * ∂vrotate∂q(model.contact_point, q3)
    # return [X Q]
    p = inv(q3) * (model.contact_point - x3)
    ∂nerf∂p = FiniteDiff.finite_difference_jacobian(p -> nerf_density(model.nerf,p), p, absstep=FDEPS, relstep=FDEPS)
    # ∂nerf∂p ./= norm(∂nerf∂p) + 1e-2
    X = -∂nerf∂p * -rotation_matrix(inv(q))
    Q = -∂nerf∂p * ∂qrotation_matrix_inv(q, model.contact_point - x) * LVᵀmat(q)
    return [X Q]
    # FiniteDiff.finite_difference_jacobian(
    #     xq -> -nerf_density(model.nerf, inv(UnitQuaternion(xq[4:7]..., false)) * (model.contact_point - xq[1:3])),
    #     [x3; vector(q3)], absstep=FDEPS, relstep=FDEPS)

end

function impulse_map(model::NerfContact, x::AbstractVector, q::UnitQuaternion, λ)
    # X = model.ainv3
    # # q * ... is a rotation by quaternion q it is equivalent to Vmat() * Lmat(q) * Rmat(q)' * Vᵀmat() * ...
    # Q = - X * q * skew(model.contact_point - vrotate(model.offset, inv(q)))
    p = inv(q) * (model.contact_point - x)
    ∂nerf∂p = FiniteDiff.finite_difference_jacobian(p -> nerf_density(model.nerf,p), p, absstep=FDEPS, relstep=FDEPS)
    ∂nerf∂p ./= norm(∂nerf∂p) + 1e-2
    X = -∂nerf∂p * -rotation_matrix(inv(q))
    Q = -∂nerf∂p * ∂qrotation_matrix_inv(q, model.contact_point - x) * LVᵀmat(q)
    # Q = FiniteDiff.finite_difference_jacobian(
    #     q -> -nerf_density(model.nerf, -model.contact_point +inv(UnitQuaternion(q..., false))*(model.offset-x)),
    #     vector(q), absstep=FDEPS, relstep=FDEPS) * LVᵀmat(q)
    return transpose([X Q])
end

function force_mapping(model::NerfContact, x::AbstractVector, q::UnitQuaternion)
    # X = FiniteDiff.finite_difference_jacobian(
    #     x -> -nerf_density(model.nerf, -model.contact_point + inv(q) * (model.offset-x)),
    #     x, absstep=FDEPS, relstep=FDEPS)
    p = inv(q) * (model.contact_point - x)
    ∂nerf∂p = FiniteDiff.finite_difference_jacobian(p -> nerf_density(model.nerf,p), p, absstep=FDEPS, relstep=FDEPS)
    ∂nerf∂p ./= norm(∂nerf∂p) + 1e-2
    X = -∂nerf∂p * -rotation_matrix(inv(q))
    return X
end

function set_matrix_vector_entries!(mechanism::Mechanism, matrix_entry::Entry, vector_entry::Entry,
    contact::ContactConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs<:NerfContact{T,N},N½}
    # ∇impulses[dual .* impulses - μ; g - s] = [diag(dual); -diag(0,1,1)]
    # ∇dual[dual .* impulses - μ; g - s] = [diag(impulses); -diag(1,0,0)]
    γ = contact.impulses[2]
    s = contact.impulses_dual[2]

    ∇s = hcat(γ, -Diagonal(sones(N½)))
    ∇γ = hcat(s, -Diagonal(szeros(N½)))
    matrix_entry.value = hcat(∇s, ∇γ)

    # [-dual .* impulses + μ; -g + s]
    vector_entry.value = vcat(-complementarityμ(mechanism, contact), -constraint(mechanism, contact))
    return
end

function get_nerf(; timestep::T=0.01, gravity=[0.0; 0.0; -9.81], cf::T=0.8, radius=0.05, mass=1.0,
        iso=65.0, scale=100.0, contact_points=[szeros(3)], contact::Bool=true,
        contact_type::Symbol=:impact) where T
    origin = Origin{T}(name=:origin)
    mass = mass
    bodies = [Sphere(radius, mass, name=:sphere)]
    joints = [JointConstraint(Floating(origin, bodies[1]), name = :floating_joint)]
    mechanism = Mechanism(origin, bodies, joints, timestep=timestep, gravity=gravity)

    if contact
        contact = [0,0,0.0]
        normal = [0,0,1.0]
        scaled_iso = iso / scale
        nerf = NerfObject133(generate_test_nerf(), scaled_iso, scale)
        body = get_body(mechanism, :sphere)
        models = [NerfContact133(body, normal, nerf; p=p) for p in contact_points]
        contacts = [ContactConstraint((model, body.id, nothing); name=Symbol(:nerf_contact, i)) for (i,model) in enumerate(models)]
        set_position!(mechanism, get_joint_constraint(mechanism, :floating_joint), [0;0;2radius;zeros(3)])
        mechanism = Mechanism(origin, bodies, joints, contacts, gravity=gravity, timestep=timestep)
    end
    return mechanism
end

function initialize_nerf!(mechanism::Mechanism; x::AbstractVector{T}=zeros(3),
        q::UnitQuaternion{T}=one(UnitQuaternion), v::AbstractVector{T}=zeros(3),
        ω::AbstractVector{T}=zeros(3)) where T
    r = collect(mechanism.bodies)[1].shape.r
    joint = get_joint_constraint(mechanism, :floating_joint)
    zero_velocity!(mechanism)
    set_position!(mechanism, joint, [x+[0,0,r] rotation_vector(q)])
    set_velocity!(mechanism, joint, [v; ω])
end

function feasibility_linesearch!(α, mechanism, contact::ContactConstraint{T,N,Nc,Cs,N½},
        vector_entry::Entry, τort, τsoc; scaling::Bool = false) where {T,N,Nc,Cs<:Union{NerfContact{T,N},ImpactContact{T,N},LinearContact{T,N}},N½}
    s = contact.impulses_dual[2]
    γ = contact.impulses[2]
    Δs = vector_entry.value[1:N½]
    Δγ = vector_entry.value[N½ .+ (1:N½)]

    αs_ort = positive_orthant_step_length(s, Δs, τ = τort)
    αγ_ort = positive_orthant_step_length(γ, Δγ, τ = τort)

    return min(α, αs_ort, αγ_ort)
end

function visualize_contact_points(mechanism::Mechanism, vis::Visualizer)
    mat = MeshPhongMaterial(color= RGBA(0.05, 0.05, 0.05, 1.0))
    zero_velocity!(mechanism)
    for (i,contact) in enumerate(mechanism.contacts)
        model = contact.model
        contact_point = model.contact_point
        obj = MeshCat.HyperEllipsoid(Point(contact_point...), 0.008*Vec(1,1,1.))
        setobject!(vis[:contact_points][Symbol(i)], obj, mat)
    end
    return nothing
end







# p = rand(Float32, 100,3) ./ 5
# render_kwargs_test = generate_test_nerf()
# density = density_query(render_kwargs_test, p)
# density_grad = density_gradient_query(render_kwargs_test, p)
#
#
# origin = Origin{Float64}(name=:origin)
# mass = 1.0
# timestep = 0.01
# gravity = 0.9
# radius = 0.10
# bodies = [Sphere(radius, mass, name=:sphere)]
# joints = [JointConstraint(Floating(origin, bodies[1]), name=:floating_joint)]
# mechanism = Mechanism(origin, bodies, joints, timestep=timestep, gravity=gravity)
#
#
# contact = [0,0,0.0]
# normal = [0,0,1.0]
# nerf = NerfObject133(render_kwargs_test, 65.0/100, 100.0)
# body = get_body(mechanism, :sphere)
#
# offset = [0,0,0.0]
# model = NerfContact133(body, normal, nerf; p=szeros(Float64, 3), offset=offset)
# contact = ContactConstraint((model, body.id, nothing); name=Symbol(:nerf_contact, 1))
# set_position!(mechanism, get_joint_constraint(mechanism, :floating_joint), [0;0;2radius;zeros(3)])
# mechanism = Mechanism(origin, bodies, joints, [contact], gravity=gravity, timestep=timestep)
# constraint(mechanism, contact)
# x2 = srand(3)*0.01
# q2 = UnitQuaternion(rand(4)...)
# x3 = srand(3)*0.01
# q3 = UnitQuaternion(rand(4)...)
# v25 = srand(3)*0.01
# ϕ25 = srand(3)*0.01
#
# constraint(mechanism, contact)
# constraint_jacobian_velocity(contact.model, x3, q3, x2, v25, q2, ϕ25, 1, timestep)
# constraint_jacobian_configuration(contact.model, x3, q3, x2, v25, q2, ϕ25, 1, timestep)
# impulse_map(contact.model, x3, q3, 0)
# force_mapping(contact.model, x3, q3)
#
#
#





# using MeshCat
# vis = Visualizer()
# open(vis)
# mesh = jldopen(joinpath(@__DIR__, "deps/bunny_mesh.jld"))["mesh"]
FDEPS = 1e-2

contact_points = vcat([[
    SVector(0.10*j,0.10*i,0.15+0.02*i^2+0.02*j^2) for i in -1:1] for j in -1:1]...)
radius = 0.35
mass = 1.0
iso = 65.0
scale = 100.0

mech = get_nerf(timestep=0.05, gravity=-0.5, radius=radius, mass=mass,
    iso=iso, scale=scale, contact_points=contact_points)
# initialize!(mech, :nerf, x=[0.00,-0.30,0.21], v=[0.00, 0.00, 0.00])
initialize!(mech, :nerf, x=[0.1,-0.0,0.27],
    q=axis_angle_to_quaternion([0,0π/2,0]), v=[0.00, 0.00, 0.00])

storage = simulate!(mech, 6.0, record=true, opts=SolverOptions(rtol=1e-3, btol=1e-3, verbose=true))
visualize(mech, storage, vis=vis)
MeshCat.setobject!(vis[:robot][:bodies]["body:1"][:mesh], mesh,
    MeshPhongMaterial(color=RGBA{Float32}(1, 0, 0, 0.95)))
visualize_contact_points(mech, vis)

storage2 = crop_s

convert_frames_to_video_and_gif("bunny_hammock")
convert_frames_to_video_and_gif("bunny_drop_point_sphere")


render_static(vis)

open(joinpath(@__DIR__, "nerf_hammock.html"), "w") do file
    write(file, static_html(vis))
end
