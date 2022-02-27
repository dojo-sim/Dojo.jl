abstract type NerfContact{T,N} <: Contact{T,N}
end

abstract type NerfObject{T}
end

mutable struct NerfObject111{T} <: NerfObject{T}
    radius::T
    std::T
    threshold::T
end

function density(nerf::NerfObject{T}, p::AbstractVector{T}) where T
    # p in Nerf frame
    std = nerf.std
    p0 = zeros(3)
    Δ = norm(p) - nerf.radius
    d = 1/(nerf.std * sqrt(2π)) * exp(-0.5*(Δ/std)^2)
    d = exp(-0.5*(Δ/std)^2)
    return d
end

function density_gradient(nerf::NerfObject{T}, p::AbstractVector{T}) where T
    # p and ∇ in Nerf frame
    ∇ = FiniteDiff.finite_difference_gradient(p -> density(nerf, p), p)
    return ∇
end

mutable struct NerfContact111{T,N} <: NerfContact{T,N}
    ainv3::Adjoint{T,SVector{3,T}} # inverse matrix
    p::SVector{3,T}
    offset::SVector{3,T} # offset of the contact point expressed in the world frame
    nerf::NerfObject{T}

    function NerfContact111(body::Body{T}, normal::AbstractVector, nerf::NerfObject;
            p = szeros(T, 3), offset::AbstractVector = szeros(T, 3)) where T
        V1, V2, V3 = orthogonal_columns(normal) # gives two plane vectors and the original normal axis
        A = [V1 V2 V3]
        Ainv = inv(A)
        ainv3 = Ainv[3,SA[1; 2; 3]]'
        new{Float64,2}(ainv3, p, offset, nerf)
    end
end

function constraint(mechanism, contact::ContactConstraint{T,N,Nc,Cs}) where {T,N,Nc,Cs<:Tuple{NerfContact{T,N}}}
    model = contact.model
    body = get_body(mechanism, contact.parent_id)
    x3, q3 = next_configuration(body.state, mechanism.timestep)
    # SVector{1,T}(model.ainv3 * (x3 + vector_rotate(model.contact_point,q3) - model.offset) - contact.impulses_dual[2][1])
    r = - model.contact_point + inv(q3) * (model.offset - x3) # position of the origin point in the frame of centered on P expressed in the body's frame
    c = model.nerf.threshold - density(model.nerf, r) # > 0
    SVector{1,T}(c - contact.impulses_dual[2][1])
end

function constraint_jacobian_velocity(model::NerfContact, x3::AbstractVector, q3::UnitQuaternion,
    x2::AbstractVector, v25::AbstractVector, q2::UnitQuaternion, ϕ25::AbstractVector, λ, timestep)
    # V = model.ainv3 * timestep
    # Ω = model.ainv3 * ∂vector_rotate∂q(model.contact_point, q3) * rotational_integrator_jacobian_velocity(q2, ϕ25, timestep)
    # return [V Ω]
    FiniteDiff.finite_difference_jacobian(
        vϕ -> -density(model.nerf,
            -model.contact_point +inv(next_orientation(q2, vϕ[SUnitRange(4,6)], timestep)) *
            (model.offset - next_position(x2, vϕ[SUnitRange(1,3)], timestep))),
        [v25; ϕ25])
end

function constraint_jacobian_configuration(model::NerfContact, x3::AbstractVector, q3::UnitQuaternion,
    x2::AbstractVector, v25::AbstractVector, q2::UnitQuaternion, ϕ25::AbstractVector, λ, timestep)
    # X = model.ainv3
    # Q = model.ainv3 * ∂vector_rotate∂q(model.contact_point, q3)
    # return [X Q]
    FiniteDiff.finite_difference_jacobian(
        xq -> -density(model.nerf, -model.contact_point +inv(UnitQuaternion(xq[4:7]..., false)) * (model.offset - xq[1:3])),
        [x3; vector(q3)])
end

function impulse_map(model::NerfContact, x::AbstractVector, q::UnitQuaternion, λ)
    # X = model.ainv3
    # # q * ... is a rotation by quaternion q it is equivalent to Vmat() * Lmat(q) * Rmat(q)' * Vᵀmat() * ...
    # Q = - X * q * skew(model.contact_point - vector_rotate(model.offset, inv(q)))
    X = FiniteDiff.finite_difference_jacobian(
        x -> -density(model.nerf, -model.contact_point +inv(q)*(model.offset-x)),
        x)
    Q = FiniteDiff.finite_difference_jacobian(
        q -> -density(model.nerf, -model.contact_point +inv(UnitQuaternion(q..., false))*(model.offset-x)),
        vector(q)) * LVᵀmat(q)
    return transpose([X Q])
end

function force_mapping(model::NerfContact, x::AbstractVector, q::UnitQuaternion)
    X = FiniteDiff.finite_difference_jacobian(
        x -> -density(model.nerf, -model.contact_point -inv(q)*(model.offset-x)),
        x)
    # X = model.ainv3
    return X
end

function set_matrix_vector_entries!(mechanism::Mechanism, matrix_entry::Entry, vector_entry::Entry,
    contact::ContactConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs<:Tuple{NerfContact{T,N}},N½}
    # ∇primal[dual .* primal - μ; g - s] = [diag(dual); -diag(0,1,1)]
    # ∇dual[dual .* primal - μ; g - s] = [diag(primal); -diag(1,0,0)]
    γ = contact.impulses[2]
    s = contact.impulses_dual[2]

    ∇s = hcat(γ, -Diagonal(sones(N½)))
    ∇γ = hcat(s, -Diagonal(szeros(N½)))
    matrix_entry.value = hcat(∇s, ∇γ)

    # [-dual .* primal + μ; -g + s]
    vector_entry.value = vcat(-complementarityμ(mechanism, contact), -constraint(mechanism, contact))
    return
end

function get_nerf(; timestep::T=0.01, gravity=[0.0; 0.0; -9.81], cf::T=0.8, radius=0.1,
        std=0.2, threshold=0.95, offsets=[szeros(3)], contact::Bool=true,
        contact_type::Symbol=:impact) where T
    origin = Origin{T}(name=:origin)
    mass = 1.0
    bodies = [Sphere(radius, mass, name=:sphere)]
    joints = [JointConstraint(Floating(origin, bodies[1]), name = :floating_joint)]
    mechanism = Mechanism(origin, bodies, joints, timestep=timestep, gravity=gravity)

    if contact
        contact = [0,0,0.0]
        normal = [0,0,1.0]
        nerf = NerfObject111(radius, std, threshold)
        body = get_body(mechanism, :sphere)
        models = [NerfContact111(body, normal, nerf; p=szeros(T, 3), offset=offset) for offset in offsets]
        contacts = [ContactConstraint((model, body.id, nothing); name=Symbol(:nerf_contact, i)) for (i,model) in enumerate(models)]
        set_position!(mechanism, get_joint_constraint(mechanism, :floating_joint), [0;0;2radius;zeros(3)])
        mechanism = Mechanism(origin, bodies, joints, contacts, gravity=gravity, timestep=timestep)
    end
    return mechanism
end

function initialize_nerf!(mechanism::Mechanism; x::AbstractVector{T}=zeros(3),
        q::UnitQuaternion{T}=one(UnitQuaternion), v::AbstractVector{T}=zeros(3),
        ω::AbstractVector{T}=zeros(3)) where T
    r = mechanism.bodies[1].shape.r
    joint = get_joint_constraint(mechanism, :floating_joint)
    zero_velocity!(mechanism)
    set_position!(mechanism, joint, [x+[0,0,r] rotation_vector(q)])
    set_minimal_velocities!(mechanism, joint, [v; ω])
end

function cone_line_search!(α, mechanism, contact::ContactConstraint{T,N,Nc,Cs,N½},
        vector_entry::Entry, τort, τsoc; scaling::Bool = false) where {T,N,Nc,Cs<:Tuple{Union{NerfContact{T,N},ImpactContact{T,N},LinearContact{T,N}}},N½}
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
    for (i,contact) in enumerate(mechanism.contacts)
        model = contact.model
        offset = model.offset
        obj = MeshCat.HyperEllipsoid(Point(offset...), 0.008*Vec(1,1,1.))
        setobject!(vis[:contact_points][Symbol(i)], obj, mat)
    end
    return nothing
end

# vis=visualizer()
# open(vis)
offsets = vcat([[SVector(0.10*j,0.10*i,0.0+0.02*i^2+0.02*j^2) for i in -0:0] for j in -0:0]...)
radius = 0.10
std = 0.03
threshold = 0.95
Δ = sqrt(-2 * std^2 * log(threshold))
radius_int = radius - Δ
radius_ext = radius + Δ

mech = get_nerf(timestep=0.001, gravity=-10.0, radius=radius,
    std=std, threshold=threshold, offsets=offsets)
# initialize!(mech, :nerf, x=[0.1,0.1,0.35], v=[0.03, 0.02, 0.1])
initialize!(mech, :nerf, x=[0.00,0.01,0.35], v=[0.00, 0.02, 0.1])
storage = simulate!(mech, 1.0, record=true, verbose=false)
visualize(mech, storage, vis=vis)
visualize_contact_points(mech, vis)

obj_int = MeshCat.HyperEllipsoid(Point(0.0, 0.0, 0.0), radius_int*Vec(1,1,1.))
obj_ext = MeshCat.HyperEllipsoid(Point(0.0, 0.0, 0.0), radius_ext*Vec(1,1,1.))
mat_int = MeshPhongMaterial(color= RGBA(0.25, 0.25, 0.25, 1.0))
mat_ext = MeshPhongMaterial(color= RGBA(0.05, 0.05, 0.05, 0.2))
setobject!(vis[:robot][:bodies]["body:1"][:ext], obj_ext, mat_ext)
setobject!(vis[:robot][:bodies]["body:1"][:int], obj_int, mat_int)



sdTorus( vec3 p, vec2 t )
{
  vec2 q = vec2(length(p.xz)-t.x,p.y);
  return length(q)-t.y;
}


nerf0 = NerfObject111(0.1, 0.2, 0.95)
p1 = [0.1, 0.2, -0.0]
density(nerf0, p1)
plot(-1:0.01:1, [density(nerf0, [0,0,i]) for i in -1:0.01:1])
convert_frames_to_video_and_gif("nerf_point_fast")



render_static(vis)

open(joinpath(@__DIR__, "nerf_hammock.html"), "w") do file
    write(file, static_html(vis))
end

a = 10
a = 10
# body0 =
# contact0 = NerfContact(body::Body{T}, normal::AbstractVector, nerf::NerfObject;
#         p = szeros(T, 3), offset::AbstractVector = szeros(T, 3))



#
# mech = get_humanoid()
# model = mech.contacts[1].constraints[1]
# q3 = UnitQuaternion(rand(4)...)
# X = model.ainv3
# Q0 = model.ainv3 * ∂vector_rotate∂q(model.contact_point, q3) * LVᵀmat(q3)
# Q1 = - X * q3 * skew(model.contact_point - vector_rotate(model.offset, inv(q3)))
#
# Q0 ./ Q1



# function closest_to_plane(nerf::NerfObject{T}, v::AbstractVector{T};
#         c0::AbstractVector{T}=v./norm(v)) where T
#     # find point c on the level set of the nerf that is the furthest in direction v
#     # this is also the closest point to the plane orthogonal to v.
#     # c0 is the initial guess
#     # p and c in Nerf frame
#     v = v./norm(v)
#     c = c0
#     for k = 1:10
#         Δc = density_gradient(nerf, c)
#         α = line_search(nerf, c, Δc)
#         c += α * Δc
#     end
#     return c
# end

# function line_search(nerf::NerfObject{T}, c::AbstractVector{T}, Δc::AbstractVector{T}) where T
#
#     return nothing
# end
