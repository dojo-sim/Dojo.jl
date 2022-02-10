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
    offset::SVector{3,T}
    nerf::NerfObject{T}

    function NerfContact111(body::Body{T}, normal::AbstractVector, nerf::NerfObject;
            p = szeros(T, 3), offset::AbstractVector = szeros(T, 3)) where T
        V1, V2, V3 = orthogonalcols(normal) # gives two plane vectors and the original normal axis
        A = [V1 V2 V3]
        Ainv = inv(A)
        ainv3 = Ainv[3,SA[1; 2; 3]]'
        new{Float64,2}(ainv3, p, offset, nerf)
    end
end

function constraint(mechanism, contact::ContactConstraint{T,N,Nc,Cs}) where {T,N,Nc,Cs<:Tuple{NerfContact{T,N}}}
    bound = contact.constraints[1]
    body = get_body(mechanism, contact.parent_id)
    x3, q3 = next_configuration(body.state, mechanism.timestep)
    # SVector{1,T}(bound.ainv3 * (x3 + vrotate(bound.p,q3) - bound.offset) - contact.primal[2][1])
    r = - bound.p - inv(q3) * x3 # position of the origin point in the frame of centered on P expressed in the body's frame
    c = bound.nerf.threshold - density(bound.nerf, r) # > 0
    SVector{1,T}(c - contact.primal[2][1])
end

@inline function constraint_jacobian_velocity(bound::NerfContact, x3::AbstractVector, q3::UnitQuaternion,
    x2::AbstractVector, v25::AbstractVector, q2::UnitQuaternion, ϕ25::AbstractVector, λ, timestep)
    # V = bound.ainv3 * timestep
    # Ω = bound.ainv3 * ∂vrotate∂q(bound.p, q3) * rotational_integrator_jacobian_velocity(q2, ϕ25, timestep)
    # return [V Ω]
    FiniteDiff.finite_difference_jacobian(
        vϕ -> -density(bound.nerf, -bound.p -inv(next_orientation(q2, vϕ[4:6], timestep)) * next_position(x2, vϕ[1:3], timestep)),
        [v25; ϕ25])
end

@inline function constraint_jacobian_configuration(bound::NerfContact, x3::AbstractVector, q3::UnitQuaternion,
    x2::AbstractVector, v25::AbstractVector, q2::UnitQuaternion, ϕ25::AbstractVector, λ, timestep)
    # X = bound.ainv3
    # Q = bound.ainv3 * ∂vrotate∂q(bound.p, q3)
    # return [X Q]
    FiniteDiff.finite_difference_jacobian(
        xq -> -density(bound.nerf, -bound.p -inv(UnitQuaternion(xq[4:7]..., false)) * xq[1:3]),
        [x3; vector(q3)])
end


@inline function impulse_map(bound::NerfContact, x::AbstractVector, q::UnitQuaternion, λ)
    # X = bound.ainv3
    # # q * ... is a rotation by quaternion q it is equivalent to Vmat() * Lmat(q) * Rmat(q)' * Vᵀmat() * ...
    # Q = - X * q * skew(bound.p - vrotate(bound.offset, inv(q)))
    X = FiniteDiff.finite_difference_jacobian(
        x -> -density(bound.nerf, -bound.p -inv(q)*x),
        x)
    Q = FiniteDiff.finite_difference_jacobian(
        q -> -density(bound.nerf, -bound.p -inv(UnitQuaternion(q..., false))*x),
        vector(q)) * LVᵀmat(q)
    return transpose([X Q])
end

@inline function force_mapping(bound::NerfContact, x::AbstractVector, q::UnitQuaternion)
    X = FiniteDiff.finite_difference_jacobian(
        x -> -density(bound.nerf, -bound.p -inv(q)*x),
        x)
    # X = bound.ainv3
    return X
end

@inline function set_matrix_vector_entries!(mechanism::Mechanism, matrix_entry::Entry, vector_entry::Entry,
    contact::ContactConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs<:Tuple{NerfContact{T,N}},N½}
    # ∇primal[dual .* primal - μ; g - s] = [diag(dual); -diag(0,1,1)]
    # ∇dual[dual .* primal - μ; g - s] = [diag(primal); -diag(1,0,0)]
    γ = contact.dual[2]
    s = contact.primal[2]

    ∇s = hcat(γ, -Diagonal(sones(N½)))
    ∇γ = hcat(s, -Diagonal(szeros(N½)))
    matrix_entry.value = hcat(∇s, ∇γ)

    # [-dual .* primal + μ; -g + s]
    vector_entry.value = vcat(-complementarityμ(mechanism, contact), -constraint(mechanism, contact))
    return
end

function get_nerf(; timestep::T=0.01, gravity=[0.0; 0.0; -9.81], cf::T=0.8, radius=0.1,
        contact::Bool=true, contact_type::Symbol=:impact) where T
    origin = Origin{T}(name=:origin)
    mass = 1.0
    bodies = [Sphere(radius, mass, name=:sphere)]
    joints = [JointConstraint(Floating(origin, bodies[1]), name = :floating_joint)]
    mechanism = Mechanism(origin, bodies, joints, timestep=timestep, gravity=gravity)

    if contact
        contact = [0,0,0.0]
        normal = [0,0,1.0]
        nerf = NerfObject111(radius, 0.05, 0.5)
        body = get_body(mechanism, :sphere)
        bound = NerfContact111(body, normal, nerf; p=szeros(T, 3), offset=szeros(T, 3))
        contacts = [ContactConstraint((bound, body.id, nothing); name=:nerf_contact)]
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

mech = get_nerf()
initialize!(mech, :nerf)
storage = simulate!(mech, 1.0, record=true, verbose=false)
visualize(mech, storage, vis=vis)

nerf0 = NerfObject111(0.1, 0.05, 0.5)
p1 = [0.05, 0.05, -0.0]
density(nerf0, p1)
plot(-1:0.01:1, [density(nerf0, [0,0,i]) for i in -1:0.01:1])


normal = [0,0,1.0]
NerfContact111(mech.bodies[1], normal, nerf0;
            p=szeros(3), offset=szeros(3))

# body0 =
# contact0 = NerfContact(body::Body{T}, normal::AbstractVector, nerf::NerfObject;
#         p = szeros(T, 3), offset::AbstractVector = szeros(T, 3))



#
# mech = get_humanoid()
# bound = mech.contacts[1].constraints[1]
# q3 = UnitQuaternion(rand(4)...)
# X = bound.ainv3
# Q0 = bound.ainv3 * ∂vrotate∂q(bound.p, q3) * LVᵀmat(q3)
# Q1 = - X * q3 * skew(bound.p - vrotate(bound.offset, inv(q3)))
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
