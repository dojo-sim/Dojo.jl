"""
    CapsuleCapsuleCollision 

    collision between two capsules

    origin_parent: radius of contact for parent body 
    origin_child: radius of contact for child body
    orientation_parent: radius of contact for parent body 
    orientation_child: radius of contact for child body
    radius_parent: radius of contact for parent body 
    radius_child: radius of contact for child body
    height_parent: radius of contact for parent body 
    height_child: radius of contact for child body
    ip: InteriorPoint
"""
mutable struct CapsuleCapsuleCollision{T,O,I,OI} <: Collision{T,O,I,OI}
    origin_parent::SVector{I,T}
    origin_child::SVector{I,T}
    orientation_parent::Quaternion{T}
    orientation_child::Quaternion{T}
    radius_parent::T
    radius_child::T
    height_parent::T
    height_child::T
    ip::InteriorPoint{T} #TODO: parameterize correctly
end

function CapsuleCapsuleCollision(radius_parent::T, height_parent::T, radius_child::T, height_child::T;
    origin_parent=zeros(3),
    orientation_parent=Quaternion(1.0, 0.0, 0.0, 0.0), 
    origin_child=zeros(3),
    orientation_child=Quaternion(1.0, 0.0, 0.0, 0.0),
    contact=:nonlinear,
    regularization=1.0e-5,
    r_func=r_capsule_capsule_func, 
    rw_func=rw_capsule_capsule_func, 
    rθ_func=rθ_capsule_capsule_func,
    opts=InteriorPointOptions(
                undercut=5.0,
                γ_reg=0.1,
                r_tol=1e-5,
                κ_tol=1e-5,  
                max_ls=50,
                ϵ_min=0.05,
                diff_sol=true,
                verbose=true)) where T

    # solver 
    nz = 4
    ny = 6
    nw = nz + ny 
    nθ = 12

    idx = IndicesOptimization(
        nw, nw, 
        [collect(1:4), collect(7:10)], [collect(1:4), collect(7:10)],
        Vector{Vector{Vector{Int}}}(), Vector{Vector{Vector{Int}}}(),
        collect(1:6),
        collect(7:10),
        Vector{Int}(),
        Vector{Vector{Int}}(),
        collect(7:10),
    )

    ip = interior_point(zeros(nw), zeros(nθ);
        s = Euclidean(nw),
        idx = idx,
        r! = r_func, 
        rz! = (rw, w, θ) -> rw_func(rw, w, θ, [regularization]), 
        rθ! = rθ_func,
        r  = zeros(idx.nΔ),
        rz = zeros(idx.nΔ, idx.nΔ),
        rθ = zeros(idx.nΔ, nθ),
        opts=opts)

    I = 3
    if contact == :nonlinear 
        O = 2 
    elseif contact == :linear 
        O = 4
    elseif contact == :impact 
        O = 0
    end

    CapsuleCapsuleCollision{T,O,I,O*I}(
        origin_parent,
        origin_child, 
        orientation_parent, 
        orientation_child, 
        radius_parent, 
        radius_child, 
        height_parent, 
        height_child,
        ip,
    )
end

include(joinpath(@__DIR__, "codegen.jl"))
path_collisions = @get_scratch!("collisions")
path_expr = joinpath(path_collisions, "capsule_capsule" * ".jld2")
@load path_expr r_capsule_capsule_func rw_capsule_capsule_func rθ_capsule_capsule_func
