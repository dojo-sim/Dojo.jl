""" 
    Rotational{T} <: Joint{T}

    constraint limiting rotational degrees of freedom

    axis: rotation axis in parent offset frame
    axis_mask1: rotation axis mask in pbody's frame
    axis_mask2: rotation axis mask in pbody's frame
    axis_mask3: rotation axis mask in pbody's frame
    orientation_offset: rotation axis offset from pbody's frame
    spring :coefficient for joint spring
    damper: coefficient for joint damper
    spring_offset: nominal joint configuration
    joint_limits: lower and upper limits on the joint configuration
    spring_type: can be :linear or :sinusoidal (currently not implemented), if linear then we need joint_limits to avoid the 360° singularity.
    input: external impulse torque
"""
mutable struct Rotational{T,Nλ,Nb,N,Nb½,N̄λ} <: Joint{T,Nλ,Nb,N,Nb½}
    axis::SVector{3,T} # rotation axis in parent offset frame
    axis_mask1::Adjoint{T,SVector{3,T}} # in pbody's frame
    axis_mask2::Adjoint{T,SVector{3,T}} # in pbody's frame
    axis_mask3::Adjoint{T,SVector{3,T}} # in pbody's frame
    orientation_offset::Quaternion{T} # in pbody's frame
    spring::T
    damper::T
    spring_offset::SVector{N̄λ,T}
    joint_limits::Vector{SVector{Nb½,T}} # lower and upper limits on the joint minimal coordinate angles
    spring_type::Symbol # the rotational springs can be :linear or :sinusoidal (currently not implemented), if linear then we need joint_limits to avoid the 360° singularity.
    input::SVector{3,T} # joint input, i.e., a torque, not an impulse
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, joint::Rotational{T,N}) where {T,N}
    summary(io, joint)
    println(io,"")
    println(io, " axis:          "*string(joint.axis))
    println(io, " axis_mask1:    "*string(joint.axis_mask1))
    println(io, " axis_mask2:    "*string(joint.axis_mask2))
    println(io, " axis_mask3:    "*string(joint.axis_mask3))
    println(io, " orientation_offset:   "*string(joint.orientation_offset))
    println(io, " spring:        "*string(joint.spring))
    println(io, " damper:        "*string(joint.damper))
    println(io, " spring_offset: "*string(joint.spring_offset))
    println(io, " joint_limits:  "*string(joint.joint_limits))
    println(io, " spring_type:   "*string(joint.spring_type))
    println(io, " input:         "*string(joint.input))
end

function Rotational{T,Nλ}(pbody::Node, cbody::Node;
        axis::AbstractVector=szeros(T,3), 
        orientation_offset::Quaternion=one(Quaternion{T}),
        spring=zero(T), 
        damper=zero(T), 
        spring_offset=szeros(T,3-Nλ),
        joint_limits=[szeros(T,0), szeros(T,0)],
        spring_type::Symbol=:linear) where {T,Nλ}
    V1, V2, V3 = orthogonal_rows(axis)
    input = zeros(T,3)
    Nb½ = length(joint_limits[1])
    Nb = 2Nb½
    N̄λ = 3 - Nλ
    N = Nλ + 2Nb
    Rotational{T,Nλ,Nb,N,Nb½,N̄λ}(axis, V1, V2, V3, orientation_offset, spring, damper, spring_offset, joint_limits, spring_type, input), pbody.id, cbody.id
end

