mutable struct Translational{T,Nλ,Nb,N,Nb½,N̄λ} <: Joint{T,Nλ,Nb,N,Nb½}
    axis::SVector{3,T} # translation axis in parent frame
    V3::Adjoint{T,SVector{3,T}} # in pbody's frame
    V12::SMatrix{2,3,T,6} # in pbody's frame
    vertices::NTuple{2,SVector{3,T}} # in pbody's & cbody's frames
    spring::T
    damper::T
    spring_offset::SVector{N̄λ,T}
    joint_limits::Vector{SVector{Nb½,T}} # lower and upper limits on the joint minimal coordinate angles
    spring_type::Symbol # the rotational springs can be :sinusoidal or :linear, if linear then we need joint_limits to avoid the 180° singularity.
    input::SVector{3,T}
end

function Translational{T,Nλ}(pbody::Node, cbody::Node;
        parent_vertex::AbstractVector=szeros(T,3), 
        child_vertex::AbstractVector=szeros(T,3), 
        axis::AbstractVector=szeros(T,3),
        spring=zero(T), 
        damper=zero(T), 
        spring_offset=szeros(T,3-Nλ),
        joint_limits=[szeros(T,0), szeros(T,0)],
        spring_type::Symbol = :linear, #TODO: re-implement sinusoidal?
    ) where {T,Nλ}

    vertices = (parent_vertex, child_vertex)
    V1, V2, V3 = orthogonal_rows(axis)
    V12 = [V1;V2]
    input = zeros(T,3)
    Nb½ = length(joint_limits[1])
    Nb = 2Nb½
    N̄λ = 3 - Nλ
    N = Nλ + 2Nb
    
    Translational{T,Nλ,Nb,N,Nb½,N̄λ}(axis, V3, V12, vertices, spring, damper, spring_offset, joint_limits, spring_type, input), pbody.id, cbody.id
end
