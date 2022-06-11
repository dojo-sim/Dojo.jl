struct DensityFieldNormalizer{T}
    position_offset::Vector{T}
    orientation_offset::Quaternion{T}
    geometry_scaling::T
    density_low::T
    density_high::T
    density_scale::T
end

function particle_normalization(x; normalizer::DensityFieldNormalizer=DensityFieldNormalizer())

    orientation_offset = normalizer.orientation_offset
    geometry_scaling = normalizer.geometry_scaling
    position_offset = normalizer.position_offset

    R = rotationmatrix(orientation_offset)
    x̄ = R * (geometry_scaling * x) + position_offset
    return x̄
end

function DensityFieldNormalizer(;nerf::Symbol=:identity)
    if nerf ==:identity
        position_offset = [0.0, 0.0, 0.0]
        orientation_offset = Quaternion(1,0,0,0.0)
        geometry_scaling = 1.0
        density_low = 0.0
        density_high = 1.0
        density_scale = 1.0
    elseif nerf ==:bunny
        position_offset = [0.0, 0.0, 0.0]
        orientation_offset = Quaternion(1,0,0,0.0)
        geometry_scaling = 1.0
        density_low = 1.0
        density_high = 105.0
        density_scale = 0.1
    elseif nerf == :bluesoap
        position_offset = [0.0, 2.0, 4.5]
        q = [1.0, 0.0, 0.0, 0.0]
        orientation_offset = Quaternion(normalize!(q))
        geometry_scaling = 2.5
        density_low = 60.0
        density_high = 90.0
        density_scale = 0.1
    elseif nerf == :halfsoap
        position_offset = [-2.5,  1.75,  4.35]
        q = [1.0, 0.0, 0.0, 0.0]
        orientation_offset = Quaternion(normalize!(q))
        geometry_scaling = 2.5
        density_low = 20.0
        density_high = 40.0
        density_scale = 0.2
    else
        @warn "unknown nerf object"
    end
    return DensityFieldNormalizer(
        position_offset,
        orientation_offset,
        geometry_scaling,
        density_low,
        density_high,
        density_scale
        )
end
