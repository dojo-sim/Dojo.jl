# gravity
get_gravity(g::T) where T <: Real = SVector{3,T}([0.0; 0.0; g])
get_gravity(g::Vector{T}) where T = SVector{3,T}(g)
get_gravity(g::SVector) = g