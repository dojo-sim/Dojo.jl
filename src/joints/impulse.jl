################################################################################
# Impulse Transform
################################################################################
function impulse_transform(relative::Symbol, joint::Rotational{T}, xa::AbstractVector,
        qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where {T}
    (relative == :parent) && (return impulse_transform_parent(joint, xa, qa, xb, qb))
    (relative == :child) && (return impulse_transform_child(joint, xa, qa, xb, qb))
end

################################################################################
# Derivatives
################################################################################
function impulse_transform_jacobian(relative::Symbol, jacobian_relative::Symbol, joint::Translational{T,Nλ,0},
        xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, p) where {T,Nλ}
    (relative == :parent) && (jacobian_relative == :parent) && (return impulse_transform_parent_jacobian_parent(joint, xa, qa, xb, qb, p))
    (relative == :parent) && (jacobian_relative == :child) && (return impulse_transform_parent_jacobian_child(joint, xa, qa, xb, qb, p))
    (relative == :child) && (jacobian_relative == :parent) && (return impulse_transform_child_jacobian_parent(joint, xa, qa, xb, qb, p))
    (relative == :child) && (jacobian_relative == :child) && (return impulse_transform_child_jacobian_child(joint, xa, qa, xb, qb, p))
    return
end
