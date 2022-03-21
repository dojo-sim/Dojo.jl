################################################################################
# Impulse Transform
################################################################################
function impulse_transform(relative::Symbol, joint::Translational{T}, xa::AbstractVector,
        qa::Quaternion, xb::AbstractVector, qb::Quaternion) where T
    X, Q = displacement_jacobian_configuration(relative, joint, xa, qa, xb, qb, attjac=true)
    Diagonal([sones(3); 0.5 * sones(3)]) * transpose([X Q])
end

function impulse_transform(relative::Symbol, joint::Rotational{T}, xa::AbstractVector,
        qa::Quaternion, xb::AbstractVector, qb::Quaternion) where T
    X = szeros(T,3,3)
    if relative == :parent
        Q = -Diagonal(sones(T,3))
        # Q = -rotation_matrix(inv(qa))
    elseif relative == :child
        Q = rotation_matrix(inv(qb) * qa)
        # Q = rotation_matrix(inv(qb))
    end
    [X; Q]
end


################################################################################
# Derivatives
################################################################################
function impulse_map_jacobian(relative::Symbol, jacobian::Symbol, joint::Joint,
        pbody::Node{T}, cbody::Node{T}, λ) where T
    # ∂(G*λ)/∂(x,q)
    p = impulse_projector(joint) * λ
    impulse_transform_jacobian(relative, jacobian,
        joint,
        current_configuration(pbody.state)...,
        current_configuration(cbody.state)...,
        p)
end
