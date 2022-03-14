################################################################################
# Control Input
################################################################################

function input_impulse!(joint::Rotational{T}, pbody::Node, cbody::Node, 
    timestep::T, clear::Bool) where T

    τ = joint.input
    xa, qa = current_configuration(pbody.state)
    xb, qb = current_configuration(cbody.state)

    pbody.state.τ2 += -τ
    cbody.state.τ2 += vector_rotate(vector_rotate(τ, qa),inv(qb))
    clear && (joint.input = szeros(T,3))
    return
end

################################################################################
# Control Jacobian
################################################################################

function input_jacobian_control(relative::Symbol, 
    joint::Rotational{T}, 
    xa::AbstractVector, qa::Quaternion,
    xb::AbstractVector, qb::Quaternion) where T

    if relative == :parent
        BFa = szeros(T, 3, 3)
        Bτa = -I
        return [BFa; Bτa]
    elseif relative == :child 
        BFb = szeros(T, 3, 3)
        Bτb = rotation_matrix(inv(qb)) * rotation_matrix(qa)
        return [BFb; Bτb]
    end
end

function input_jacobian_configuration(relative::Symbol, 
    joint::Rotational{T}, 
    xa::AbstractVector, qa::Quaternion,
    xb::AbstractVector, qb::Quaternion) where T

    τ = joint.input

    if relative == :parent
        FaXa = szeros(T,3,3)
        FaQa = szeros(T,3,4)
        τaXa = szeros(T,3,3)
        τaQa = szeros(T,3,4)
        FbXa = szeros(T,3,3)
        FbQa = szeros(T,3,4)
        τbXa = szeros(T,3,3)
        τbQa = rotation_matrix(inv(qb)) * ∂rotation_matrix∂q(qa, τ)#*LVᵀmat(qa)
        return FaXa, FaQa, τaXa, τaQa, FbXa, FbQa, τbXa, τbQa
    elseif relative == :child 
        FaXb = szeros(T,3,3)
        FaQb = szeros(T,3,4)
        τaXb = szeros(T,3,3)
        τaQb = szeros(T,3,4)
        FbXb = szeros(T,3,3)
        FbQb = szeros(T,3,4)
        τbXb = szeros(T,3,3)
        τbQb = ∂rotation_matrix_inv∂q(qb, rotation_matrix(qa)*τ)#*LVᵀmat(qb)
        return FaXb, FaQb, τaXb, τaQb, FbXb, FbQb, τbXb, τbQb
    end
end