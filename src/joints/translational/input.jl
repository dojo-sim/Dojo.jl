################################################################################
# Control Input
################################################################################

function input_impulse!(joint::Translational{T}, 
    pbody::Node, cbody::Node,
    timestep::T, clear::Bool) where T

    xa, qa = current_configuration(pbody.state)
    xb, qb = current_configuration(cbody.state)
  
    input = joint.input
    Ta = impulse_transform(:parent, joint, xa, qa, xb, qb)
    Tb = impulse_transform(:child, joint, xa, qa, xb, qb)
    Faw = Ta[1:3,1:3] * input
    τaa = Ta[4:6,1:3] * input
    Fbw = Tb[1:3,1:3] * input
    τbb = Tb[4:6,1:3] * input

    pbody.state.F2 += Faw
    pbody.state.τ2 += τaa/2
    cbody.state.F2 += Fbw
    cbody.state.τ2 += τbb/2
    
    clear && (joint.input = szeros(T,3))
    return
end

################################################################################
# Control Jacobian
################################################################################

function input_jacobian_control(relative::Symbol, 
    joint::Translational, 
    xa::AbstractVector, qa::Quaternion,
    xb::AbstractVector, qb::Quaternion) where T

    Ta = impulse_transform(relative, joint, xa, qa, xb, qb)
    X = Ta[1:3,1:3]
    Q = 0.5 * Ta[4:6,1:3]
    return [X; Q]
end

function input_jacobian_configuration(relative::Symbol, 
    joint::Translational{T}, 
    xa::AbstractVector, qa::Quaternion,
    xb::AbstractVector, qb::Quaternion) where T

    # d[Faw;2τaa]/d[xa,qa]
    ∇aa = impulse_transform_jacobian(:parent, relative, joint, xa, qa, xb, qb, joint.input)
    FaXa = ∇aa[1:3, 1:3]
    FaQa = ∇aa[1:3, 4:6]
    τaXa = 0.5 * ∇aa[4:6, 1:3]
    τaQa = 0.5 * ∇aa[4:6, 4:6]

    # d[Fbw;2τbb]/d[xa,qa]
    ∇ba = impulse_transform_jacobian(:child, relative, joint, xa, qa, xb, qb, joint.input)
    FbXa = ∇ba[1:3,1:3]
    FbQa = ∇ba[1:3,4:6]
    τbXa = 0.5 * ∇ba[4:6,1:3]
    τbQa = 0.5 * ∇ba[4:6,4:6]

    return FaXa, FaQa, τaXa, τaQa, FbXa, FbQa, τbXa, τbQa
end
