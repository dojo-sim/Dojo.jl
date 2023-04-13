################################################################################
# Control Input
################################################################################

function input_impulse!(joint::Translational{T}, 
    pbody::Node, cbody::Node,
    input_scaling::T, clear::Bool) where T

    xa, qa = current_configuration(pbody.state)
    xb, qb = current_configuration(cbody.state)
  
    input = joint.input * input_scaling
    Ta = impulse_transform(:parent, joint, xa, qa, xb, qb)
    Tb = impulse_transform(:child, joint, xa, qa, xb, qb)
    JFaw = Ta[SA[1;2;3],SA[1;2;3]] * input
    Jτaa = Ta[SA[4;5;6],SA[1;2;3]] * input
    JFbw = Tb[SA[1;2;3],SA[1;2;3]] * input
    Jτbb = Tb[SA[4;5;6],SA[1;2;3]] * input

    pbody.state.JF2 += JFaw
    pbody.state.Jτ2 += Jτaa/2
    cbody.state.JF2 += JFbw
    cbody.state.Jτ2 += Jτbb/2
    
    clear && (joint.input = szeros(T,3))
    return
end

################################################################################
# Control Jacobian
################################################################################

function input_jacobian_control(relative::Symbol, 
    joint::Translational, 
    xa::AbstractVector, qa::Quaternion,
    xb::AbstractVector, qb::Quaternion,
    input_scaling)

    Ta = impulse_transform(relative, joint, xa, qa, xb, qb)
    X = Ta[SA[1;2;3],SA[1;2;3]]
    Q = 0.5 * Ta[SA[4;5;6],SA[1;2;3]]
    return [X; Q] * input_scaling
end

function input_jacobian_configuration(relative::Symbol, 
    joint::Translational{T}, 
    xa::AbstractVector, qa::Quaternion,
    xb::AbstractVector, qb::Quaternion) where T

    # d[Faw;2τaa]/d[xa,qa]
    ∇aa = impulse_transform_jacobian(:parent, relative, joint, xa, qa, xb, qb, joint.input)
    FaXa = ∇aa[SA[1;2;3], SA[1;2;3]]
    FaQa = ∇aa[SA[1;2;3], SA[4;5;6]]
    τaXa = 0.5 * ∇aa[SA[4;5;6], SA[1;2;3]]
    τaQa = 0.5 * ∇aa[SA[4;5;6], SA[4;5;6]]

    # d[Fbw;2τbb]/d[xa,qa]
    ∇ba = impulse_transform_jacobian(:child, relative, joint, xa, qa, xb, qb, joint.input)
    FbXa = ∇ba[SA[1;2;3],SA[1;2;3]]
    FbQa = ∇ba[SA[1;2;3],SA[4;5;6]]
    τbXa = 0.5 * ∇ba[SA[4;5;6],SA[1;2;3]]
    τbQa = 0.5 * ∇ba[SA[4;5;6],SA[4;5;6]]

    return FaXa, FaQa, τaXa, τaQa, FbXa, FbQa, τbXa, τbQa
end
