function set_input!(joint::JointConstraint{T,N,Nc}, input::AbstractVector) where {T,N,Nc}
    @assert length(input) == input_dimension(joint)
    # translational
    r_idx = SUnitRange(joint.minimal_index[1][1], joint.minimal_index[1][2])
    length(r_idx) > 0 && set_input!(joint.translational, input[r_idx])
    # rotational
    r_idx = SUnitRange(joint.minimal_index[2][1], joint.minimal_index[2][2])
    length(r_idx) > 0 && set_input!(joint.rotational, input[r_idx])
    return
end

function add_input!(joint::JointConstraint{T,N,Nc}, input::AbstractVector) where {T,N,Nc}
    @assert length(input) == input_dimension(joint)
    add_input!(joint.translational, input[SUnitRange(joint.minimal_index[1][1], joint.minimal_index[1][2])])
    add_input!(joint.rotational, input[SUnitRange(joint.minimal_index[2][1], joint.minimal_index[2][2])])
    return
end

@generated function input_jacobian_control(mechanism, joint::JointConstraint{T,N,Nc}, body::Body) where {T,N,Nc}
    relative = :(body.id == joint.parent_id ? :parent : :child)
    pbody = :(get_body(mechanism, joint.parent_id))
    cbody = :(get_body(mechanism, joint.child_id))
    rot = :(input_jacobian_control($relative, joint.translational, $pbody, $cbody, mechanism.input_scaling))
    tra = :(input_jacobian_control($relative, joint.rotational, $pbody, $cbody, mechanism.input_scaling))
    return :(hcat($rot, $tra))
end

function input_impulse!(joint::JointConstraint{T,N,Nc}, mechanism, clear::Bool=true) where {T,N,Nc}
    pbody = get_body(mechanism, joint.parent_id)
    cbody = get_body(mechanism, joint.child_id)
    input_impulse!(joint.translational, pbody, cbody, mechanism.input_scaling, clear)
    input_impulse!(joint.rotational, pbody, cbody, mechanism.input_scaling, clear)
    return
end

function input_dimension(joint::JointConstraint{T,N,Nc};
    ignore_floating_base::Bool=false) where {T,N,Nc}
    ignore_floating_base && (N == 0) && return 0
    N̄ = 0
    N̄ = input_dimension(joint.translational) + input_dimension(joint.rotational)
    return N̄
end