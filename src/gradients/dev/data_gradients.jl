################################################################################
# Data Jacobians
################################################################################
function joint_constraint_jacobian_body_data(mechanism::Mechanism, joint::JointConstraint{T,N},
        body::Body{T}) where {T,N}
    Nd = data_dim(body)
    ∇m = szeros(T,N,1)
    ∇J = szeros(T,N,6)
    ∇v15 = szeros(T,N,3)
    ∇ϕ15 = szeros(T,N,3)
    ∇z2 = -constraint_jacobian_configuration(mechanism, joint, body) *
        integrator_jacobian_configuration(body, mechanism.timestep, attjac=true)
    ∇g = [∇m ∇J ∇v15 ∇ϕ15 ∇z2]
    return ∇g
end

function joint_constraint_jacobian_joint_data(mechanism::Mechanism, joint::JointConstraint{T,N}) where {T,N}
    Nd = data_dim(joint)
    return szeros(T,N,Nd)
end


function body_constraint_jacobian_body_data(mechanism::Mechanism, body::Body{T}) where T
    Δt = mechanism.timestep
    Nd = data_dim(body)
    N = 6
    x1, v15, q1, ϕ15 = previous_configuration_velocity(body.state)
    x2, v25, q2, ϕ25 = current_configuration_velocity(body.state)
    x3, q3 = next_configuration(body.state, Δt)
    # Mass
    gravity = mechanism.gravity
    ∇m = [1 / Δt * (x2 - x1) + Δt/2 * gravity - 1 / Δt * (x3 - x2) + Δt/2 * gravity;
          szeros(T,3,1)]
    # Inertia
    ∇J = 4 / Δt * LVᵀmat(q2)' * LVᵀmat(q1) * ∂inertia(VLᵀmat(q1) * vector(q2))
    ∇J += 4 / Δt * LVᵀmat(q2)' * Tmat() * RᵀVᵀmat(q3) * ∂inertia(VLᵀmat(q2) * vector(q3))
    ∇J = [szeros(T,3,6); ∇J]

    # initial conditions: v15, ϕ15
    ∇v15 = body.mass * SMatrix{3,3,T,9}(Diagonal(sones(T,3)))
    ∇q1 = -4 / Δt * LVᵀmat(q2)' * ∂qLVᵀmat(body.inertia * VLᵀmat(q1) * vector(q2))
    ∇q1 += -4 / Δt * LVᵀmat(q2)' * LVᵀmat(q1) * body.inertia * ∂qVLᵀmat(vector(q2))
    ∇ϕ15 = ∇q1 * ∂integrator∂ϕ(q2, -ϕ15, Δt)
    ∇15 = [∇v15 szeros(T,3,3);
           szeros(T,3,3) ∇ϕ15]

    # current configuration: z2 = x2, q2
    # manipulator's equation contribution
    ∇tra_x2 = - 2 / Δt * body.mass * SMatrix{3,3,T,9}(Diagonal(sones(T,3)))
    ∇tra_q2 = szeros(T,3,3)
    ∇rot_x2 = szeros(T,3,3)
    ∇rot_q2 = -4 / Δt * VLᵀmat(q2) * LVᵀmat(q1) * body.inertia * VLᵀmat(q1)
    ∇rot_q2 += -4 / Δt * VLᵀmat(q2) * Tmat() * RᵀVᵀmat(q3) * body.inertia * ∂qVLᵀmat(vector(q3))
    ∇rot_q2 += -4 / Δt * ∂qVLᵀmat(LVᵀmat(q1) * body.inertia * VLᵀmat(q1) * vector(q2) + Tmat() * RᵀVᵀmat(q3) * body.inertia * VLᵀmat(q2) * vector(q3))
    ∇rot_q2 *= LVᵀmat(q2)
    ∇z2 = [∇tra_x2 ∇tra_q2;
           ∇rot_x2 ∇rot_q2]
    # TODO
    # joint constraints impulses contribution
    ∇z2 += impulse_map_child_jacobian_child(joint, szeros(3), one(UnitQuaternion), current_configuration(body.state))

    # TODO
    # # constact constraints impulses contribution
    # ∇z2 +=
    # TODO
    # # spring and damper impulses contribution
    # ∇z2 +=
    return [∇m ∇J ∇15 ∇z2]
end

function body_constraint_jacobian_joint_data(mechanism::Mechanism{T}, body::Body{T},
    joint::JointConstraint{T}) where {T}
    Δt = mechanism.timestep
    Nd = data_dim(joint)
    N = 6
    x1, v15, q1, ϕ15 = previous_configuration_velocity(body.state)
    x2, v25, q2, ϕ25 = current_configuration_velocity(body.state)
    x3, q3 = next_configuration(body.state, Δt)
    ∇u = Diagonal(SVector{6,T}(1,1,1,2,2,2)) * input_jacobian_control(mechanism, joint, body)
    ∇spring = szeros(T,6,1)#springforce(mechanism, joint, body, unitary=true)
    ∇damper = szeros(T,6,1)#damperforce(mechanism, joint, body, unitary=true)
    return [∇u ∇spring ∇damper]
end

function body_constraint_jacobian_contact_data(mechanism::Mechanism, body::Body{T},
    contact::ContactConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs<:Tuple{NonlinearContact{T,N}},N½}
    Nd = data_dim(contact)
    bound = contact.constraints[1]
    offset = bound.offset
    x3, q3 = next_configuration(body.state, mechanism.timestep)
    γ = contact.γsol[2]

    ∇cf = szeros(T,3,1)

    X = force_mapping(bound)
    # this what we differentiate: Qᵀγ = - skew(p - vrotate(offset, inv(q3))) * VRmat(q3) * LᵀVᵀmat(q3) * X' * γ
    ∇p = - ∂pskew(VRmat(q3) * LᵀVᵀmat(q3) * X' * γ)
    ∇off = - ∂pskew(VRmat(q3) * LᵀVᵀmat(q3) * X' * γ) * -∂vrotate∂p(offset, inv(q3))

    ∇X = szeros(T,3,Nd)
    ∇Q = [∇cf ∇p ∇off]
    return [∇X; ∇Q]
end


function contact_constraint_jacobian_contact_data(mechanism::Mechanism, contact::ContactConstraint{T,N,Nc,Cs,N½},
        body::Body{T}) where {T,N,Nc,Cs<:Tuple{NonlinearContact{T,N}},N½}
    Nd = data_dim(contact)
    bound = contact.constraints[1]
    p = bound.p
    offset = bound.offset
    x2, v25, q2, ϕ25 = current_configuration_velocity(body.state)
    x3, q3 = next_configuration(body.state, mechanism.timestep)
    s = contact.ssol[2]
    γ = contact.γsol[2]

    ∇cf = SA[0,γ[1],0,0]
    ∇off = [-bound.ainv3; szeros(T,1,3); -bound.Bx * skew(vrotate(ϕ25, q3))]
    ∇p = [bound.ainv3 * ∂vrotate∂p(bound.p, q3); szeros(T,1,3); bound.Bx * skew(vrotate(ϕ25, q3)) * ∂vrotate∂p(bound.p, q3)]

    ∇compμ = szeros(T,N½,Nd)
    ∇g = [∇cf ∇p ∇off]
    return [∇compμ; ∇g]
end

function contact_constraint_jacobian_body_data(mechanism::Mechanism, contact::ContactConstraint{T,N,Nc,Cs,N½},
        body::Body{T}) where {T,N,Nc,Cs,N½}
    Nd = data_dim(body)
    ∇compμ = szeros(T,N½,Nd)
    ∇m = szeros(T,N½,1)
    ∇J = szeros(T,N½,6)
    ∇z1 = szeros(T,N½,6)
    ∇z3 = constraint_jacobian_configuration(mechanism, contact, body)
    ∇z2 = ∇z3 * ∂i∂z(body, mechanism.timestep, attjac=true) # 4x7 * 7x6 = 4x6
    ∇g = [∇m ∇J ∇z1 ∇z2]
    return [∇compμ; ∇g]
end

################################################################################
# System Data Jacobians
################################################################################
function create_data_matrix(joints::Vector{<:JointConstraint}, bodies::Vector{<:Body},
        contacts::Vector{<:ContactConstraint})
    nodes = [joints; bodies; contacts]
    A = data_adjacency_matrix(joints, bodies, contacts)
    dimrow = length.(nodes)
    dimcol = data_dim.(nodes)
    D = data_matrix(A, dimrow, dimcol)
    return D
end

function jacobian_data!(data_system::System, mech::Mechanism)
    jacobian_contact_data!(data_system, mech::Mechanism)
    jacobian_body_data!(data_system, mech::Mechanism)
    jacobian_joint_data!(data_system, mech::Mechanism)
    return nothing
end

function jacobian_contact_data!(data_system::System, mechanism::Mechanism{T}) where {T}
    # ∂body∂ineqcdata
    for contact in mechanism.contacts
        pbody = getbody(mechanism, contact.parentid)
        data_system.matrix_entries[pbody.id, contact.id].value += body_constraint_jacobian_contact_data(mechanism, pbody, contact)
    end
    # ∂contact∂contactdata
    for contact in mechanism.contacts
        pbody = get_body(mechanism, contact.parentid)
        data_system.matrix_entries[contact.id, contact.id].value += contact_constraint_jacobian_contact_data(mechanism, contact, pbody)
    end
    return nothing
end

function jacobian_joint_data!(data_system::System, mechanism::Mechanism{T}) where T
    # ∂joint∂jointdata
    for joint in mechanism.joints
        data_system.matrix_entries[joint.id, joint.id].value += joint_constraint_jacobian_joint_data(mechanism, joint)
    end
    # ∂body∂jointdata
    # TODO adapt this to handle cycles
    for body in mechanism.bodies
        for joint in mechanism.joints
            if (body.id == joint.parentid) || (body.id ∈ joint.childids)
                data_system.matrix_entries[body.id, joint.id].value += body_constraint_jacobian_joint_data(mechanism, body, joint)
            end
        end
    end
    return nothing
end

function jacobian_body_data!(data_system::System, mechanism::Mechanism{T}) where T
    # ∂joint∂bodydata
    # TODO adapt this to handle cycles
    for body in mechanism.bodies
        for joint in mechanism.joints
            if (body.id == joint.parentid) || (body.id ∈ joint.childids)
                data_system.matrix_entries[joint.id, body.id].value += joint_constraint_jacobian_body_data(mechanism, joint, body)
            end
        end
    end
    # ∂body∂bodydata
    for body in mechanism.bodies
        data_system.matrix_entries[body.id, body.id].value += body_constraint_jacobian_body_data(mechanism, body)
    end
    # ∂contact∂bodydata
    for contact in mechanism.contacts
        pbody = get_body(mechanism, contact.parentid)
        data_system.matrix_entries[contact.id, pbody.id].value += contact_constraint_jacobian_body_data(mechanism, contact, pbody)
    end
    return nothing
end

#
# function ∂body∂z_local(body::Body{T}, Δt::T; attjac::Bool=true) where T
#     state = body.state
#     q2 = state.q2[1]
#     # ϕ25 = state.ϕsol[2]
#     Z3 = szeros(T,3,3)
#     Z34 = szeros(T,3,4)
#     ZT = attjac ? szeros(T,6,6) : szeros(T,6,7)
#     ZR = szeros(T,7,6)
#
#     x1, q1 = previous_configuration(state)
#     x2, q2 = current_configuration(state)
#     x3, q3 = next_configuration(state, Δt)
#
#     AposT = [-I Z3]
#     # AvelT = [Z3 -I*body.mass] # solving for impulses
#
#     AposR = [-∂integrator∂q(q2, ϕ25, Δt, attjac = attjac) szeros(4,3)]
#
#     rot_q1(q) = -4 / Δt * LVᵀmat(q2)' * Lmat(UnitQuaternion(q..., false)) * Vᵀmat() * body.inertia * Vmat() * Lmat(UnitQuaternion(q..., false))' * vector(q2)
#     rot_q2(q) = -4 / Δt * LVᵀmat(UnitQuaternion(q..., false))' * Tmat() * Rmat(getq3(UnitQuaternion(q..., false), state.ϕsol[2], Δt))' * Vᵀmat() * body.inertia * Vmat() * Lmat(UnitQuaternion(q..., false))' * vector(getq3(UnitQuaternion(q..., false), state.ϕsol[2], Δt)) + -4 / Δt * LVᵀmat(UnitQuaternion(q..., false))' * Lmat(getq3(UnitQuaternion(q..., false), -state.ϕ15, Δt)) * Vᵀmat() * body.inertia * Vmat() * Lmat(getq3(UnitQuaternion(q..., false), -state.ϕ15, Δt))' * q
#
#     # dynR_ϕ15 = -1.0 * FiniteDiff.finite_difference_jacobian(rot_q1, vector(q1)) * ∂integrator∂ϕ(q2, -state.ϕ15, Δt)
#     dynR_q2 = FiniteDiff.finite_difference_jacobian(rot_q2, vector(q2))
#     AvelR = attjac ? [dynR_q2 * LVᵀmat(q2) dynR_ϕ15] : [dynR_q2 dynR_ϕ15]
#
#     return [[AposT;AvelT] ZT;
#              ZR [AposR;AvelR]]
# end
