################################################################################
# Data Jacobians
################################################################################
function ∂joint∂body_data(mechanism::Mechanism, joint::JointConstraint{T,N},
        body::Body{T}) where {T,N}
    Nd = data_dim(body)
    ∇m = szeros(T,N,1)
    ∇J = szeros(T,N,6)
    ∇z1 = szeros(T,N,6)
    ∇z2 = constraint_jacobian_configuration(mechanism, joint, body) * ∂i∂z(body, mechanism.timestep, attjac=true)
    ∇g = [∇m ∇J ∇z1 ∇z2]
    return ∇g
end

function ∂joint∂joint_data(mechanism::Mechanism, joint::JointConstraint{T,N}) where {T,N}
    Nd = data_dim(joint)
    return szeros(T,N,Nd)
end


function ∂body∂body_data(mechanism::Mechanism, body::Body{T}) where T
    Nd = data_dim(body)
    N = 6
    x1, v15, q1, ϕ15 = fullargs1(body.state)
    x2, v25, q2, ϕ25 = current_configuration_velocity(body.state)
    x3, q3 = next_configuration(body.state, timestep)
    # Mass
    ezg = SA{T}[0; 0; -mechanism.g]
    ∇m = [- 1 / timestep * (x2 - x1) + timestep/2 * ezg + 1 / timestep * (x3 - x2) + timestep/2 * ezg;
          szeros(T,3,1)]
    # Inertia
    ∇J = -4 / timestep * LVᵀmat(q2)' * LVᵀmat(q1) * ∂inertia(VLᵀmat(q1) * vector(q2))
    ∇J += -4 / timestep * LVᵀmat(q2)' * Tmat() * RᵀVᵀmat(q3) * ∂inertia(VLᵀmat(q2) * vector(q3))
    ∇J = [szeros(T,3,6); ∇J]

    # configuration 1: z1 = x1, q1
    ∇x1 = 1 / timestep * body.m * SMatrix{3,3,T,9}(Diagonal(sones(T,3)))
    ∇q1 = -4 / timestep * LVᵀmat(q2)' * ∂qLVᵀmat(body.J * VLᵀmat(q1) * vector(q2))
    ∇q1 += -4 / timestep * LVᵀmat(q2)' * LVᵀmat(q1) * body.J * ∂qVLᵀmat(vector(q2))
    ∇q1 *= LVᵀmat(q1) # attjac
    ∇z1 = [∇x1 szeros(T,3,3);
           szeros(T,3,3) ∇q1]

    # configuration 2: z2 = x2, q2
    # TODO
    ∇z2 = szeros(T,N,6)

    return [∇m ∇J ∇z1 ∇z2]
end

function ∂body∂joint_data(mechanism::Mechanism{T}, joint::JointConstraint{T},
        body::Body{T}) where T
    Nd = data_dim(joint)
    N = 6
    x1, v15, q1, ϕ15 = fullargs1(body.state)
    x2, v25, q2, ϕ25 = current_configuration_velocity(body.state)
    x3, q3 = next_configuration(body.state, timestep)
    ∇u = Diagonal(SVector{6,T}(-1,-1,-1,-2,-2,-2)) * input_jacobian_control(mechanism, joint, body)
    ∇spring = apply_spring(mechanism, joint, body, unitary=true)
    ∇damper = apply_damper(mechanism, joint, body, unitary=true)
    # TODO
    nu = control_dimension(joint)
    ∇spring_offset = szeros(T,N,nu)
    return [∇u ∇spring ∇damper ∇spring_offset]
end

function ∂body∂contact_data(mechanism::Mechanism, contact::ContactConstraint{T,N,Nc,Cs,N½},
        body::Body{T}) where {T,N,Nc,Cs<:Tuple{NonlinearContact{T,N}},N½}
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


function ∂contact∂contact_data(mechanism::Mechanism, contact::ContactConstraint{T,N,Nc,Cs,N½},
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

function ∂contact∂body_data(mechanism::Mechanism, contact::ContactConstraint{T,N,Nc,Cs,N½},
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
function ∂contact_data!(data_system::System, mechanism::Mechanism{T}) where T
    # ∂body∂contactdata
    for contact in mechanism.contacts
        pbody = get_body(mechanism, contact.parentid)
        data_system.matrix_entries[pbody.id, contact.id].value += ∂body∂contact_data(mechanism, contact, pbody)
    end
    # ∂contact∂contactdata
    for contact in mechanism.contacts
        pbody = get_body(mechanism, contact.parentid)
        data_system.matrix_entries[contact.id, contact.id].value += ∂contact∂contact_data(mechanism, contact, pbody)
    end
    return nothing
end

function ∂joint_data!(data_system::System, mechanism::Mechanism{T}) where T
    # ∂joint∂jointdata
    for joint in mechanism.joints
        data_system.matrix_entries[joint.id, joint.id].value += ∂joint∂joint_data(mechanism, joint)
    end
    # ∂body∂jointdata
    # TODO adapt this to handle cycles
    for body in mechanism.bodies
        for joint in mechanism.joints
            if (body.id == joint.parentid) || (body.id ∈ joint.childids)
                data_system.matrix_entries[body.id, joint.id].value += ∂body∂joint_data(mechanism, joint, body)
            end
        end
    end
    return nothing
end

function ∂body_data!(data_system::System, mechanism::Mechanism{T}) where T
    # ∂joint∂bodydata
    # TODO adapt this to handle cycles
    for body in mechanism.bodies
        for joint in mechanism.joints
            if (body.id == joint.parentid) || (body.id ∈ joint.childids)
                data_system.matrix_entries[joint.id, body.id].value += ∂joint∂body_data(mechanism, joint, body)
            end
        end
    end
    # ∂body∂bodydata
    for body in mechanism.bodies
        data_system.matrix_entries[body.id, body.id].value += ∂body∂body_data(mechanism, body)
    end
    # ∂contact∂bodydata
    for contact in mechanism.contacts
        pbody = get_body(mechanism, contact.parentid)
        data_system.matrix_entries[contact.id, pbody.id].value += ∂contact∂body_data(mechanism, contact, pbody)
    end
    return nothing
end
