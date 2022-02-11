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


function body_constraint_jacobian_body_data(mechanism::Mechanism, body::Body{T}) where T
    Δt = mechanism.timestep
    N = 6
    x1, v15, q1, ϕ15 = previous_configuration_velocity(body.state)
    x2, v25, q2, ϕ25 = current_configuration_velocity(body.state)
    x3, q3 = next_configuration(body.state, Δt)
    # Mass
    gravity = mechanism.gravity
    ∇m = [1 / Δt * (x2 - x1) + Δt/2 * gravity - 1 / Δt * (x3 - x2) + Δt/2 * gravity;
          szeros(T,3,1)]
    # Inertia
    ∇J = 2 / Δt * LVᵀmat(q2)' * LVᵀmat(q1) * ∂inertia(VLᵀmat(q1) * vector(q2))
    ∇J += 2 / Δt * LVᵀmat(q2)' * Tmat() * RᵀVᵀmat(q3) * ∂inertia(VLᵀmat(q2) * vector(q3))
    ∇J = [szeros(T,3,6); ∇J]

    # initial conditions: v15, ϕ15
    ∇v15 = body.mass * SMatrix{3,3,T,9}(Diagonal(sones(T,3)))
    ∇q1 = -2 / Δt * LVᵀmat(q2)' * ∂qLVᵀmat(body.inertia * VLᵀmat(q1) * vector(q2))
    ∇q1 += -2 / Δt * LVᵀmat(q2)' * LVᵀmat(q1) * body.inertia * ∂qVLᵀmat(vector(q2))
    ∇ϕ15 = ∇q1 * rotational_integrator_jacobian_velocity(q2, -ϕ15, Δt)
    ∇15 = [∇v15 szeros(T,3,3);
           szeros(T,3,3) ∇ϕ15]

    # current configuration: z2 = x2, q2
    # manipulator's equation contribution
    ∇tra_x2 = - 2 / Δt * body.mass * SMatrix{3,3,T,9}(Diagonal(sones(T,3)))
    ∇tra_q2 = szeros(T,3,3)
    ∇rot_x2 = szeros(T,3,3)
    ∇rot_q2 = -2 / Δt * VLᵀmat(q2) * LVᵀmat(q1) * body.inertia * VLᵀmat(q1)
    ∇rot_q2 += -2 / Δt * VLᵀmat(q2) * Tmat() * RᵀVᵀmat(q3) * body.inertia * ∂qVLᵀmat(vector(q3))
    ∇rot_q2 += -2 / Δt * ∂qVLᵀmat(LVᵀmat(q1) * body.inertia * VLᵀmat(q1) * vector(q2) + Tmat() * RᵀVᵀmat(q3) * body.inertia * VLᵀmat(q2) * vector(q3))
    ∇rot_q2 *= LVᵀmat(q2)
    ∇z2 = [∇tra_x2 ∇tra_q2;
           ∇rot_x2 ∇rot_q2]
    # @show ∇z2
    # TODO
    # # constact constraints impulses contribution
    # @warn "000"
    return [∇m ∇J ∇15 0.0000000000*∇z2] #TODO not sure why we need to zero out this block, maybe finite diff is not correct and we try to match finite diff
end

function body_constraint_jacobian_body_data(mechanism::Mechanism, bodya::Node{T},
        bodyb::Node{T}, joint::JointConstraint{T,N,Nc}) where {T,N,Nc}
    # Jacobian of Bodya's dynamics constraints wrt Bodyb's data (x2b, q2b)
    # this comes from the fact that the Joint Constraint force mapping of Bodya
    # depends on Bodyb's data (x2b, q2b)
    # This is the same for spring and damper forces.
    timestep = mechanism.timestep

    ∇z2_aa = szeros(T,6,6)
    ∇z2_ab = szeros(T,6,6)
    # joint constraints impulses contribution
    for i = 1:Nc
        λ = getλJoint(joint, i)
        if bodyb.id == joint.child_id
            ∇z2_aa += impulse_map_parent_jacobian_parent([joint.translational, joint.rotational][i],
                bodya, bodyb, λ)
            ∇z2_ab += impulse_map_parent_jacobian_child([joint.translational, joint.rotational][i],
                bodya, bodyb, λ)
        elseif bodya.id == joint.child_id
            ∇z2_aa += impulse_map_child_jacobian_child([joint.translational, joint.rotational][i],
                bodyb, bodya, λ)
            ∇z2_ab += impulse_map_child_jacobian_parent([joint.translational, joint.rotational][i],
                bodyb, bodya, λ)
        end
    end
    # spring and damper impulses contribution
    if joint.spring
        for i = 1:Nc
            λ = getλJoint(joint, i)
            if bodyb.id == joint.child_id
                ∇z2_aa += spring_parent_jacobian_configuration_parent(
                    [joint.translational, joint.rotational][i], bodya, bodyb, timestep)
                ∇z2_ab += spring_parent_jacobian_configuration_child(
                    [joint.translational, joint.rotational][i], bodya, bodyb, timestep)
            elseif bodya.id == joint.child_id
                ∇z2_aa += spring_child_jacobian_configuration_child(
                    [joint.translational, joint.rotational][i], bodyb, bodya, timestep)
                ∇z2_ab += spring_child_jacobian_configuration_parent(
                    [joint.translational, joint.rotational][i], bodyb, bodya, timestep)
            end
        end
    end
    if joint.damper
        for i = 1:Nc
            λ = getλJoint(joint, i)
            if bodyb.id == joint.child_id
                ∇z2_aa += damper_parent_jacobian_configuration_parent(
                    [joint.translational, joint.rotational][i], bodya, bodyb, timestep)
                ∇z2_ab += damper_parent_jacobian_configuration_child(
                    [joint.translational, joint.rotational][i], bodya, bodyb, timestep)
            elseif bodya.id == joint.child_id
                ∇z2_aa += damper_child_jacobian_configuration_child(
                    [joint.translational, joint.rotational][i], bodyb, bodya, timestep)
                ∇z2_ab += damper_child_jacobian_configuration_parent(
                    [joint.translational, joint.rotational][i], bodyb, bodya, timestep)
            end
        end
    end

    # TODO
    # # contact constraints impulses contribution
    # ∇z2 +=
    # TODO
    return [szeros(T,6,13) ∇z2_aa], [szeros(T,6,13) ∇z2_ab]
end

function body_constraint_jacobian_body_data(mechanism::Mechanism, body::Node{T},
        contact::ContactConstraint{T,N,Nc}) where {T,N,Nc}
    # Jacobian of the Body's dynamics constraints wrt the Body's data (x2, q2)
    # this comes from the fact that the Contact Constraint force mapping depends
    # on the Body's data (x2, q2)
    # contact constraints impulses contribution
    ∇z3 = impulse_map_jacobian_configuration(mechanism, body, contact)
    ∇z2 = ∇z3 * integrator_jacobian_configuration(body, mechanism.timestep)
    return [szeros(T,6,13) ∇z2]
end

function body_constraint_jacobian_joint_data(mechanism::Mechanism{T}, body::Body{T},
    joint::JointConstraint{T}) where {T}
    Δt = mechanism.timestep
    Nd = data_dim(joint)
    N = 6
    x1, v15, q1, ϕ15 = previous_configuration_velocity(body.state)
    x2, v25, q2, ϕ25 = current_configuration_velocity(body.state)
    x3, q3 = next_configuration(body.state, Δt)
    # ∇u = Diagonal(SVector{6,T}(1,1,1,2,2,2)) * input_jacobian_control(mechanism, joint, body)
    ∇u = Diagonal(SVector{6,T}(1,1,1,1,1,1)) * input_jacobian_control(mechanism, joint, body)
    ∇spring = joint.spring ? apply_spring(mechanism, joint, body, unitary=true) : szeros(T,6,1)
    ∇damper = joint.damper ? apply_damper(mechanism, joint, body, unitary=true) : szeros(T,6,1)
    return [∇u ∇spring ∇damper]
end

function body_constraint_jacobian_contact_data(mechanism::Mechanism, body::Body{T},
    contact::ContactConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs<:NonlinearContact{T,N},N½}
    Nd = data_dim(contact)
    model = contact.model
    offset = model.offset
    x3, q3 = next_configuration(body.state, mechanism.timestep)
    γ = contact.impulses[2]

    ∇friction_coefficient = szeros(T,3,1)

    X = force_mapping(model, x3, q3)
    # this what we differentiate: Qᵀγ = - skew(p - vrotate(offset, inv(q3))) * VRmat(q3) * LᵀVᵀmat(q3) * X' * γ
    ∇p = - ∂pskew(VRmat(q3) * LᵀVᵀmat(q3) * X' * γ)
    ∇off = - ∂pskew(VRmat(q3) * LᵀVᵀmat(q3) * X' * γ) * -∂vrotate∂p(offset, inv(q3))

    ∇X = szeros(T,3,Nd)
    ∇Q = -[∇friction_coefficient ∇off ∇p]
    return [∇X; ∇Q]
end


function contact_constraint_jacobian_contact_data(mechanism::Mechanism, contact::ContactConstraint{T,N,Nc,Cs,N½},
        body::Body{T}) where {T,N,Nc,Cs<:NonlinearContact{T,N},N½}
    Nd = data_dim(contact)
    model = contact.model
    p = model.contact_point
    offset = model.offset
    x2, v25, q2, ϕ25 = current_configuration_velocity(body.state)
    x3, q3 = next_configuration(body.state, mechanism.timestep)
    s = contact.impulses_dual[2]
    γ = contact.impulses[2]

    ∇friction_coefficient = SA[0,γ[1],0,0]
    ∇off = [-model.surface_normal_projector; szeros(T,1,3); -model.surface_projector * skew(vrotate(ϕ25, q3))]
    ∇p = [model.surface_normal_projector * ∂vrotate∂p(model.contact_point, q3); szeros(T,1,3); model.surface_projector * skew(vrotate(ϕ25, q3)) * ∂vrotate∂p(model.contact_point, q3)]

    ∇compμ = szeros(T,N½,Nd)
    ∇g = -[∇friction_coefficient ∇off ∇p]
    return [∇compμ; ∇g]
end

function contact_constraint_jacobian_body_data(mechanism::Mechanism, contact::ContactConstraint{T,N,Nc,Cs,N½},
        body::Body{T}) where {T,N,Nc,Cs,N½}
    Nd = data_dim(body)
    ∇compμ = szeros(T,N½,Nd)
    ∇m = szeros(T,N½,1)
    ∇J = szeros(T,N½,6)
    ∇v15 = szeros(T,N½,3)
    ∇ϕ15 = szeros(T,N½,3)
    ∇z3 = - constraint_jacobian_configuration(mechanism, contact, body) # minus sign coming from res = [-compμ; -constraint]
    ∇z2 = ∇z3 * integrator_jacobian_configuration(body, mechanism.timestep, attjac=true) # 4x7 * 7x6 = 4x6
    ∇g = [∇m ∇J ∇v15 ∇ϕ15 ∇z2]
    return [∇compμ; ∇g]
end

################################################################################
# System Data Jacobians
################################################################################
function data_adjacency_matrix(joints::Vector{<:JointConstraint}, bodies::Vector{<:Body}, contacts::Vector{<:ContactConstraint})
    # mode can be impulses or data depending on whi
    nodes = [joints; bodies; contacts]
    n = length(nodes)
    A = zeros(Bool, n, n)

    for node1 in nodes
        for node2 in nodes
            T1 = typeof(node1)
            T2 = typeof(node2)
            if T1 <: Body
                if T2 <: Body
                    (node1.id == node2.id) && (A[node1.id, node2.id] = 1) # self loop
                    linked = length(indirect_link0(node1.id, node2.id, [joints; contacts])) > 0
                    linked && (A[node1.id, node2.id] = 1) # linked through a common joint
                elseif T2 <: JointConstraint
                    (node1.id == node2.parent_id || node1.id == node2.child_id) && (A[node1.id, node2.id] = 1) # linked
                elseif T2 <: ContactConstraint
                    (node1.id == node2.parent_id || node1.id == node2.child_id) && (A[node1.id, node2.id] = 1) # linked
                end
            elseif T1 <: JointConstraint
                if T2 <: Body
                    (node2.id == node1.parent_id || node2.id == node1.child_id) && (A[node1.id, node2.id] = 1) # linked
                end
            elseif T1 <: ContactConstraint
                if T2 <: Body
                    (node2.id == node1.parent_id || node2.id == node1.child_id) && (A[node1.id, node2.id] = 1) # linked
                elseif T2 <: ContactConstraint
                    (node1.id == node2.id) && (A[node1.id, node2.id] = 1) # self loop
                end
            end
        end
    end
    A = convert(Matrix{Int64}, A)
    return A
end

function indirect_link0(id1, id2, nodes::Vector{S}) where {S<:Node}
    ids = zeros(Int, 0)
    for node in nodes
        parent_id = node.parent_id
        (parent_id == nothing) && (parent_id = 0) #handle the origin's corner case
        linked = (id1 == node.child_id) && (id2 == parent_id)
        linked |= (id2 == node.child_id) && (id1 == parent_id)
        linked && push!(ids, node.id)
    end
    return ids
    # mech = get_halfcheetah()
    # @test indirect_link0(8,14,mech.joints) == [2]
    # @test indirect_link0(14,7,mech.joints) == []
    # @test indirect_link0(7,7,mech.joints) == []
end

function create_data_matrix(joints::Vector{<:JointConstraint}, bodies::Vector{B},
        contacts::Vector{<:ContactConstraint}; force_static::Bool=false) where {T,B<:Body{T}}
    nodes = [joints; bodies; contacts]
    A = data_adjacency_matrix(joints, bodies, contacts)
    dimrow = length.(nodes)
    dimcol = data_dim.(nodes)

    N = length(dimrow)
    static = force_static || (all(dimrow.<=10) && all(dimcol.<=10))
    data_matrix = spzeros(Entry,N,N)

    for i = 1:N
        for j = 1:N
            if A[i,j] == 1
                data_matrix[i,j] = Entry{T}(dimrow[i], dimcol[j], static = static)
            end
        end
    end
    return data_matrix
end

function jacobian_data!(data_matrix::SparseMatrixCSC, mechanism::Mechanism)
    jacobian_contact_data!(data_matrix, mechanism)
    jacobian_body_data!(data_matrix, mechanism)
    jacobian_joint_data!(data_matrix, mechanism)
    return nothing
end

function jacobian_contact_data!(data_matrix::SparseMatrixCSC, mechanism::Mechanism{T}) where {T}
    # ∂body∂ineqcdata
    for contact in mechanism.contacts
        pbody = get_body(mechanism, contact.parent_id)
        data_matrix[pbody.id, contact.id].value += body_constraint_jacobian_contact_data(mechanism, pbody, contact)
    end
    # ∂contact∂contactdata
    for contact in mechanism.contacts
        pbody = get_body(mechanism, contact.parent_id)
        data_matrix[contact.id, contact.id].value += contact_constraint_jacobian_contact_data(mechanism, contact, pbody)
    end
    return nothing
end

function jacobian_joint_data!(data_matrix::SparseMatrixCSC, mechanism::Mechanism{T}) where T
    # ∂body∂jointdata
    # TODO adapt this to handle cycles
    for body in mechanism.bodies
        for joint in mechanism.joints
            if (body.id == joint.parent_id) || (body.id == joint.child_id)
                data_matrix[body.id, joint.id].value += body_constraint_jacobian_joint_data(mechanism, body, joint)
            end
        end
    end
    return nothing
end

function jacobian_body_data!(data_matrix::SparseMatrixCSC, mechanism::Mechanism{T}) where T
    # ∂joint∂bodydata
    # TODO adapt this to handle cycles
    for body in mechanism.bodies
        for joint in mechanism.joints
            if (body.id == joint.parent_id) || (body.id == joint.child_id)
                data_matrix[joint.id, body.id].value += joint_constraint_jacobian_body_data(mechanism, joint, body)
            end
        end
    end
    # ∂body∂bodydata
    for body1 in mechanism.bodies
        data_matrix[body1.id, body1.id].value += body_constraint_jacobian_body_data(mechanism, body1)

        for body2 in [mechanism.bodies; mechanism.origin]
            joint_links = indirect_link0(body1.id, body2.id, mechanism.joints)
            joints = [get_joint_constraint(mechanism, id) for id in joint_links]
            for joint in joints
                ∇11, ∇12 = body_constraint_jacobian_body_data(mechanism, body1, body2, joint)
                (typeof(body1) <: Body) && (data_matrix[body1.id, body1.id].value += ∇11)
                (typeof(body1) <: Body && typeof(body2) <: Body) && (data_matrix[body1.id, body2.id].value += ∇12)
            end
            # pretty sure this is useless because contact is never linked to two bodies
            # contact_links = indirect_link0(body1.id, body2.id, mechanism.contacts)
            # contacts = [get_joint_constraint(mechanism, id) for id in contact_links]
            # for contact in contacts
                # data_matrix[body1.id, body2.id].value += body_constraint_jacobian_body_data(mechanism, body1, body2, contact)
            # end
        end
    end
    for contact in mechanism.contacts
        body = get_body(mechanism, contact.parent_id)
        data_matrix[body.id, body.id].value += body_constraint_jacobian_body_data(mechanism, body, contact)
    end
    # ∂contact∂bodydata
    for contact in mechanism.contacts
        pbody = get_body(mechanism, contact.parent_id)
        data_matrix[contact.id, pbody.id].value += contact_constraint_jacobian_body_data(mechanism, contact, pbody)
    end
    return nothing
end
