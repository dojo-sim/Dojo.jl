################################################################################
# Data Jacobians
################################################################################
function joint_constraint_jacobian_body_data(mechanism::Mechanism, joint::JointConstraint{T,N}, body::Body{T}) where {T,N}
    Nd = data_dim(body)
    ∇m = szeros(T,N,1)
    ∇J = szeros(T,N,6)
    ∇v15 = szeros(T,N,3)
    ∇ω15 = szeros(T,N,3)
    ∇z2 = -constraint_jacobian_configuration(mechanism, joint, body) *
        integrator_jacobian_configuration(body, mechanism.timestep, attjac=true)
    ∇g = [∇m ∇J ∇v15 ∇ω15 ∇z2]
    return ∇g
end

function body_constraint_jacobian_body_data(mechanism::Mechanism, body::Body{T}) where T
    Δt = mechanism.timestep
    N = 6
    x1, v15, q1, ω15 = previous_configuration_velocity(body.state)
    x2, v25, q2, ω25 = current_configuration_velocity(body.state)
    x3, q3 = next_configuration(body.state, Δt)
    # Mass
    gravity=mechanism.gravity
    ∇m = [1 / Δt * (x2 - x1) + Δt/2 * gravity - 1 / Δt * (x3 - x2) + Δt/2 * gravity;
          szeros(T,3,1)]
    # Inertia
    ∇J = 2 / Δt * LVᵀmat(q2)' * LVᵀmat(q1) * ∂Jp∂J(VLᵀmat(q1) * vector(q2))
    ∇J += 2 / Δt * LVᵀmat(q2)' * Tmat() * RᵀVᵀmat(q3) * ∂Jp∂J(VLᵀmat(q2) * vector(q3))
    ∇J = [szeros(T,3,6); ∇J]

    # initial conditions: v15, ω15
    ∇v15 = body.mass * SMatrix{3,3,T,9}(Diagonal(sones(T,3)))
    ∇q1 = -2 / Δt * LVᵀmat(q2)' * ∂LVᵀmat∂q(body.inertia * VLᵀmat(q1) * vector(q2))
    ∇q1 += -2 / Δt * LVᵀmat(q2)' * LVᵀmat(q1) * body.inertia * ∂VLᵀmat∂q(vector(q2))
    ∇ω15 = ∇q1 * rotational_integrator_jacobian_velocity(q2, -ω15, Δt)
    ∇15 = [∇v15 szeros(T,3,3);
           szeros(T,3,3) ∇ω15]

    # current configuration: z2 = x2, q2
    # manipulator's equation contribution
    ∇tra_x2 = - 2 / Δt * body.mass * SMatrix{3,3,T,9}(Diagonal(sones(T,3)))
    ∇tra_q2 = szeros(T,3,3)
    ∇rot_x2 = szeros(T,3,3)
    ∇rot_q2 = -2 / Δt * VLᵀmat(q2) * LVᵀmat(q1) * body.inertia * VLᵀmat(q1)
    ∇rot_q2 += -2 / Δt * VLᵀmat(q2) * Tmat() * RᵀVᵀmat(q3) * body.inertia * ∂VLᵀmat∂q(vector(q3))
    ∇rot_q2 += -2 / Δt * ∂VLᵀmat∂q(LVᵀmat(q1) * body.inertia * VLᵀmat(q1) * vector(q2) + Tmat() * RᵀVᵀmat(q3) * body.inertia * VLᵀmat(q2) * vector(q3))
    ∇rot_q2 *= LVᵀmat(q2)
    ∇z2 = [∇tra_x2 ∇tra_q2;
           ∇rot_x2 ∇rot_q2]
    # @show ∇z2
    # TODO
    # # contact constraints impulses contribution
    return [∇m ∇J ∇15 0.0000000000*∇z2] #TODO not sure why we need to zero out this block, maybe finite diff is not correct and we try to match finite diff
end

function body_constraint_jacobian_body_data(mechanism::Mechanism, pbody::Node{T},
        cbody::Node{T}, joint::JointConstraint{T,N,Nc}) where {T,N,Nc}
    # Jacobian of pbody's dynamics constraints wrt cbody's data (x2b, q2b)
    # this comes from the fact that the Joint Constraint force mapping of pbody
    # depends on cbody's data (x2b, q2b)
    # This is the same for spring and damper forces.
    timestep= mechanism.timestep

    ∇z2_aa = szeros(T,6,6)
    ∇z2_ab = szeros(T,6,6)
    # joint constraints impulses contribution
    for i = 1:Nc
        λ = get_joint_impulses(joint, i)
        if cbody.id == joint.child_id
            ∇z2_aa += impulse_map_jacobian(:parent, :parent, (joint.translational, joint.rotational)[i],
                pbody, cbody, λ)
            ∇z2_ab += impulse_map_jacobian(:parent, :child, (joint.translational, joint.rotational)[i],
                pbody, cbody, λ)
        elseif pbody.id == joint.child_id
            ∇z2_aa += impulse_map_jacobian(:child, :child, (joint.translational, joint.rotational)[i],
                cbody, pbody, λ)
            ∇z2_ab += impulse_map_jacobian(:child, :parent, (joint.translational, joint.rotational)[i],
                cbody, pbody, λ)
        end
    end
    # spring and damper impulses contribution
    if joint.spring
        for i = 1:Nc
            λ = get_joint_impulses(joint, i)
            if cbody.id == joint.child_id
                ∇z2_aa += spring_jacobian_configuration(
                    :parent, :parent,
                    (joint.translational, joint.rotational)[i], pbody, cbody, timestep)
                ∇z2_ab += spring_jacobian_configuration(
                    :parent, :child,
                    (joint.translational, joint.rotational)[i], pbody, cbody, timestep)
            elseif pbody.id == joint.child_id
                ∇z2_aa += spring_jacobian_configuration(
                    :child, :child,
                    (joint.translational, joint.rotational)[i], cbody, pbody, timestep)
                ∇z2_ab += spring_jacobian_configuration(
                    :child, :parent,
                    (joint.translational, joint.rotational)[i], cbody, pbody, timestep)
            end
        end
    end
    if joint.damper
        for i = 1:Nc
            λ = get_joint_impulses(joint, i)
            if cbody.id == joint.child_id
                ∇z2_aa += damper_jacobian_configuration(
                    :parent, :parent,
                    (joint.translational, joint.rotational)[i], pbody, cbody, timestep)
                ∇z2_ab += damper_jacobian_configuration(
                    :parent, :child,
                    (joint.translational, joint.rotational)[i], pbody, cbody, timestep)
            elseif pbody.id == joint.child_id
                ∇z2_aa += damper_jacobian_configuration(
                    :child, :child,
                    (joint.translational, joint.rotational)[i], cbody, pbody, timestep)
                ∇z2_ab += damper_jacobian_configuration(
                    :child, :parent,
                    (joint.translational, joint.rotational)[i], cbody, pbody, timestep)
            end
        end
    end
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
    x1, v15, q1, ω15 = previous_configuration_velocity(body.state)
    x2, v25, q2, ω25 = current_configuration_velocity(body.state)
    x3, q3 = next_configuration(body.state, Δt)
    # ∇u = Diagonal(SVector{6,T}(1,1,1,2,2,2)) * input_jacobian_control(mechanism, joint, body)
    ∇u = Diagonal(SVector{6,T}(1,1,1,1,1,1)) * input_jacobian_control(mechanism, joint, body)
    ∇spring = joint.spring ? spring_impulses(mechanism, joint, body, unitary=true) : szeros(T,6,1)
    ∇damper = joint.damper ? damper_impulses(mechanism, joint, body, unitary=true) : szeros(T,6,1)
    return [∇u ∇spring ∇damper]
end

function body_constraint_jacobian_contact_data(mechanism::Mechanism, body::Body{T},
        contact::RigidContactConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs<:NonlinearContact{T,N},N½}
    Nd = data_dim(contact)
    model = contact.model
    contact_radius = model.collision.contact_radius
    offset = model.collision.contact_normal' * contact_radius
    xp3, qp3 = next_configuration(get_body(mechanism, contact.parent_id).state, mechanism.timestep)
    xc3, qc3 = next_configuration(get_body(mechanism, contact.child_id).state, mechanism.timestep)

    γ = contact.impulses[2]

    ∇friction_coefficient = szeros(T,3,1)

    X = force_mapping(:parent, model, xp3, qp3, xc3, qc3)
    ∇p = - ∂skew∂p(VRmat(qp3) * LᵀVᵀmat(qp3) * X * γ)
    ∇contact_radius = - ∂skew∂p(VRmat(qp3) * LᵀVᵀmat(qp3) * X * γ) * -rotation_matrix(inv(qp3)) * model.collision.contact_normal'

    ∇X = szeros(T,3,Nd)
    ∇Q = -[∇friction_coefficient ∇contact_radius ∇p]
    return [∇X; ∇Q]
end

function contact_constraint_jacobian_contact_data(mechanism::Mechanism, contact::RigidContactConstraint{T,N,Nc,Cs,N½}, body::Body{T}) where {T,N,Nc,Cs<:NonlinearContact{T,N},N½}
    Nd = data_dim(contact)
    model = contact.model

    xp3, vp25, qp3, ωp25 = next_configuration_velocity(get_body(mechanism, contact.parent_id).state, mechanism.timestep)
    xc3, vc25, qc3, ωc25 = next_configuration_velocity(get_body(mechanism, contact.child_id).state, mechanism.timestep)

    s = contact.impulses_dual[2]
    γ = contact.impulses[2]

    ∇friction_coefficient = SA[0,γ[1],0,0]
    ∇contact_radius = [-model.collision.contact_normal; szeros(T,1,3); -model.collision.contact_tangent * skew(vector_rotate(ωp25, qp3))] * model.collision.contact_normal'
    ∇p = [model.collision.contact_normal * rotation_matrix(qp3); szeros(T,1,3); model.collision.contact_tangent * skew(vector_rotate(ωp25, qp3)) * rotation_matrix(qp3)]

    ∇compμ = szeros(T,N½,Nd)
    ∇g = -[∇friction_coefficient ∇contact_radius ∇p]

    return [∇compμ; ∇g]
end

function contact_constraint_jacobian_body_data(mechanism::Mechanism, contact::RigidContactConstraint{T,N,Nc,Cs,N½}, body::Body{T}) where {T,N,Nc,Cs,N½}
    Nd = data_dim(body)
    ∇compμ = szeros(T,N½,Nd)
    ∇m = szeros(T,N½,1)
    ∇J = szeros(T,N½,6)
    ∇v15 = szeros(T,N½,3)
    ∇ω15 = szeros(T,N½,3)
    ∇z3 = - constraint_jacobian_configuration(mechanism, contact, body) # minus sign coming from res = [-compμ; -constraint]
    ∇z2 = ∇z3 * integrator_jacobian_configuration(body, mechanism.timestep, attjac=true) # 4x7 * 7x6 = 4x6
    ∇g = [∇m ∇J ∇v15 ∇ω15 ∇z2]
    return [∇compμ; ∇g]
end

################################################################################
# System Data Jacobians
################################################################################
function data_adjacency_matrix(joints::Vector{<:JointConstraint}, bodies::Vector{<:Body},
        contacts::Vector{<:ContactConstraint})
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
                    linked = length(indirect_link(node1.id, node2.id, [joints; contacts])) > 0
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

function indirect_link(id1, id2, nodes::Vector{S}) where {S<:Node}
    ids = zeros(Int, 0)
    for node in nodes
        parent_id = node.parent_id
        (parent_id == nothing) && (parent_id = 0) #handle the origin's corner case
        linked = (id1 == node.child_id) && (id2 == parent_id)
        linked |= (id2 == node.child_id) && (id1 == parent_id)
        linked && push!(ids, node.id)
    end
    return ids
end

function create_data_matrix(joints::Vector{<:JointConstraint}, bodies::Vector{B}, contacts::Vector{<:ContactConstraint};
        force_static::Bool=false) where {T,B<:Body{T}}
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
    for pbody in mechanism.bodies
        data_matrix[pbody.id, pbody.id].value += body_constraint_jacobian_body_data(mechanism, pbody)

        for cbody in [mechanism.bodies; mechanism.origin]
            joint_links = indirect_link(pbody.id, cbody.id, mechanism.joints)
            joints = [get_joint(mechanism, id) for id in joint_links]
            for joint in joints
                ∇11, ∇12 = body_constraint_jacobian_body_data(mechanism, pbody, cbody, joint)
                (typeof(pbody) <: Body) && (data_matrix[pbody.id, pbody.id].value += ∇11)
                (typeof(pbody) <: Body && typeof(cbody) <: Body) && (data_matrix[pbody.id, cbody.id].value += ∇12)
            end
            # pretty sure this is useless because contact is never linked to two bodies
            # contact_links = indirect_link(pbody.id, cbody.id, mechanism.contacts)
            # contacts = [get_joint(mechanism, id) for id in contact_links]
            # for contact in contacts
                # data_matrix[pbody.id, cbody.id].value += body_constraint_jacobian_body_data(mechanism, pbody, cbody, contact)
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
