function joint_constraint_jacobian(mechanism::Mechanism{T,Nn,Ne,Nb}) where {T,Nn,Ne,Nb}
    timestep = mechanism.timestep
    joints = mechanism.joints
    nc = sum(length.(joints))
    Gl = zeros(T,nc,13Nb)

    oneindc = 0
    for joint in joints
        ind1 = 1
        ind2 = 0

        parentid = joint.parentid
        if parentid !== nothing
            parentind = parentid - Ne
            pbody = get_body(mechanism,parentid)
            pstate = pbody.state

            for (i,childid) in enumerate(joint.childids)
                childind = childid - Ne
                cbody = get_body(mechanism,childid)
                cstate = cbody.state
                element = joint.constraints[i]

                ind2 += length(element)
                range = oneindc+ind1:oneindc+ind2

                pcol13 = offset_range(parentind,13)
                ccol13 = offset_range(childind,13)

                p = constraint_jacobian_parent(element, pbody, cbody, joint.λsol[2], timestep) # x3
                c = constraint_jacobian_child(element, pbody, cbody, joint.λsol[2], timestep) # x3

                pGlx = p[:, 1:3]
                pGlq = p[:, 4:7]
                cGlx = c[:, 1:3]
                cGlq = c[:, 4:7]
                @show size(p)
                @show size(c)
                @show element
                @show length(element)
                @show size(pGlx)
                @show range
                @show pcol13[1:3]
                Gl[range, pcol13[1:3]] = pGlx
                Gl[range, pcol13[7:10]] = pGlq
                Gl[range,ccol13[1:3]] = cGlx
                Gl[range,ccol13[7:10]] = cGlq
                ind1 = ind2+1
            end
        else
            for (i,childid) in enumerate(joint.childids)
                childind = childid - Ne
                cbody = get_body(mechanism,childid)
                cstate = cbody.state
                element = joint.constraints[i]
                ind2 += length(element)
                range = oneindc+ind1:oneindc+ind2

                ccol13 = offset_range(childind,13)

                c = constraint_jacobian_child(element, mechanism.origin, cbody, joint.λsol[2], timestep) # x3

                Gl[range,ccol13[1:3]] = c[:, 1:3]
                Gl[range,ccol13[7:10]] = c[:, 4:7]
                ind1 = ind2+1
            end
        end
        oneindc += length(joint)
    end
    return Gl
end

function joint_dynamics_jacobian(mechanism::Mechanism{T,Nn,Ne,Nb}) where {T,Nn,Ne,Nb}
    timestep = mechanism.timestep
    J = zeros(T,6Nb,13Nb)

    for joint in mechanism.joints
        parentid = joint.parentid

        if parentid !== nothing
            parentind = parentid - Ne
            pbody = get_body(mechanism, parentid)
            pstate = pbody.state
            prow6 = offset_range(parentind,6)
            pcol13 = offset_range(parentind,13)

            for (i, childid) in enumerate(joint.childids)
                childind = childid - Ne
                cbody = get_body(mechanism, childid)
                cstate = cbody.state
                element = joint.constraints[i]
                crow6 = offset_range(childind,6)
                ccol13 = offset_range(childind,13)
                λ = getλJoint(joint, i)

                Aaa = zeros(T,6,13)
                Aab = zeros(T,6,13)
                Aba = zeros(T,6,13)
                Abb = zeros(T,6,13)

                xa, qa = current_configuration(pstate)
                xb, qb = current_configuration(cstate)

                if typeof(impulse_map_parent(element, xa, qa, xb, qb, joint.λsol[2])) <: AbstractArray && length(element) > 0
                    XX = FiniteDiff.finite_difference_jacobian(w -> -impulse_map_parent(element, w, qa, xb, qb, joint.λsol[2])[1:3, :] * λ, xa)
                    XQ = FiniteDiff.finite_difference_jacobian(w -> -impulse_map_parent(element, xa, UnitQuaternion(w..., false), xb, qb, joint.λsol[2])[1:3, :] * λ, [qa.w; qa.x; qa.y; qa.z])
                    QX = FiniteDiff.finite_difference_jacobian(w -> -impulse_map_parent(element, w, qa, xb, qb, joint.λsol[2])[4:6,:] * λ, xa)
                    QQ = FiniteDiff.finite_difference_jacobian(w -> -impulse_map_parent(element, xa, UnitQuaternion(w..., false), xb, qb, joint.λsol[2])[4:6, :] * λ, [qa.w; qa.x; qa.y; qa.z])

                    Aaa[1:3,1:3] = XX
                    Aaa[1:3,7:10] = XQ
                    Aaa[4:6,1:3] = QX
                    Aaa[4:6,7:10] = QQ
                end

                if typeof(impulse_map_parent(element, xa, qa, xb, qb, joint.λsol[2])) <: AbstractArray && length(element) > 0
                    XX = FiniteDiff.finite_difference_jacobian(w -> -impulse_map_parent(element, xa, qa, w, qb, joint.λsol[2])[1:3, :] * λ, xb)
                    XQ = FiniteDiff.finite_difference_jacobian(w -> -impulse_map_parent(element, xa, qa, xb, UnitQuaternion(w..., false), joint.λsol[2])[1:3, :] * λ, [qb.w; qb.x; qb.y; qb.z])
                    QX = FiniteDiff.finite_difference_jacobian(w -> -impulse_map_parent(element, xa, qa, w, qb, joint.λsol[2])[4:6, :] * λ, xb)
                    QQ = FiniteDiff.finite_difference_jacobian(w -> -impulse_map_parent(element, xa, qa, xb, UnitQuaternion(w..., false), joint.λsol[2])[4:6, :] * λ, [qb.w; qb.x; qb.y; qb.z])

                    Aab[1:3,1:3] = XX
                    Aab[1:3,7:10] = XQ
                    Aab[4:6,1:3] = QX
                    Aab[4:6,7:10] = QQ
                end

                if typeof(impulse_map_child(element, xa, qa, xb, qb, joint.λsol[2])) <: AbstractArray && length(element) > 0
                    XX = FiniteDiff.finite_difference_jacobian(w -> -impulse_map_child(element, w, qa, xb, qb, joint.λsol[2])[1:3, :] * λ, xa)
                    XQ = FiniteDiff.finite_difference_jacobian(w -> -impulse_map_child(element, xa, UnitQuaternion(w..., false), xb, qb, joint.λsol[2])[1:3, :] * λ, [qa.w; qa.x; qa.y; qa.z])
                    QX = FiniteDiff.finite_difference_jacobian(w -> -impulse_map_child(element, w, qa, xb, qb, joint.λsol[2])[4:6, :] * λ, xa)
                    QQ = FiniteDiff.finite_difference_jacobian(w -> -impulse_map_child(element, xa, UnitQuaternion(w..., false), xb, qb, joint.λsol[2])[4:6, :] * λ, [qa.w; qa.x; qa.y; qa.z])

                    Aba[1:3,1:3] = XX
                    Aba[1:3,7:10] = XQ
                    Aba[4:6,1:3] = QX
                    Aba[4:6,7:10] = QQ
                end

                if typeof(impulse_map_child(element, xa, qa, xb, qb, joint.λsol[2])) <: AbstractArray && length(element) > 0
                    XX = FiniteDiff.finite_difference_jacobian(w -> -impulse_map_child(element, xa, qa, w, qb, joint.λsol[2])[1:3, :] * λ, xb)
                    XQ = FiniteDiff.finite_difference_jacobian(w -> -impulse_map_child(element, xa, qa, xb, UnitQuaternion(w..., false), joint.λsol[2])[1:3, :] * λ, [qb.w; qb.x; qb.y; qb.z])
                    QX = FiniteDiff.finite_difference_jacobian(w -> -impulse_map_child(element, xa, qa, w, qb, joint.λsol[2])[4:6, :] * λ, xb)
                    QQ = FiniteDiff.finite_difference_jacobian(w -> -impulse_map_child(element, xa, qa, xb, UnitQuaternion(w..., false), joint.λsol[2])[4:6, :] * λ, [qb.w; qb.x; qb.y; qb.z])

                    Abb[1:3,1:3] = XX
                    Abb[1:3,7:10] = XQ
                    Abb[4:6,1:3] = QX
                    Abb[4:6,7:10] = QQ
                end

                J[prow6,pcol13] += Aaa
                J[prow6,ccol13] += Aab
                J[crow6,pcol13] += Aba
                J[crow6,ccol13] += Abb
            end
        else
            for (i, childid) in enumerate(joint.childids)
                childind = childid - Ne
                cbody = get_body(mechanism, childid)
                cstate = cbody.state
                element = joint.constraints[i]
                crow6 = offset_range(childind,6)
                ccol13 = offset_range(childind,13)
                λ = getλJoint(joint, i)

                Abb = zeros(T,6,13)

                xa, qa = current_configuration(mechanism.origin.state)
                xb, qb = current_configuration(cstate)

                if typeof(impulse_map_child(element, xa, qb, xb, qb, joint.λsol[2])) <: AbstractArray && length(element) > 0

                    XX = FiniteDiff.finite_difference_jacobian(w -> -impulse_map_child(element, xa, qa, w, qb, joint.λsol[2])[1:3, :] * λ, xb)
                    XQ = FiniteDiff.finite_difference_jacobian(w -> -impulse_map_child(element, xa, qa, xb, UnitQuaternion(w..., false), joint.λsol[2])[1:3, :] * λ, [qb.w; qb.x; qb.y; qb.z])
                    QX = FiniteDiff.finite_difference_jacobian(w -> -impulse_map_child(element, xa, qa, w, qb, joint.λsol[2])[4:6, :] * λ, xb)
                    QQ = FiniteDiff.finite_difference_jacobian(w -> -impulse_map_child(element, xa, qa, xb, UnitQuaternion(w..., false), joint.λsol[2])[4:6, :] * λ, [qb.w; qb.x; qb.y; qb.z])

                    Abb[1:3,1:3] = XX
                    Abb[1:3,7:10] = XQ
                    Abb[4:6,1:3] = QX
                    Abb[4:6,7:10] = QQ
                end

                J[crow6,ccol13] += Abb
            end
        end
    end
    return J
end

function springapply_damperjacobian(mechanism::Mechanism{T,Nn,Ne,Nb}) where {T,Nn,Ne,Nb}
    timestep = mechanism.timestep
    J = zeros(T,6Nb,13Nb)

    for joint in mechanism.joints
        parentid = joint.parentid

        if parentid !== nothing
            parentind = parentid - Ne
            pbody = get_body(mechanism, parentid)
            pstate = pbody.state
            prow6 = offset_range(parentind,6)
            pcol13 = offset_range(parentind,13)

            for (i, childid) in enumerate(joint.childids)
                childind = childid - Ne
                cbody = get_body(mechanism, childid)
                cstate = cbody.state
                element = joint.constraints[i]
                crow6 = offset_range(childind,6)
                ccol13 = offset_range(childind,13)
                λ = getλJoint(joint, i)

                Aaa = zeros(T,6,13)
                Aab = zeros(T,6,13)
                Aba = zeros(T,6,13)
                Abb = zeros(T,6,13)

                Aaa[:, [1:3; 7:10]] -= spring_parent_jacobian_configuration_parent(element, pbody, cbody, timestep, attjac = false)
                Aaa[:, [1:3; 7:10]] -= damper_parent_jacobian_configuration_parent(element, pbody, cbody, timestep, attjac = false)

                Aab[:, [1:3; 7:10]] -= spring_parent_jacobian_configuration_child(element, pbody, cbody, timestep, attjac = false)
                Aab[:, [1:3; 7:10]] -= damper_parent_jacobian_configuration_child(element, pbody, cbody, timestep, attjac = false)

                Aba[:, [1:3; 7:10]] -= spring_child_jacobian_configuraion_parent(element, pbody, cbody, timestep, attjac = false)
                Aba[:, [1:3; 7:10]] -= damper_child_jacobian_configuration_parent(element, pbody, cbody, timestep, attjac = false)

                Abb[:, [1:3; 7:10]] -= spring_child_jacobian_configuration_child(element, pbody, cbody, timestep, attjac = false)
                Abb[:, [1:3; 7:10]] -= damper_child_jacobian_configuration_child(element, pbody, cbody, timestep, attjac = false)

                J[prow6,pcol13] += Aaa
                J[prow6,ccol13] += Aab
                J[crow6,pcol13] += Aba
                J[crow6,ccol13] += Abb
            end
        else
            for (i, childid) in enumerate(joint.childids)
                pbody = mechanism.origin
                childind = childid - Ne
                cbody = get_body(mechanism, childid)
                cstate = cbody.state
                element = joint.constraints[i]
                crow6 = offset_range(childind,6)
                ccol13 = offset_range(childind,13)
                λ = getλJoint(joint, i)

                Abb = zeros(T,6,13)

                Abb[:, [1:3; 7:10]] -= spring_child_jacobian_configuration_child(element, pbody, cbody, timestep, attjac = false)
                Abb[:, [1:3; 7:10]] -= damper_child_jacobian_configuration_child(element, pbody, cbody, timestep, attjac = false)

                J[crow6,ccol13] += Abb
            end
        end
    end
    return J
end

function ∂body∂z(body::Body{T}, timestep::T; attjac::Bool = true) where T
    state = body.state
    q2 = state.q2[1]
    ϕ25 = state.ϕsol[2]
    Z3 = szeros(T,3,3)
    Z34 = szeros(T,3,4)
    ZT = attjac ? szeros(T,6,6) : szeros(T,6,7)
    ZR = szeros(T,7,6)

    x1, q1 = previous_configuration(state)
    x2, q2 = current_configuration(state)
    x3, q3 = next_configuration(state, timestep)

    AposT = [-I Z3]
    AvelT = [Z3 -I*body.m] # solving for impulses

    AposR = [-∂integrator∂q(q2, ϕ25, timestep, attjac = attjac) szeros(4,3)]

    rot_q1(q) = -4 / timestep * LVᵀmat(q2)' * Lmat(UnitQuaternion(q..., false)) * Vᵀmat() * body.J * Vmat() * Lmat(UnitQuaternion(q..., false))' * vector(q2)
    rot_q2(q) = -4 / timestep * LVᵀmat(UnitQuaternion(q..., false))' * Tmat() * Rmat(next_orientation(UnitQuaternion(q..., false), state.ϕsol[2], timestep))' * Vᵀmat() * body.J * Vmat() * Lmat(UnitQuaternion(q..., false))' * vector(next_orientation(UnitQuaternion(q..., false), state.ϕsol[2], timestep)) + -4 / timestep * LVᵀmat(UnitQuaternion(q..., false))' * Lmat(next_orientation(UnitQuaternion(q..., false), -state.ϕ15, timestep)) * Vᵀmat() * body.J * Vmat() * Lmat(next_orientation(UnitQuaternion(q..., false), -state.ϕ15, timestep))' * q

    dynR_ϕ15 = -1.0 * FiniteDiff.finite_difference_jacobian(rot_q1, vector(q1)) * ∂integrator∂ϕ(q2, -state.ϕ15, timestep)
    dynR_q2 = FiniteDiff.finite_difference_jacobian(rot_q2, vector(q2))
    AvelR = attjac ? [dynR_q2 * LVᵀmat(q2) dynR_ϕ15] : [dynR_q2 dynR_ϕ15]

    return [[AposT;AvelT] ZT;
             ZR [AposR;AvelR]]
end

function ∂body∂u(body::Body{T}, timestep) where T
    Z3 = szeros(T,3,3)
    Z43 = szeros(T,4,3)

    BposT = [Z3 Z3]
    BvelT = [-I Z3]
    BposR = [Z43 Z43]
    BvelR = [Z3 -2.0 * I]
    return [BposT;BvelT;BposR;BvelR]
end

function dynamics_jacobian(mechanism::Mechanism{T,Nn,Ne,Nb}, jointids) where {T,Nn,Ne,Nb}
    timestep = mechanism.timestep
    nu = control_dimension(mechanism)

    # get state linearization
    Fz = zeros(T,6Nb,13Nb)
    Fu = zeros(T,6Nb,6Nb)

    Bcontrol = zeros(T,6Nb,nu)

    for (i, body) in enumerate(mechanism.bodies)
        col6 = offset_range(i,6)
        col13 = offset_range(i,13)

        Fzi = ∂body∂z(body, timestep, attjac = false)[[4:6; 11:13],:]
        Fui = ∂body∂u(body, timestep)[[4:6; 11:13],:]

        Fz[col6,col13] = Fzi
        Fu[col6,col6] = Fui
    end

    n1 = 1
    n2 = 0

    for id in jointids
        joint = get_joint_constraint(mechanism, id)
        n2 += control_dimension(joint)

        parentid = joint.parentid
        if parentid != nothing
            parentind = parentid - Ne
            col6 = offset_range(parentind,6)
            Bcontrol[col6,n1:n2] = input_jacobian_control_parent(mechanism, joint, get_body(mechanism, parentid))
        end
        for childid in joint.childids
            childind = childid - Ne
            col6 = offset_range(childind,6)
            Bcontrol[col6,n1:n2] = input_jacobian_control_child(mechanism, joint, get_body(mechanism, childid))
        end

        n1 = n2+1
    end
    return Fz, Fu * Bcontrol
end

function contact_dynamics_jacobian(mechanism::Mechanism{T,Nn,Ne,Nb}) where {T,Nn,Ne,Nb}
    timestep = mechanism.timestep
    J = zeros(6Nb, 13Nb)

    offr = 0
    offc = 0
    for body in mechanism.bodies
        for contact in mechanism.contacts
            if contact.parentid == body.id
                bound = contact.constraints[1]
                bound_type = typeof(bound)
                x3, q3 = next_configuration(body.state, timestep)

                function d(vars)
                    x = vars[1:3]
                    q = UnitQuaternion(vars[4:7]..., false)
                    return impulse_map(bound, x, q, nothing)' * contact.γsol[2]
                end

                if bound_type <: NonlinearContact
                    # J[offr .+ (1:6), offc .+ [1:3; 7:10]] -= _dN(x3, vector(q3), contact.γsol[2][1:1], bound.p)
                    # J[offr .+ (1:6), offc .+ [1:3; 7:10]] -= _dB(x3, vector(q3), contact.γsol[2][2:4], bound.p)
                    J[offr .+ (1:6), offc .+ [1:3; 7:10]] -= FiniteDiff.finite_difference_jacobian(d, [x3; vector(q3)])
                elseif bound_type <: LinearContact
                    J[offr .+ (1:6), offc .+ [1:3; 7:10]] -= FiniteDiff.finite_difference_jacobian(d, [x3; vector(q3)])
                elseif bound_type <: ImpactContact
                    J[offr .+ (1:6), offc .+ [1:3; 7:10]] -= FiniteDiff.finite_difference_jacobian(d, [x3; vector(q3)])
                end
            end
        end
        offr += 6
        offc += 13
    end
    return J
end

function contact_constraint_jacobian(mechanism::Mechanism{T,Nn,Ne,Nb}) where {T,Nn,Ne,Nb}
    timestep = mechanism.timestep
    contacts = mechanism.contacts
    ncontacts = contact_dimension(mechanism)
    J = zeros(ncontacts, 13Nb)

    offr = 0
    for contact in contacts
        bound = contact.constraints[1]
        body = get_body(mechanism, contact.parentid)
        N½ = Int(length(contact)/2)
        x2, v25, q2, ϕ25 = current_configuration_velocity(body.state)
        x3, q3 = next_configuration(body.state, timestep)
        ibody = findfirst(x -> x.id == contact.parentid, mechanism.bodies)
        bound_type = typeof(bound)

        function d(vars)
            q3 = UnitQuaternion(vars..., false)
            # Bxmat = bound.Bx
            # Bqmat = Bxmat * ∂vrotate∂q(bound.p, q) * LVᵀmat(q)
            # return Bqmat * ϕ25
            vp = v25 + skew(vrotate(ϕ25, q3)) * (vrotate(bound.p, q3) - bound.offset)
            return bound.Bx * vp
        end
        if bound_type <: NonlinearContact
            J[offr + N½ .+ (1:1), (ibody-1)*13 .+ (1:3)] =  bound.ainv3
            J[offr + N½ .+ (1:1), (ibody-1)*13 .+ (7:10)] = bound.ainv3 * (VLmat(q3) * Lmat(UnitQuaternion(bound.p)) * Tmat() + VRᵀmat(q3) * Rmat(UnitQuaternion(bound.p)))#) * Rmat(quaternion_map(ϕ25, timestep)*timestep/2)*LVᵀmat(q2)
            J[offr + N½ .+ (3:4), (ibody-1)*13 .+ (7:10)] = FiniteDiff.finite_difference_jacobian(d, vector(q3))#dBω(vector(q3), ϕ25, bound.p)
        elseif bound_type <: LinearContact
            J[offr + N½ .+ (1:1), (ibody-1)*13 .+ (1:3)] =  bound.ainv3
            J[offr + N½ .+ (1:1), (ibody-1)*13 .+ (7:10)] = bound.ainv3 * ∂vrotate∂q(bound.p, q3)
            J[offr + N½ .+ (3:6), (ibody-1)*13 .+ (7:10)] = FiniteDiff.finite_difference_jacobian(d, vector(q3))
        elseif bound_type <: ImpactContact
            J[offr + N½ .+ (1:1), (ibody-1)*13 .+ (1:3)] =  bound.ainv3
            J[offr + N½ .+ (1:1), (ibody-1)*13 .+ (7:10)] = bound.ainv3 * (VLmat(q3) * Lmat(UnitQuaternion(bound.p)) * Tmat() + VRᵀmat(q3) * Rmat(UnitQuaternion(bound.p)))#) * Rmat(quaternion_map(ϕ25, timestep)*timestep/2)*LVᵀmat(q2)
        end
        offr += length(contact)
    end
    return J
end

function full_data_matrix(mechanism::Mechanism{T,Nn,Ne,Nb}; attjac::Bool = true) where {T,Nn,Ne,Nb}
    mechanism = deepcopy(mechanism)
    timestep = mechanism.timestep
    system = mechanism.system
    joints = mechanism.joints
    contacts = mechanism.contacts
    jointids = getfield.(joints, :id)

    resdims = [length(system.vector_entries[i].value) for i=1:Nn]
    joint_dimensions = length.(joints)
    contact_dimensions = length.(contacts)
    nu = control_dimension(mechanism)
    njoints = joint_dimension(mechanism)
    ncontacts = contact_dimension(mechanism)
    nic = attjac ? 12Nb : 13Nb # initial conditions x2, v1, q2, ϕ15

    Fz, Fu = dynamics_jacobian(mechanism, jointids)
    data = get_data(mechanism)
    solution = get_solution(mechanism)
    G = attitude_jacobian(data, Nb)[1:13Nb,1:12Nb]
    H = integrator_jacobian(data, solution, timestep, Nb, njoints, attjac = attjac)[1:13Nb,1:nic]

    B = joint_constraint_jacobian(mechanism) * H
    D = contact_dynamics_jacobian(mechanism) * H
    E = contact_constraint_jacobian(mechanism) * H

    C = Fz + joint_dynamics_jacobian(mechanism) + springapply_damperjacobian(mechanism)
    attjac && (C = C * G)

    A = zeros(sum(resdims), data_dimension(mechanism, attjac = attjac))
    A[1:njoints, 1:nic] += B
    A[njoints .+ (1:6Nb), 1:nic] += C
    A[njoints .+ (1:6Nb), 1:nic] += D
    A[njoints .+ (1:6Nb), nic .+ (1:nu)] += Fu
    A[njoints + 6Nb .+ (1:ncontacts), 1:nic] += E
    return A
end

function attitude_jacobian(data::AbstractVector, Nb::Int)
    G = zeros(0,0)
    for i = 1:Nb
        x2, v15, q2, ϕ15 = unpack_data(data[13*(i-1) .+ (1:13)])
        q2 = UnitQuaternion(q2..., false)
        G = cat(G, I(6), LVᵀmat(q2), I(3), dims = (1,2))
    end
    ndata = length(data)
    nu = ndata - size(G)[1]
    G = cat(G, I(nu), dims = (1,2))
    return G
end

function integrator_jacobian(data::AbstractVector, sol::AbstractVector, timestep, Nb::Int, njoints::Int; attjac::Bool = true)
    H = zeros(0,0)
    for i = 1:Nb
        x2, v15, q2, ϕ15 = unpack_data(data[13*(i-1) .+ (1:13)])
        ϕ25 = sol[njoints + 6*(i-1) + 3 .+ SVector{3,Int}(1:3)]
        q2 = UnitQuaternion(q2..., false)
        H = cat(H, I(6), ∂integrator∂q(q2, ϕ25, timestep, attjac = attjac), I(3), dims = (1,2))
    end
    ndata = length(data)
    nu = ndata - size(H)[1]
    H = cat(H, I(nu), dims = (1,2))
    return H
end

function getλJoint(joint::JointConstraint{T,N,Nc}, i::Int) where {T,N,Nc}
    n1 = 1
    for j = 1:i-1
        n1 += length(joint.constraints[j])
    end
    n2 = n1 - 1 + length(joint.constraints[i])

    λi = SVector{n2-n1+1,T}(joint.λsol[2][n1:n2])
    return λi
end
