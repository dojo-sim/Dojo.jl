function joint_constraint_jacobian(mechanism::Mechanism{T,Nn,Ne,Nb}) where {T,Nn,Ne,Nb}
    Δt = mechanism.Δt
    eqcs = mechanism.eqconstraints
    nc = sum(length.(eqcs))
    Gl = zeros(T,nc,13Nb)

    oneindc = 0
    for eqc in eqcs
        ind1 = 1
        ind2 = 0

        parentid = eqc.parentid
        if parentid !== nothing
            parentind = parentid - Ne
            pbody = getbody(mechanism,parentid)
            pstate = pbody.state

            for (i,childid) in enumerate(eqc.childids)
                childind = childid - Ne
                cbody = getbody(mechanism,childid)
                cstate = cbody.state
                joint = eqc.constraints[i]

                ind2 += length(joint)
                range = oneindc+ind1:oneindc+ind2

                pcol13 = offsetrange(parentind,13)
                ccol13 = offsetrange(childind,13)

                pXl, pQl = ∂g∂posa(joint, pbody, cbody, eqc.λsol[2], Δt) # x3
                cXl, cQl = ∂g∂posb(joint, pbody, cbody, eqc.λsol[2], Δt) # x3

                pGlx = pXl
                pGlq = pQl
                cGlx = cXl
                cGlq = cQl
                Gl[range, pcol13[1:3]] = pGlx
                Gl[range, pcol13[7:10]] = pGlq
                Gl[range,ccol13[1:3]] = cGlx
                Gl[range,ccol13[7:10]] = cGlq
                ind1 = ind2+1
            end
        else
            for (i,childid) in enumerate(eqc.childids)
                childind = childid - Ne
                cbody = getbody(mechanism,childid)
                cstate = cbody.state
                joint = eqc.constraints[i]
                ind2 += length(joint)
                range = oneindc+ind1:oneindc+ind2

                ccol13 = offsetrange(childind,13)

                cXl, cQl = ∂g∂posb(joint, mechanism.origin, cbody, eqc.λsol[2], Δt) # x3

                cGlx = cXl
                cGlq = cQl
                # if typeof(joint) <: Translational
                #     @show A
                #     # @show cGlx
                #     @show cGlq
                # end
                Gl[range,ccol13[1:3]] = cGlx
                Gl[range,ccol13[7:10]] = cGlq
                ind1 = ind2+1
            end
        end
        oneindc += length(eqc)
    end
    return Gl
end

function joint_dynamics_jacobian(mechanism::Mechanism{T,Nn,Ne,Nb}) where {T,Nn,Ne,Nb}
    Δt = mechanism.Δt
    J = zeros(T,6Nb,13Nb)

    for eqc in mechanism.eqconstraints
        parentid = eqc.parentid

        if parentid !== nothing
            parentind = parentid - Ne
            pbody = getbody(mechanism, parentid)
            pstate = pbody.state
            prow6 = offsetrange(parentind,6)
            pcol13 = offsetrange(parentind,13)

            for (i, childid) in enumerate(eqc.childids)
                childind = childid - Ne
                cbody = getbody(mechanism, childid)
                cstate = cbody.state
                joint = eqc.constraints[i]
                crow6 = offsetrange(childind,6)
                ccol13 = offsetrange(childind,13)
                λ = getλJoint(eqc, i)

                Aaa = zeros(T,6,13)
                Aab = zeros(T,6,13)
                Aba = zeros(T,6,13)
                Abb = zeros(T,6,13)

                xa, qa = posargs2(pstate)
                xb, qb = posargs2(cstate)

                # XX, XQ, QX, QQ = ∂2g∂posaa(joint, posargs2(pstate)..., posargs2(cstate)...)
                if typeof(∂g∂ʳposa(joint, xa, qa, xb, qb, eqc.λsol[2])) <: AbstractArray && length(joint) > 0
                    XX = FiniteDiff.finite_difference_jacobian(w -> -transpose(∂g∂ʳposa(joint, w, qa, xb, qb, eqc.λsol[2])[:, 1:3]) * λ, xa)
                    XQ = FiniteDiff.finite_difference_jacobian(w -> -transpose(∂g∂ʳposa(joint, xa, UnitQuaternion(w..., false), xb, qb, eqc.λsol[2])[:, 1:3]) * λ, [qa.w; qa.x; qa.y; qa.z])
                    QX = FiniteDiff.finite_difference_jacobian(w -> -transpose(∂g∂ʳposa(joint, w, qa, xb, qb, eqc.λsol[2])[:, 4:6]) * λ, xa)
                    QQ = FiniteDiff.finite_difference_jacobian(w -> -transpose(∂g∂ʳposa(joint, xa, UnitQuaternion(w..., false), xb, qb, eqc.λsol[2])[:, 4:6]) * λ, [qa.w; qa.x; qa.y; qa.z])

                    Aaa[1:3,1:3] = XX
                    Aaa[1:3,7:10] = XQ
                    Aaa[4:6,1:3] = QX
                    Aaa[4:6,7:10] = QQ
                end

                # XX, XQ, QX, QQ = ∂2g∂posab(joint, posargs2(pstate)..., posargs2(cstate)...)

                if typeof(∂g∂ʳposa(joint, xa, qa, xb, qb, eqc.λsol[2])) <: AbstractArray && length(joint) > 0
                    XX = FiniteDiff.finite_difference_jacobian(w -> -transpose(∂g∂ʳposa(joint, xa, qa, w, qb, eqc.λsol[2])[:, 1:3]) * λ, xb)
                    XQ = FiniteDiff.finite_difference_jacobian(w -> -transpose(∂g∂ʳposa(joint, xa, qa, xb, UnitQuaternion(w..., false), eqc.λsol[2])[:, 1:3]) * λ, [qb.w; qb.x; qb.y; qb.z])
                    QX = FiniteDiff.finite_difference_jacobian(w -> -transpose(∂g∂ʳposa(joint, xa, qa, w, qb, eqc.λsol[2])[:, 4:6]) * λ, xb)
                    QQ = FiniteDiff.finite_difference_jacobian(w -> -transpose(∂g∂ʳposa(joint, xa, qa, xb, UnitQuaternion(w..., false), eqc.λsol[2])[:, 4:6]) * λ, [qb.w; qb.x; qb.y; qb.z])

                    Aab[1:3,1:3] = XX
                    Aab[1:3,7:10] = XQ
                    Aab[4:6,1:3] = QX
                    Aab[4:6,7:10] = QQ
                end

                # XX, XQ, QX, QQ = ∂2g∂posba(joint, posargs2(pstate)..., posargs2(cstate)...)

                if typeof(∂g∂ʳposb(joint, xa, qa, xb, qb, eqc.λsol[2])) <: AbstractArray && length(joint) > 0
                    XX = FiniteDiff.finite_difference_jacobian(w -> -transpose(∂g∂ʳposb(joint, w, qa, xb, qb, eqc.λsol[2])[:, 1:3]) * λ, xa)
                    XQ = FiniteDiff.finite_difference_jacobian(w -> -transpose(∂g∂ʳposb(joint, xa, UnitQuaternion(w..., false), xb, qb, eqc.λsol[2])[:, 1:3]) * λ, [qa.w; qa.x; qa.y; qa.z])
                    QX = FiniteDiff.finite_difference_jacobian(w -> -transpose(∂g∂ʳposb(joint, w, qa, xb, qb, eqc.λsol[2])[:, 4:6]) * λ, xa)
                    QQ = FiniteDiff.finite_difference_jacobian(w -> -transpose(∂g∂ʳposb(joint, xa, UnitQuaternion(w..., false), xb, qb, eqc.λsol[2])[:, 4:6]) * λ, [qa.w; qa.x; qa.y; qa.z])

                    Aba[1:3,1:3] = XX
                    Aba[1:3,7:10] = XQ
                    Aba[4:6,1:3] = QX
                    Aba[4:6,7:10] = QQ
                end

                # XX, XQ, QX, QQ = ∂2g∂posbb(joint, posargs2(pstate)..., posargs2(cstate)...)

                if typeof(∂g∂ʳposb(joint, xa, qa, xb, qb, eqc.λsol[2])) <: AbstractArray && length(joint) > 0
                    XX = FiniteDiff.finite_difference_jacobian(w -> -transpose(∂g∂ʳposb(joint, xa, qa, w, qb, eqc.λsol[2])[:, 1:3]) * λ, xb)
                    XQ = FiniteDiff.finite_difference_jacobian(w -> -transpose(∂g∂ʳposb(joint, xa, qa, xb, UnitQuaternion(w..., false), eqc.λsol[2])[:, 1:3]) * λ, [qb.w; qb.x; qb.y; qb.z])
                    QX = FiniteDiff.finite_difference_jacobian(w -> -transpose(∂g∂ʳposb(joint, xa, qa, w, qb, eqc.λsol[2])[:, 4:6]) * λ, xb)
                    QQ = FiniteDiff.finite_difference_jacobian(w -> -transpose(∂g∂ʳposb(joint, xa, qa, xb, UnitQuaternion(w..., false), eqc.λsol[2])[:, 4:6]) * λ, [qb.w; qb.x; qb.y; qb.z])

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
            for (i, childid) in enumerate(eqc.childids)
                childind = childid - Ne
                cbody = getbody(mechanism, childid)
                cstate = cbody.state
                joint = eqc.constraints[i]
                crow6 = offsetrange(childind,6)
                ccol13 = offsetrange(childind,13)
                λ = getλJoint(eqc, i)

                Abb = zeros(T,6,13)

                xb, qb = posargs2(cstate)

                if typeof(∂g∂ʳposb(joint, xb, qb, eqc.λsol[2])) <: AbstractArray && length(joint) > 0

                    XX = FiniteDiff.finite_difference_jacobian(w -> -transpose(∂g∂ʳposb(joint, w, qb, eqc.λsol[2])[:, 1:3]) * λ, xb)
                    XQ = FiniteDiff.finite_difference_jacobian(w -> -transpose(∂g∂ʳposb(joint, xb, UnitQuaternion(w..., false), eqc.λsol[2])[:, 1:3]) * λ, [qb.w; qb.x; qb.y; qb.z])
                    QX = FiniteDiff.finite_difference_jacobian(w -> -transpose(∂g∂ʳposb(joint, w, qb, eqc.λsol[2])[:, 4:6]) * λ, xb)
                    QQ = FiniteDiff.finite_difference_jacobian(w -> -transpose(∂g∂ʳposb(joint, xb, UnitQuaternion(w..., false), eqc.λsol[2])[:, 4:6]) * λ, [qb.w; qb.x; qb.y; qb.z])

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

function spring_damper_jacobian(mechanism::Mechanism{T,Nn,Ne,Nb}) where {T,Nn,Ne,Nb}
    Δt = mechanism.Δt
    J = zeros(T,6Nb,13Nb)

    for eqc in mechanism.eqconstraints
        parentid = eqc.parentid

        if parentid !== nothing
            parentind = parentid - Ne
            pbody = getbody(mechanism, parentid)
            pstate = pbody.state
            prow6 = offsetrange(parentind,6)
            pcol13 = offsetrange(parentind,13)

            for (i, childid) in enumerate(eqc.childids)
                childind = childid - Ne
                cbody = getbody(mechanism, childid)
                cstate = cbody.state
                joint = eqc.constraints[i]
                crow6 = offsetrange(childind,6)
                ccol13 = offsetrange(childind,13)
                λ = getλJoint(eqc, i)

                Aaa = zeros(T,6,13)
                Aab = zeros(T,6,13)
                Aba = zeros(T,6,13)
                Abb = zeros(T,6,13)

                Aaa[:, [1:3; 7:10]] -= ∂springforcea∂posa(joint, pbody, cbody, Δt, attjac = false)
                Aaa[:, [1:3; 7:10]] -= ∂damperforcea∂posa(joint, pbody, cbody, Δt, attjac = false)

                Aab[:, [1:3; 7:10]] -= ∂springforcea∂posb(joint, pbody, cbody, Δt, attjac = false)
                Aab[:, [1:3; 7:10]] -= ∂damperforcea∂posb(joint, pbody, cbody, Δt, attjac = false)

                Aba[:, [1:3; 7:10]] -= ∂springforceb∂posa(joint, pbody, cbody, Δt, attjac = false)
                Aba[:, [1:3; 7:10]] -= ∂damperforceb∂posa(joint, pbody, cbody, Δt, attjac = false)

                Abb[:, [1:3; 7:10]] -= ∂springforceb∂posb(joint, pbody, cbody, Δt, attjac = false)
                Abb[:, [1:3; 7:10]] -= ∂damperforceb∂posb(joint, pbody, cbody, Δt, attjac = false)

                J[prow6,pcol13] += Aaa
                J[prow6,ccol13] += Aab
                J[crow6,pcol13] += Aba
                J[crow6,ccol13] += Abb
            end
        else
            for (i, childid) in enumerate(eqc.childids)
                pbody = mechanism.origin
                childind = childid - Ne
                cbody = getbody(mechanism, childid)
                cstate = cbody.state
                joint = eqc.constraints[i]
                crow6 = offsetrange(childind,6)
                ccol13 = offsetrange(childind,13)
                λ = getλJoint(eqc, i)

                Abb = zeros(T,6,13)

                Abb[:, [1:3; 7:10]] -= ∂springforceb∂posb(joint, pbody, cbody, Δt, attjac = false)
                Abb[:, [1:3; 7:10]] -= ∂damperforceb∂posb(joint, pbody, cbody, Δt, attjac = false)

                J[crow6,ccol13] += Abb
            end
        end
    end
    return J
end

function control_datamat(mechanism::Mechanism{T,Nn,Ne,Nb}) where {T,Nn,Ne,Nb}
    Δt = mechanism.Δt
    Fzu = zeros(T,13Nb,12Nb)
    # Fzu = zeros(T,13Nb,13Nb)

    for eqc in mechanism.eqconstraints
        parentid = eqc.parentid
        for (i,childid) in enumerate(eqc.childids)
            childind = childid - Ne
            joint = eqc.constraints[i]
            if parentid != nothing
                parentind = parentid - Ne
                pbody = getbody(mechanism, parentid)
                cbody = getbody(mechanism, childid)

                FaXa, FaQa, τaXa, τaQa, FbXa, FbQa, τbXa, τbQa = ∂Fτ∂posa(joint, pbody.state, cbody.state, Δt)
                FaXb, FaQb, τaXb, τaQb, FbXb, FbQb, τbXb, τbQb = ∂Fτ∂posb(joint, pbody.state, cbody.state, Δt)

                xa, qa = posargs2(pbody.state)
                Ma = [I zeros(3,3); zeros(4,3) LVᵀmat(qa)]
                # @show size(Ma)
                xb, qb = posargs2(cbody.state)
                Mb = [I zeros(3,3); zeros(4,3) LVᵀmat(qb)]

                cola6 = offsetrange(parentind,6)
                colb6 = offsetrange(childind,6)
                # @show cola6
                # @show colb6
                rowav = offsetrange(parentind,3,13,2)
                rowaω = offsetrange(parentind,3,13,4).+1
                rowbv = offsetrange(childind,3,13,2)
                rowbω = offsetrange(childind,3,13,4).+1

                Fzu[[rowav; rowaω],cola6] = [FaXa FaQa; τaXa τaQa] * Ma
                # @show [FaXa FaQa; τaXa τaQa]
                Fzu[[rowbv; rowbω],cola6] = [FbXa FbQa; τbXa τbQa] * Ma
                # @show [FbXa FbQa; τbXa τbQa]
                Fzu[[rowav; rowaω],colb6] = [FaXb FaQb; τaXb τaQb] * Mb
                # @show [FaXb FaQb; τaXb τaQb]
                Fzu[[rowbv; rowbω],colb6] = [FbXb FbQb; τbXb τbQb] * Mb
                # @show [FbXb FbQb; τbXb τbQb]
            else
                pbody = mechanism.origin
                cbody = getbody(mechanism, childid)
                FaXb, FaQb, τaXb, τaQb, FbXb, FbQb, τbXb, τbQb = ∂Fτ∂posb(joint, cbody.state, Δt)

                colb6 = offsetrange(childind,6)
                # @show colb6
                rowbv = offsetrange(childind,3,13,2)
                rowbω = offsetrange(childind,3,13,4).+1

                xb, qb = posargs2(cbody.state)
                Mb = [I zeros(3,3); zeros(4,3) LVᵀmat(qb)]

                Fzu[[rowbv; rowbω], colb6] = [FbXb FbQb; τbXb τbQb] * Mb
                # @show [FbXb FbQb; τbXb τbQb]

            end
        end
    end
    return Fzu
end

function ∂body∂z(body::Body{T}, Δt::T; attjac::Bool = true) where T
    state = body.state
    q2 = state.q2[1]
    ϕ25 = state.ϕsol[2]
    Z3 = szeros(T,3,3)
    Z34 = szeros(T,3,4)
    ZT = attjac ? szeros(T,6,6) : szeros(T,6,7)
    ZR = szeros(T,7,6)

    x1, q1 = posargs1(state)
    x2, q2 = posargs2(state)
    x3, q3 = posargs3(state, Δt)

    AposT = [-I Z3]
    AvelT = [Z3 -I*body.m] # solving for impulses

    AposR = [-∂integrator∂q(q2, ϕ25, Δt, attjac = attjac) szeros(4,3)]
    
    rot_q1(q) = -4 / Δt * LVᵀmat(q2)' * Lmat(UnitQuaternion(q..., false)) * Vᵀmat() * body.J * Vmat() * Lmat(UnitQuaternion(q..., false))' * vector(q2)
    rot_q2(q) = -4 / Δt * LVᵀmat(UnitQuaternion(q..., false))' * Tmat() * Rmat(getq3(UnitQuaternion(q..., false), state.ϕsol[2], Δt))' * Vᵀmat() * body.J * Vmat() * Lmat(UnitQuaternion(q..., false))' * vector(getq3(UnitQuaternion(q..., false), state.ϕsol[2], Δt)) + -4 / Δt * LVᵀmat(UnitQuaternion(q..., false))' * Lmat(getq3(UnitQuaternion(q..., false), -state.ϕ15, Δt)) * Vᵀmat() * body.J * Vmat() * Lmat(getq3(UnitQuaternion(q..., false), -state.ϕ15, Δt))' * q

    dynR_ϕ15 = -1.0 * FiniteDiff.finite_difference_jacobian(rot_q1, vector(q1)) * ∂integrator∂ϕ(q2, -state.ϕ15, Δt)
    dynR_q2 = FiniteDiff.finite_difference_jacobian(rot_q2, vector(q2))
    AvelR = attjac ? [dynR_q2 * LVᵀmat(q2) dynR_ϕ15] : [dynR_q2 dynR_ϕ15]

    # AposR = attjac ? [-Rmat(ωbar(state.ϕ15, Δt)*Δt/2)*LVᵀmat(state.q1) -Lmat(state.q1)*derivωbar(state.ϕ15, Δt)*Δt/2] : [-Rmat(ωbar(state.ϕ15, Δt)*Δt/2) -Lmat(state.q1)*derivωbar(state.ϕ15, Δt)*Δt/2]
    # J = body.J
    # ω1 = state.ϕ15
    # sq1 = sqrt(4 / Δt^2 - ω1' * ω1)
    # ω1func = -skewplusdiag(-ω1, sq1) * J + J * ω1 * (ω1' / sq1) - skew(J * ω1)

    # AvelR = [(attjac ? Z3 : Z34) ω1func*Δt] # solving for impulses

    return [[AposT;AvelT] ZT;
             ZR [AposR;AvelR]]
end

function ∂body∂u(body::Body{T}, Δt) where T
    Z3 = szeros(T,3,3)
    Z43 = szeros(T,4,3)

    BposT = [Z3 Z3]
    BvelT = [-I Z3]
    BposR = [Z43 Z43]
    BvelR = [Z3 -2.0 * I]
    return [BposT;BvelT;BposR;BvelR]
end

function dynamics_jacobian(mechanism::Mechanism{T,Nn,Ne,Nb}, eqcids) where {T,Nn,Ne,Nb}
    Δt = mechanism.Δt
    nu = controldim(mechanism)

    # get state linearization
    Fz = zeros(T,6Nb,13Nb)
    Fu = zeros(T,6Nb,6Nb)

    Bcontrol = zeros(T,6Nb,nu)

    for (i, body) in enumerate(mechanism.bodies)
        col6 = offsetrange(i,6)
        col13 = offsetrange(i,13)

        Fzi = ∂body∂z(body, Δt, attjac = false)[[4:6; 11:13],:]
        Fui = ∂body∂u(body, Δt)[[4:6; 11:13],:]

        Fz[col6,col13] = Fzi
        Fu[col6,col6] = Fui
    end

    ############################################################################
    # @warn "control datamat seems wrong"
    # control datamat is always zero and seems to be filled only for the first half
    # Fz = Fz * attitudejacobian(getdata(mechanism), Nb)[1:13Nb,1:12Nb]
    # Fz += control_datamat(mechanism)
    # @show norm(control_datamat(mechanism))
    ############################################################################

    n1 = 1
    n2 = 0

    for id in eqcids
        eqc = geteqconstraint(mechanism, id)
        n2 += controldim(eqc)

        parentid = eqc.parentid
        if parentid != nothing
            parentind = parentid - Ne
            col6 = offsetrange(parentind,6)
            Bcontrol[col6,n1:n2] = ∂Fτ∂ua(mechanism, eqc, getbody(mechanism, parentid))
        end
        for childid in eqc.childids
            childind = childid - Ne
            col6 = offsetrange(childind,6)
            Bcontrol[col6,n1:n2] = ∂Fτ∂ub(mechanism, eqc, getbody(mechanism, childid))
        end

        n1 = n2+1
    end
    return Fz, Fu * Bcontrol
end

function contact_dynamics_jacobian(mechanism::Mechanism{T,Nn,Ne,Nb}) where {T,Nn,Ne,Nb}
    Δt = mechanism.Δt
    J = zeros(6Nb, 13Nb)

    offr = 0
    offc = 0
    for body in mechanism.bodies
        for ineqc in mechanism.ineqconstraints
            if ineqc.parentid == body.id
                bound = ineqc.constraints[1]
                bound_type = typeof(bound)
                x3, q3 = posargs3(body.state, Δt)

                function d(vars)
                    x = vars[1:3]
                    q = UnitQuaternion(vars[4:7]..., false)
                    return ∂g∂ʳpos(bound, x, q, nothing)' * ineqc.γsol[2]
                end

                if bound_type <: ContactBound
                    # J[offr .+ (1:6), offc .+ [1:3; 7:10]] -= _dN(x3, vector(q3), ineqc.γsol[2][1:1], bound.p)
                    # J[offr .+ (1:6), offc .+ [1:3; 7:10]] -= _dB(x3, vector(q3), ineqc.γsol[2][2:4], bound.p)
                    J[offr .+ (1:6), offc .+ [1:3; 7:10]] -= FiniteDiff.finite_difference_jacobian(d, [x3; vector(q3)])
                elseif bound_type <: LinearContactBound
                    J[offr .+ (1:6), offc .+ [1:3; 7:10]] -= FiniteDiff.finite_difference_jacobian(d, [x3; vector(q3)])
                elseif bound_type <: ImpactBound
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
    Δt = mechanism.Δt
    ineqcs = mechanism.ineqconstraints
    nineqcs = ineqcdim(mechanism)
    J = zeros(nineqcs, 13Nb)

    offr = 0
    for ineqc in ineqcs
        bound = ineqc.constraints[1]
        body = getbody(mechanism, ineqc.parentid)
        N½ = Int(length(ineqc)/2)
        x2, v25, q2, ϕ25 = fullargssol(body.state)
        x3, q3 = posargs3(body.state, Δt)
        ibody = findfirst(x -> x == body.id, mechanism.bodies.keys)
        bound_type = typeof(bound)

        function d(vars)
            q3 = UnitQuaternion(vars..., false)
            # Bxmat = bound.Bx
            # Bqmat = Bxmat * ∂vrotate∂q(bound.p, q) * LVᵀmat(q)
            # return Bqmat * ϕ25
            vp = v25 + skew(vrotate(ϕ25, q3)) * (vrotate(bound.p, q3) - bound.offset)
            return bound.Bx * vp
        end
        if bound_type <: ContactBound
            J[offr + N½ .+ (1:1), (ibody-1)*13 .+ (1:3)] =  bound.ainv3
            J[offr + N½ .+ (1:1), (ibody-1)*13 .+ (7:10)] = bound.ainv3 * (VLmat(q3) * Lmat(UnitQuaternion(bound.p)) * Tmat() + VRᵀmat(q3) * Rmat(UnitQuaternion(bound.p)))#) * Rmat(ωbar(ϕ25, Δt)*Δt/2)*LVᵀmat(q2)
            J[offr + N½ .+ (3:4), (ibody-1)*13 .+ (7:10)] = FiniteDiff.finite_difference_jacobian(d, vector(q3))#dBω(vector(q3), ϕ25, bound.p)
        elseif bound_type <: LinearContactBound
            J[offr + N½ .+ (1:1), (ibody-1)*13 .+ (1:3)] =  bound.ainv3
            J[offr + N½ .+ (1:1), (ibody-1)*13 .+ (7:10)] = bound.ainv3 * ∂vrotate∂q(bound.p, q3)
            J[offr + N½ .+ (3:6), (ibody-1)*13 .+ (7:10)] = FiniteDiff.finite_difference_jacobian(d, vector(q3))
        elseif bound_type <: ImpactBound
            J[offr + N½ .+ (1:1), (ibody-1)*13 .+ (1:3)] =  bound.ainv3
            J[offr + N½ .+ (1:1), (ibody-1)*13 .+ (7:10)] = bound.ainv3 * (VLmat(q3) * Lmat(UnitQuaternion(bound.p)) * Tmat() + VRᵀmat(q3) * Rmat(UnitQuaternion(bound.p)))#) * Rmat(ωbar(ϕ25, Δt)*Δt/2)*LVᵀmat(q2)
        end
        offr += length(ineqc)
    end
    return J
end

function dBω(q, ω, p)
    q₁, q₂, q₃, q₄ = q
    ω₁, ω₂, ω₃ = ω
    p₁, p₂, p₃ = p
    dB = zeros(2, 4)

    dB[1, 1] = ω₁*(4.0p₂*q₃ + 4.0p₃*q₄) + ω₂*(4.0p₃*q₁ - (4.0p₁*q₃)) + ω₃*(-4.0p₁*q₄ - (4.0p₂*q₁))
    dB[2, 1] = ω₁*(-4.0p₂*q₂ - (4.0p₃*q₁)) + ω₂*(4.0p₁*q₂ + 4.0p₃*q₄) + ω₃*(4.0p₁*q₁ - (4.0p₂*q₄))

    dB[1, 2] = ω₁*(4.0p₂*q₄ - (4.0p₃*q₃)) + ω₂*(4.0p₃*q₂ - (4.0p₁*q₄)) + ω₃*(4.0p₁*q₃ - (4.0p₂*q₂))
    dB[2, 2] = ω₁*(4.0p₃*q₂ - (4.0p₂*q₁)) + ω₂*(4.0p₁*q₁ + 4.0p₃*q₃) + ω₃*(-4.0p₁*q₂ - (4.0p₂*q₃))

    dB[1, 3] = ω₁*(4.0p₂*q₁ - (4.0p₃*q₂)) + ω₂*(-4.0p₁*q₁ - (4.0p₃*q₃)) + ω₃*(4.0p₁*q₂ + 4.0p₂*q₃)
    dB[2, 3] = ω₁*(4.0p₂*q₄ - (4.0p₃*q₃)) + ω₂*(4.0p₃*q₂ - (4.0p₁*q₄)) + ω₃*(4.0p₁*q₃ - (4.0p₂*q₂))

    dB[1, 4] = ω₁*(4.0p₂*q₂ + 4.0p₃*q₁) + ω₂*(-4.0p₁*q₂ - (4.0p₃*q₄)) + ω₃*(4.0p₂*q₄ - (4.0p₁*q₁))
    dB[2, 4] = ω₁*(4.0p₂*q₃ + 4.0p₃*q₄) + ω₂*(4.0p₃*q₁ - (4.0p₁*q₃)) + ω₃*(-4.0p₁*q₄ - (4.0p₂*q₁))

    return dB
end

function _dB(x, q, b, p)
    p₁, p₂, p₃ = p
    z₁, z₂, z₃ = x
    z₄, z₅, z₆, z₇ = q
    b₁, b₂, b₃ = b
    dBb = zeros(6, 7)

    dBb[4,4] = b₂*(4.0p₂*z₆ + 4.0p₃*z₇) + b₃*(-4.0p₂*z₅ - (4.0p₃*z₄))
    dBb[4,5] = b₂*(4.0p₂*z₇ - (4.0p₃*z₆)) + b₃*(4.0p₃*z₅ - (4.0p₂*z₄))
    dBb[4,6] = b₂*(4.0p₂*z₄ - (4.0p₃*z₅)) + b₃*(4.0p₂*z₇ - (4.0p₃*z₆))
    dBb[4,7] = b₂*(4.0p₂*z₅ + 4.0p₃*z₄) + b₃*(4.0p₂*z₆ + 4.0p₃*z₇)

    dBb[5,4] = b₂*(4.0p₃*z₄ - (4.0p₁*z₆)) + b₃*(4.0p₁*z₅ + 4.0p₃*z₇)
    dBb[5,5] = b₂*(4.0p₃*z₅ - (4.0p₁*z₇)) + b₃*(4.0p₁*z₄ + 4.0p₃*z₆)
    dBb[5,6] = b₂*(-4.0p₁*z₄ - (4.0p₃*z₆)) + b₃*(4.0p₃*z₅ - (4.0p₁*z₇))
    dBb[5,7] = b₂*(-4.0p₁*z₅ - (4.0p₃*z₇)) + b₃*(4.0p₃*z₄ - (4.0p₁*z₆))

    dBb[6,4] = b₂*(-4.0p₁*z₇ - (4.0p₂*z₄)) + b₃*(4.0p₁*z₄ - (4.0p₂*z₇))
    dBb[6,5] = b₂*(4.0p₁*z₆ - (4.0p₂*z₅)) + b₃*(-4.0p₁*z₅ - (4.0p₂*z₆))
    dBb[6,6] = b₂*(4.0p₁*z₅ + 4.0p₂*z₆) + b₃*(4.0p₁*z₆ - (4.0p₂*z₅))
    dBb[6,7] = b₂*(4.0p₂*z₇ - (4.0p₁*z₄)) + b₃*(-4.0p₁*z₇ - (4.0p₂*z₄))

    return dBb
end

dB(x, q, b, p) = _dB(x, q, b, p) * [I zeros(3,3); zeros(4,3) G(q)]

function _dN(x, q, γ, p)
    p₁, p₂, p₃ = p
    z₁, z₂, z₃ = x
    γ₁ = γ[1]
    z₄, z₅, z₆, z₇ = q
    dNγ = zeros(6, 7)

    dNγ[4,4] = γ₁*(4.0p₂*z₄ - (4.0p₃*z₅))
    dNγ[4,5] = γ₁*(-4.0p₂*z₅ - (4.0p₃*z₄))
    dNγ[4,6] = γ₁*(-4.0p₂*z₆ - (4.0p₃*z₇))
    dNγ[4,7] = γ₁*(4.0p₂*z₇ - (4.0p₃*z₆))

    dNγ[5,4] = γ₁*(-4.0p₁*z₄ - (4.0p₃*z₆))
    dNγ[5,5] = γ₁*(4.0p₁*z₅ + 4.0p₃*z₇)
    dNγ[5,6] = γ₁*(4.0p₁*z₆ - (4.0p₃*z₄))
    dNγ[5,7] = γ₁*(4.0p₃*z₅ - (4.0p₁*z₇))

    dNγ[6,4] = γ₁*(4.0p₁*z₅ + 4.0p₂*z₆)
    dNγ[6,5] = γ₁*(4.0p₁*z₄ - (4.0p₂*z₇))
    dNγ[6,6] = γ₁*(4.0p₁*z₇ + 4.0p₂*z₄)
    dNγ[6,7] = γ₁*(4.0p₁*z₆ - (4.0p₂*z₅))

    return dNγ
end

dN(x, q, γ, p) = _dN(x, q, γ, p) * [I zeros(3,3); zeros(4,3) G(q)]

function _dG(joint::Rotational{T,N}, x, q, γ, p) where {T,N}
    p₁, p₂, p₃ = p
    z₁, z₂, z₃ = x
    γ₁ = γ[1]
    z₄, z₅, z₆, z₇ = q
    dGγ = zeros(N, 7)

    dGγ[4,4] = γ₁*(4.0p₂*z₄ - (4.0p₃*z₅))
    dGγ[4,5] = γ₁*(-4.0p₂*z₅ - (4.0p₃*z₄))
    dGγ[4,6] = γ₁*(-4.0p₂*z₆ - (4.0p₃*z₇))
    dGγ[4,7] = γ₁*(4.0p₂*z₇ - (4.0p₃*z₆))

    dGγ[5,4] = γ₁*(-4.0p₁*z₄ - (4.0p₃*z₆))
    dGγ[5,5] = γ₁*(4.0p₁*z₅ + 4.0p₃*z₇)
    dGγ[5,6] = γ₁*(4.0p₁*z₆ - (4.0p₃*z₄))
    dGγ[5,7] = γ₁*(4.0p₃*z₅ - (4.0p₁*z₇))

    dGγ[6,4] = γ₁*(4.0p₁*z₅ + 4.0p₂*z₆)
    dGγ[6,5] = γ₁*(4.0p₁*z₄ - (4.0p₂*z₇))
    dGγ[6,6] = γ₁*(4.0p₁*z₇ + 4.0p₂*z₄)
    dGγ[6,7] = γ₁*(4.0p₁*z₆ - (4.0p₂*z₅))

    return dGγ
end

function _dG(joint::Translational{T,N}, x, q, γ, p) where {T,N}
    p₁, p₂, p₃ = p
    z₁, z₂, z₃ = x
    γ₁ = γ[1]
    z₄, z₅, z₆, z₇ = q
    dGγ = zeros(N, 7)

    dGγ[4,4] = γ₁*(4.0p₂*z₄ - (4.0p₃*z₅))
    dGγ[4,5] = γ₁*(-4.0p₂*z₅ - (4.0p₃*z₄))
    dGγ[4,6] = γ₁*(-4.0p₂*z₆ - (4.0p₃*z₇))
    dGγ[4,7] = γ₁*(4.0p₂*z₇ - (4.0p₃*z₆))

    dGγ[5,4] = γ₁*(-4.0p₁*z₄ - (4.0p₃*z₆))
    dGγ[5,5] = γ₁*(4.0p₁*z₅ + 4.0p₃*z₇)
    dGγ[5,6] = γ₁*(4.0p₁*z₆ - (4.0p₃*z₄))
    dGγ[5,7] = γ₁*(4.0p₃*z₅ - (4.0p₁*z₇))

    dGγ[6,4] = γ₁*(4.0p₁*z₅ + 4.0p₂*z₆)
    dGγ[6,5] = γ₁*(4.0p₁*z₄ - (4.0p₂*z₇))
    dGγ[6,6] = γ₁*(4.0p₁*z₇ + 4.0p₂*z₄)
    dGγ[6,7] = γ₁*(4.0p₁*z₆ - (4.0p₂*z₅))

    return dGγ
end

dG(joint::Joint, x, q, γ, p) = _dG(joint::Joint, x, q, γ, p) * [I zeros(3,3); zeros(4,3) G(q)]

function full_data_matrix(mechanism::Mechanism{T,Nn,Ne,Nb}; attjac::Bool = true) where {T,Nn,Ne,Nb}
    mechanism = deepcopy(mechanism)
    Δt = mechanism.Δt
    system = mechanism.system
    eqcs = mechanism.eqconstraints
    ineqcs = mechanism.ineqconstraints
    eqcids = getfield.(eqcs, :id)

    resdims = [length(system.vector_entries[i].value) for i=1:Nn]
    eqcdims = length.(eqcs)
    ineqcdims = length.(ineqcs)
    nu = controldim(mechanism)
    neqcs = eqcdim(mechanism)
    nineqcs = ineqcdim(mechanism)
    nic = attjac ? 12Nb : 13Nb # initial conditions x2, v1, q2, ϕ15

    Fz, Fu = dynamics_jacobian(mechanism, eqcids)
    data = getdata(mechanism)
    solution = getsolution(mechanism)
    G = attitudejacobian(data, Nb)[1:13Nb,1:12Nb]
    H = integratorjacobian(data, solution, Δt, Nb, neqcs, attjac = attjac)[1:13Nb,1:nic]

    B = joint_constraint_jacobian(mechanism) * H
    D = contact_dynamics_jacobian(mechanism) * H
    E = contact_constraint_jacobian(mechanism) * H

    C = Fz + joint_dynamics_jacobian(mechanism) + spring_damper_jacobian(mechanism)
    attjac && (C = C * G)

    A = zeros(sum(resdims), datadim(mechanism, attjac = attjac))
    A[1:neqcs, 1:nic] += B
    A[neqcs .+ (1:6Nb), 1:nic] += C
    A[neqcs .+ (1:6Nb), 1:nic] += D
    A[neqcs .+ (1:6Nb), nic .+ (1:nu)] += Fu
    A[neqcs + 6Nb .+ (1:nineqcs), 1:nic] += E
    return A
end

function G(q::AbstractVector)
    qs = q[1]
    qv = q[2:4]
    m = zeros(eltype(q), 4, 3)
    m[1,:] += -qv
    m[2:end,:] += qs * I + skew(qv)
    return m
end

function attitudejacobian(data::AbstractVector, Nb::Int)
    G = zeros(0,0)
    for i = 1:Nb
        x2, v15, q2, ϕ15 = unpackdata(data[13*(i-1) .+ (1:13)])
        q2 = UnitQuaternion(q2..., false)
        G = cat(G, I(6), LVᵀmat(q2), I(3), dims = (1,2))
    end
    ndata = length(data)
    nu = ndata - size(G)[1]
    G = cat(G, I(nu), dims = (1,2))
    return G
end

function integratorjacobian(data::AbstractVector, sol::AbstractVector, Δt, Nb::Int, neqcs::Int; attjac::Bool = true)
    H = zeros(0,0)
    for i = 1:Nb
        x2, v15, q2, ϕ15 = unpackdata(data[13*(i-1) .+ (1:13)])
        ϕ25 = sol[neqcs + 6*(i-1) + 3 .+ SVector{3,Int}(1:3)]
        q2 = UnitQuaternion(q2..., false)
        H = cat(H, I(6), ∂integrator∂q(q2, ϕ25, Δt, attjac = attjac), I(3), dims = (1,2))
    end
    ndata = length(data)
    nu = ndata - size(H)[1]
    H = cat(H, I(nu), dims = (1,2))
    return H
end

function getλJoint(eqc::EqualityConstraint{T,N,Nc}, i::Int) where {T,N,Nc}
    n1 = 1
    for j = 1:i-1
        n1 += length(eqc.constraints[j])
    end
    n2 = n1 - 1 + length(eqc.constraints[i])

    λi = SVector{n2-n1+1,T}(eqc.λsol[2][n1:n2])
    return λi
end

function sensitivities(mechanism, sol, data)
    setdata!(mechanism, data)
    setsolution!(mechanism, sol)
    setentries!(mechanism)
    datamat = full_data_matrix(mechanism)
    solmat = full_matrix(mechanism.system)
    sensi = -1.0 * (solmat \ datamat)
    return sensi
end

function jvp(mechanism, sol, data, v)
    sensi = sensitivities(mechanism, sol, data)
    return sensi * v
end

################################################################################
# Index and Dimensions
################################################################################

function Fz_indices(Nb::Int)
    return vcat([13*(i-1) .+ [4:6; 11:13] for i = 1:Nb]...)
end

function datadim(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}; attjac::Bool = true) where {T,Nn,Ne,Nb,Ni}
    d = 0
    d += 12Nb
    !attjac && (d += Nb)
    d += controldim(mechanism)
    return d
end

function soldim(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}) where {T,Nn,Ne,Nb,Ni}
    d = 0
    d += 6Nb
    d += eqcdim(mechanism)
    d += ineqcdim(mechanism)
    return d
end

function controldim(eqc::EqualityConstraint{T,N,Nc,Cs}; ignore_floating_base::Bool = false) where {T,N,Nc,Cs} 
    ignore_floating_base && (N == 0) && return 0

    N̄ = 0
    for (i, joint) in enumerate(eqc.constraints) 
        N̄ += controldim(joint) 
    end

    return N̄
end

function controldim(joint::Joint{T,N}) where {T,N}
    return 3 - N
end

function controldim(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}; ignore_floating_base::Bool = false) where {T,Nn,Ne,Nb,Ni}
    nu = 0
    for eqc in mechanism.eqconstraints
        nu += controldim(eqc, ignore_floating_base = ignore_floating_base)
    end
    return nu
end

function minCoordDim(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}) where {T,Nn,Ne,Nb,Ni}
    nx = 0
    free_rot_base = false # we are going to check if the link attached to the base has free orientation
    nx = 2 * controldim(mechanism, ignore_floating_base = false)
    free_rot_base && (nx += 1)
    return nx
end

maxCoordDim(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}) where {T,Nn,Ne,Nb,Ni} = 13Nb

function ineqcdim(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}) where {T,Nn,Ne,Nb,Ni}
    nineqcs = 0
    for ineqc in mechanism.ineqconstraints
        nineqcs += length(ineqc)
    end
    return nineqcs
end

function eqcdim(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}) where {T,Nn,Ne,Nb,Ni}
    neqcs = 0
    for eqc in mechanism.eqconstraints
        neqcs += length(eqc)
    end
    return neqcs
end
using Pkg
