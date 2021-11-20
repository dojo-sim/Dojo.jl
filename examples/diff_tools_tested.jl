function joint_datamat(mechanism::Mechanism{T,Nn,Ne,Nb}) where {T,Nn,Ne,Nb}
    Δt = mechanism.Δt
    eqcs = mechanism.eqconstraints
    nc = sum(length.(eqcs))
    Gl = zeros(T,nc,12Nb)

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

                pcol3a12 = offsetrange(parentind,3,12,1)
                pcol3c12 = offsetrange(parentind,3,12,3)

                ccol3a12 = offsetrange(childind,3,12,1)
                ccol3c12 = offsetrange(childind,3,12,3)

                pXl, pQl = ∂g∂posa(joint, pbody, cbody, Δt) # x3
                cXl, cQl = ∂g∂posb(joint, pbody, cbody, Δt) # x3

                A = constraintmat(joint)
                pGlx = A * pXl
                pGlq = A * pQl
                cGlx = A * cXl
                cGlq = A * cQl

                Gl[range,pcol3a12] = pGlx
                Gl[range,pcol3c12] = pGlq * ∂integrator∂q(pstate.qsol[2], pstate.ϕsol[2], Δt)

                Gl[range,ccol3a12] = cGlx
                Gl[range,ccol3c12] = cGlq * ∂integrator∂q(cstate.qsol[2], cstate.ϕsol[2], Δt)
                ind1 = ind2+1
            end
        else
            for (i,childid) in enumerate(eqc.childids)
                childind = childid - Ne
                cbody = getbody(mechanism,childid)
                cstate = cbody.state

                ind2 += length(eqc.constraints[i])
                range = oneindc+ind1:oneindc+ind2

                ccol3a12 = offsetrange(childind,3,12,1)
                ccol3b12 = offsetrange(childind,3,12,2)
                ccol3c12 = offsetrange(childind,3,12,3)
                ccol3d12 = offsetrange(childind,3,12,4)

                cXl, cQl =  ∂g∂posb(eqc.constraints[i], mechanism.origin, cbody, Δt) # x3
                mat = constraintmat(eqc.constraints[i])
                cGlx = mat * cXl
                cGlq = mat * cQl

                Gl[range,ccol3a12] = cGlx
                Gl[range,ccol3c12] = cGlq * ∂integrator∂q(cstate.qsol[2], cstate.ϕsol[2], Δt)
                ind1 = ind2+1
            end
        end
        oneindc += length(eqc)
    end
    return Gl
end

function joint_jacobian_datamat(mechanism::Mechanism{T,Nn,Ne,Nb}) where {T,Nn,Ne,Nb}
    Δt = mechanism.Δt
    FfzG = zeros(T,Nb*13,Nb*13)

    for eqc in mechanism.eqconstraints
        parentid = eqc.parentid

        if parentid !== nothing
            parentind = parentid - Ne
            pbody = getbody(mechanism, parentid)
            pstate = pbody.state
            pcol13 = offsetrange(parentind,13)

            for (i, childid) in enumerate(eqc.childids)
                childind = childid - Ne
                cbody = getbody(mechanism, childid)
                cstate = cbody.state
                joint = eqc.constraints[i]
                ccol13 = offsetrange(childind,13)

                n1 = 1
                n2 = 0
                for j=1:i-1
                    n1 += length(eqc.constraints[j])
                    n2 += length(eqc.constraints[j])
                end
                n2 += length(joint)
                λ = eqc.λsol[2][n1:n2]

                Aaa = zeros(T,13,13)
                Aab = zeros(T,13,13)
                Aba = zeros(T,13,13)
                Abb = zeros(T,13,13)

                xa, qa = posargs2(pstate)
                xb, qb = posargs2(cstate)

                # XX, XQ, QX, QQ = ∂2g∂posaa(joint, posargs2(pstate)..., posargs2(cstate)...)
                λ = SVector{length(λ)}(λ)
                if typeof(∂g∂ʳposa(joint, xa, qa, xb, qb)) <: AbstractArray && length(joint) > 0
                    XX = FiniteDiff.finite_difference_jacobian(w -> -transpose(∂g∂ʳposa(joint, w, qa, xb, qb)[1:3, 1:3]) * constraintmat(joint)' * λ, xa)
                    XQ = FiniteDiff.finite_difference_jacobian(w -> -transpose(∂g∂ʳposa(joint, xa, UnitQuaternion(w..., false), xb, qb)[1:3, 1:3]) * constraintmat(joint)' * λ, [qa.w; qa.x; qa.y; qa.z])
                    QX = FiniteDiff.finite_difference_jacobian(w -> -transpose(∂g∂ʳposa(joint, w, qa, xb, qb)[1:3, 4:6]) * constraintmat(joint)' * λ, xa)
                    QQ = FiniteDiff.finite_difference_jacobian(w -> -transpose(∂g∂ʳposa(joint, xa, UnitQuaternion(w..., false), xb, qb)[1:3, 4:6]) * constraintmat(joint)' * λ, [qa.w; qa.x; qa.y; qa.z])

                    Aaa[4:6,1:3] = XX
                    Aaa[4:6,7:10] = XQ
                    Aaa[11:13,1:3] = QX
                    Aaa[11:13,7:10] = QQ
                end

                # XX, XQ, QX, QQ = ∂2g∂posab(joint, posargs2(pstate)..., posargs2(cstate)...)

                if typeof(∂g∂ʳposa(joint, xa, qa, xb, qb)) <: AbstractArray && length(joint) > 0
                    XX = FiniteDiff.finite_difference_jacobian(w -> -transpose(∂g∂ʳposa(joint, xa, qa, w, qb)[1:3, 1:3]) * constraintmat(joint)' * λ, xb)
                    XQ = FiniteDiff.finite_difference_jacobian(w -> -transpose(∂g∂ʳposa(joint, xa, qa, xb, UnitQuaternion(w..., false))[1:3, 1:3]) * constraintmat(joint)' * λ, [qb.w; qb.x; qb.y; qb.z])
                    QX = FiniteDiff.finite_difference_jacobian(w -> -transpose(∂g∂ʳposa(joint, xa, qa, w, qb)[1:3, 4:6]) * constraintmat(joint)' * λ, xb)
                    QQ = FiniteDiff.finite_difference_jacobian(w -> -transpose(∂g∂ʳposa(joint, xa, qa, xb, UnitQuaternion(w..., false))[1:3, 4:6]) * constraintmat(joint)' * λ, [qb.w; qb.x; qb.y; qb.z])

                    Aab[4:6,1:3] = XX
                    Aab[4:6,7:10] = XQ
                    Aab[11:13,1:3] = QX
                    Aab[11:13,7:10] = QQ
                end

                # XX, XQ, QX, QQ = ∂2g∂posba(joint, posargs2(pstate)..., posargs2(cstate)...)

                if typeof(∂g∂ʳposb(joint, xa, qa, xb, qb)) <: AbstractArray && length(joint) > 0
                    XX = FiniteDiff.finite_difference_jacobian(w -> -transpose(∂g∂ʳposb(joint, w, qa, xb, qb)[1:3, 1:3]) * constraintmat(joint)' * λ, xa)
                    XQ = FiniteDiff.finite_difference_jacobian(w -> -transpose(∂g∂ʳposb(joint, xa, UnitQuaternion(w..., false), xb, qb)[1:3, 1:3]) * constraintmat(joint)' * λ, [qa.w; qa.x; qa.y; qa.z])
                    QX = FiniteDiff.finite_difference_jacobian(w -> -transpose(∂g∂ʳposb(joint, w, qa, xb, qb)[1:3, 4:6]) * constraintmat(joint)' * λ, xa)
                    QQ = FiniteDiff.finite_difference_jacobian(w -> -transpose(∂g∂ʳposb(joint, xa, UnitQuaternion(w..., false), xb, qb)[1:3, 4:6]) * constraintmat(joint)' * λ, [qa.w; qa.x; qa.y; qa.z])

                    Aba[4:6,1:3] = XX
                    Aba[4:6,7:10] = XQ
                    Aba[11:13,1:3] = QX
                    Aba[11:13,7:10] = QQ
                end

                # XX, XQ, QX, QQ = ∂2g∂posbb(joint, posargs2(pstate)..., posargs2(cstate)...)

                if typeof(∂g∂ʳposb(joint, xa, qa, xb, qb)) <: AbstractArray && length(joint) > 0
                    XX = FiniteDiff.finite_difference_jacobian(w -> -transpose(∂g∂ʳposb(joint, xa, qa, w, qb)[1:3, 1:3]) * constraintmat(joint)' * λ, xb)
                    XQ = FiniteDiff.finite_difference_jacobian(w -> -transpose(∂g∂ʳposb(joint, xa, qa, xb, UnitQuaternion(w..., false))[1:3, 1:3]) * constraintmat(joint)' * λ, [qb.w; qb.x; qb.y; qb.z])
                    QX = FiniteDiff.finite_difference_jacobian(w -> -transpose(∂g∂ʳposb(joint, xa, qa, w, qb)[1:3, 4:6]) * constraintmat(joint)' * λ, xb)
                    QQ = FiniteDiff.finite_difference_jacobian(w -> -transpose(∂g∂ʳposb(joint, xa, qa, xb, UnitQuaternion(w..., false))[1:3, 4:6]) * constraintmat(joint)' * λ, [qb.w; qb.x; qb.y; qb.z])

                    Abb[4:6,1:3] = XX
                    Abb[4:6,7:10] = XQ
                    Abb[11:13,1:3] = QX
                    Abb[11:13,7:10] = QQ
                end

                FfzG[pcol13,pcol13] += Aaa
                FfzG[pcol13,ccol13] += Aab
                FfzG[ccol13,pcol13] += Aba
                FfzG[ccol13,ccol13] += Abb
            end
        else
            for (i, childid) in enumerate(eqc.childids)
                childind = childid - Ne
                cbody = getbody(mechanism, childid)
                cstate = cbody.state
                joint = eqc.constraints[i]
                ccol13 = offsetrange(childind,13)

                n1 = 1
                n2 = 0
                for i=1:i-1
                    n1 += length(eqc.constraints[i])
                    n2 += length(eqc.constraints[i])
                end
                n2 += length(joint)
                λ = eqc.λsol[2][n1:n2]

                Abb = zeros(T,13,13)

                xb, qb = posargs2(cstate)

                if typeof(∂g∂ʳposb(joint, xb, qb)) <: AbstractArray && length(joint) > 0

                    XX = FiniteDiff.finite_difference_jacobian(w -> -transpose(∂g∂ʳposb(joint, w, qb)[1:3, 1:3]) * constraintmat(joint)' * λ, xb)
                    XQ = FiniteDiff.finite_difference_jacobian(w -> -transpose(∂g∂ʳposb(joint, xb, UnitQuaternion(w..., false))[1:3, 1:3]) * constraintmat(joint)' * λ, [qb.w; qb.x; qb.y; qb.z])
                    QX = FiniteDiff.finite_difference_jacobian(w -> -transpose(∂g∂ʳposb(joint, w, qb)[1:3, 4:6]) * constraintmat(joint)' * λ, xb)
                    QQ = FiniteDiff.finite_difference_jacobian(w -> -transpose(∂g∂ʳposb(joint, xb, UnitQuaternion(w..., false))[1:3, 4:6]) * constraintmat(joint)' * λ, [qb.w; qb.x; qb.y; qb.z])

                    Abb[4:6,1:3] = XX
                    Abb[4:6,7:10] = XQ
                    Abb[11:13,1:3] = QX
                    Abb[11:13,7:10] = QQ
                end

                FfzG[ccol13,ccol13] += Abb
            end
        end
    end

    return FfzG
end

function spring_damper_datamat(mechanism::Mechanism{T,Nn,Ne,Nb}) where {T,Nn,Ne,Nb}
    Δt = mechanism.Δt
    FfzG = zeros(T,13Nb,12Nb)

    for eqc in mechanism.eqconstraints
        parentid = eqc.parentid

        if parentid !== nothing
            parentind = parentid - Ne
            pbody = getbody(mechanism, parentid)
            pstate = pbody.state
            pcol12 = offsetrange(parentind,12)
            pcol13 = offsetrange(parentind,13)

            for (i, childid) in enumerate(eqc.childids)
                childind = childid - Ne
                cbody = getbody(mechanism, childid)
                cstate = cbody.state
                joint = eqc.constraints[i]
                ccol12 = offsetrange(childind,12)
                ccol13 = offsetrange(childind,13)

                n1 = 1
                n2 = 0
                for j=1:i-1
                    n1 += length(eqc.constraints[j])
                    n2 += length(eqc.constraints[j])
                end
                n2 += length(joint)
                λ = eqc.λsol[2][n1:n2]

                Aaa = zeros(T,13,12)
                Aab = zeros(T,13,12)
                Aba = zeros(T,13,12)
                Abb = zeros(T,13,12)

                Aaa[[4:6; 11:13], [1:3; 7:9]] -= ∂springforcea∂posa(joint, pbody, cbody, Δt)
                Aaa[[4:6; 11:13], [1:3; 7:9]] -= ∂damperforcea∂posa(joint, pbody, cbody, Δt)

                Aab[[4:6; 11:13], [1:3; 7:9]] -= ∂springforcea∂posb(joint, pbody, cbody, Δt)
                Aab[[4:6; 11:13], [1:3; 7:9]] -= ∂damperforcea∂posb(joint, pbody, cbody, Δt)

                Aba[[4:6; 11:13], [1:3; 7:9]] -= ∂springforceb∂posa(joint, pbody, cbody, Δt)
                Aba[[4:6; 11:13], [1:3; 7:9]] -= ∂damperforceb∂posa(joint, pbody, cbody, Δt)

                Abb[[4:6; 11:13], [1:3; 7:9]] -= ∂springforceb∂posb(joint, pbody, cbody, Δt)
                Abb[[4:6; 11:13], [1:3; 7:9]] -= ∂damperforceb∂posb(joint, pbody, cbody, Δt)

                FfzG[pcol13,pcol12] += Aaa
                FfzG[pcol13,ccol12] += Aab
                FfzG[ccol13,pcol12] += Aba
                FfzG[ccol13,ccol12] += Abb
            end
        else
            for (i, childid) in enumerate(eqc.childids)
                pbody = mechanism.origin
                childind = childid - Ne
                cbody = getbody(mechanism, childid)
                cstate = cbody.state
                joint = eqc.constraints[i]
                ccol12 = offsetrange(childind,12)
                ccol13 = offsetrange(childind,13)

                n1 = 1
                n2 = 0
                for i=1:i-1
                    n1 += length(eqc.constraints[i])
                    n2 += length(eqc.constraints[i])
                end
                n2 += length(joint)
                λ = eqc.λsol[2][n1:n2]

                Abb = zeros(T,13,12)

                Abb[[4:6; 11:13], [1:3; 7:9]] -= ∂springforceb∂posb(joint, pbody, cbody, Δt)
                Abb[[4:6; 11:13], [1:3; 7:9]] -= ∂damperforceb∂posb(joint, pbody, cbody, Δt)

                FfzG[ccol13,ccol12] += Abb
            end
        end
    end

    return FfzG
end

function control_datamat(mechanism::Mechanism{T,Nn,Ne,Nb}) where {T,Nn,Ne,Nb}
    Fzu = zeros(T,Nb*13,Nb*12)
    Δt = mechanism.Δt
    for eqc in mechanism.eqconstraints
        parentid = eqc.parentid
        for (i,childid) in enumerate(eqc.childids)
            childind = childid - Ne
            joint = eqc.constraints[i]
            if parentid !== nothing
                parentind = parentid - Ne
                pbody = getbody(mechanism, parentid)
                cbody = getbody(mechanism, childid)

                FaXa, FaQa, τaXa, τaQa, FbXa, FbQa, τbXa, τbQa = ∂Fτ∂posa(joint, pbody.state, cbody.state, mechanism.Δt)
                FaXb, FaQb, τaXb, τaQb, FbXb, FbQb, τbXb, τbQb = ∂Fτ∂posb(joint, pbody.state, cbody.state, mechanism.Δt)


                xa, qa = posargs2(pbody.state)
                Ma = [I zeros(3,3); zeros(4,3) LVᵀmat(qa)]

                xb, qb = posargs2(cbody.state)
                Mb = [I zeros(3,3); zeros(4,3) LVᵀmat(qb)]

                cola6 = offsetrange(parentind,6)
                colb6 = offsetrange(childind,6)
                rowav = offsetrange(parentind,3,13,2)
                rowaω = offsetrange(parentind,3,13,4).+1
                rowbv = offsetrange(childind,3,13,2)
                rowbω = offsetrange(childind,3,13,4).+1

                Fzu[[rowav; rowaω],cola6] = [FaXa FaQa; τaXa τaQa] * Ma

                Fzu[[rowbv; rowbω],cola6] = [FbXa FbQa; τbXa τbQa] * Ma

                Fzu[[rowav; rowaω],colb6] = [FaXb FaQb; τaXb τaQb] * Mb

                Fzu[[rowbv; rowbω],colb6] = [FbXb FbQb; τbXb τbQb] * Mb

            else
                pbody = mechanism.origin
                cbody = getbody(mechanism, childid)
                FaXb, FaQb, τaXb, τaQb, FbXb, FbQb, τbXb, τbQb = ∂Fτ∂posb(joint, cbody.state, mechanism.Δt)

                colb6 = offsetrange(childind,6)
                rowbv = offsetrange(childind,3,13,2)
                rowbω = offsetrange(childind,3,13,4).+1

                x2, q2 = posargs2(cbody.state)
                M = [I zeros(3,3); zeros(4,3) LVᵀmat(q2)]

                Fzu[[rowbv; rowbω], colb6] = [FbXb FbQb; τbXb τbQb] * M

            end
        end
    end
    return Fzu
end

function data_lineardynamics(mechanism::Mechanism{T,Nn,Ne,Nb}, eqcids) where {T,Nn,Ne,Nb}
    Δt = mechanism.Δt
    bodies = mechanism.bodies
    eqcs = mechanism.eqconstraints

    nc = 0
    for eqc in eqcs
        nc += length(eqc)
    end

    nu = 0
    for id in eqcids
        eqc = geteqconstraint(mechanism, id)
        nu += controldim(eqc)
    end



    # get state linearization
    Fz = zeros(T,Nb*13,Nb*12)
    Fu = zeros(T,Nb*13,Nb*6)
    Fλ = zeros(T,Nb*13,nc)
    # Ffz = zeros(T,Nb*13,Nb*13)

    Bcontrol = zeros(T,Nb*6,nu)

    for (ind,body) in enumerate(bodies)
        col6 = offsetrange(ind,6)
        col12 = offsetrange(ind,12)
        col13 = offsetrange(ind,13)

        Fzi = ∂F∂z(body, Δt)
        Fui = ∂F∂u(body, Δt)
        # Ffzi, invFfzquati = ∂F∂fz(body, Δt)

        Fz[col13,col12] = Fzi
        Fu[col13,col6] = Fui
        # Ffz[col13,col13] = Ffzi
    end

    Fz += control_datamat(mechanism)

    n1 = 1
    n2 = 0

    for id in eqcids
        eqc = geteqconstraint(mechanism, id)
        n2 += controldim(eqc)

        parentid = eqc.parentid
        if parentid !== nothing
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

function full_data_matrix(mechanism::Mechanism{T,Nn,Ne,Nb}) where {T,Nn,Ne,Nb}
    mechanism = deepcopy(mechanism)
    system = mechanism.system
    bodies = mechanism.bodies
    eqcs = mechanism.eqconstraints
    ineqcs = mechanism.ineqconstraints
    eqcids = getfield.(eqcs, :id)
    nbodies = length(bodies)

    resdims = [length(system.vector_entries[i].value) for i=1:Nn]
    eqcdims = length.(eqcs)
    ineqcdims = length.(ineqcs)
    bodydims = 6 * ones(Int, nbodies)
    Fz, Fu = data_lineardynamics(mechanism, eqcids)
    data = getdata(mechanism)

    A = zeros(sum(resdims), datadim(mechanism))
    A[1:sum(eqcdims), 1:12Nb] += joint_datamat(mechanism)
    A[sum(eqcdims) .+ (1:sum(bodydims)), 1:12Nb] += Fz[Fz_indices(length(bodies)),:]
    A[sum(eqcdims) .+ (1:sum(bodydims)), 1:12Nb] += joint_jacobian_datamat(mechanism)[Fz_indices(length(bodies)), :] * attitudejacobian(data, Nb)[1:13Nb,1:12Nb]
    A[sum(eqcdims) .+ (1:sum(bodydims)), 1:12Nb] += spring_damper_datamat(mechanism)[Fz_indices(length(bodies)), :]

    offr = 0
    offc = 0
    for body in collect(bodies)
        for ineqc in ineqcs

            if ineqc.parentid == body.id
                bnd = ineqc.constraints[1]
                bnd_type = typeof(bnd)
                Δt = mechanism.Δt
                x3, q3 = posargs3(body.state, Δt)
                x2, v2, q2, ω2 = fullargssol(body.state)
                M = [I zeros(3,3); zeros(4,3) Rmat(ωbar(ω2, Δt)*Δt/2)*LVᵀmat(q2)]
                G = [I zeros(3,3); zeros(4,3) LVᵀmat(q2)]

                function d(vars)
                    x = vars[1:3]
                    q = UnitQuaternion(vars[4:7]..., false)
                    return ∂g∂ʳpos(bnd, x, q)' * ineqc.γsol[2]
                end

                if bnd_type <: ContactBound
                    A[sum(eqcdims) + offr .+ (1:6), offc .+ [1:3; 7:9]] -= _dN(x3, [q3.w, q3.x, q3.y, q3.z], ineqc.γsol[2][1:1], bnd.p) * M
                    A[sum(eqcdims) + offr .+ (1:6), offc .+ [1:3; 7:9]] -= _dB(x3, [q3.w, q3.x, q3.y, q3.z], ineqc.γsol[2][2:4], bnd.p) * M
                elseif bnd_type <: LinearContactBound
                    A[sum(eqcdims) + offr .+ (1:6), offc .+ [1:3; 7:9]] -= FiniteDiff.finite_difference_jacobian(d, [x3; q3.w; q3.x; q3.y; q3.z]) * M
                elseif bnd_type <: ImpactBound
                    A[sum(eqcdims) + offr .+ (1:6), offc .+ [1:3; 7:9]] -= FiniteDiff.finite_difference_jacobian(d, [x3; q3.w; q3.x; q3.y; q3.z]) * M
                end
            end
        end
        offr += 6
        offc += 12
    end

    nu = isempty(eqcs) ? 0 : sum(controldim.(eqcs))
    A[sum(eqcdims) .+ (1:6Nb), 12Nb .+ (1:nu)] += Fu[Fz_indices(Nb), :]

    offr = 0
    for ineqc in ineqcs
        bnd = ineqc.constraints[1]
        body = getbody(mechanism, ineqc.parentid)
        Δt = mechanism.Δt
        N½ = Int(length(ineqc)/2)
        x2, v2, q2, ω2 = fullargssol(body.state)
        x2, q2 = posargs2(body.state)
        x3, q3 = posargs3(body.state, Δt)
        ibody = findfirst(x -> x == body.id, mechanism.bodies.keys)
        bnd_type = typeof(bnd)
        if bnd_type <: ContactBound
            A[sum(eqcdims) + sum(bodydims) + offr + N½ .+ (1:1), (ibody-1)*12 .+ (1:3)] = bnd.ainv3
            A[sum(eqcdims) + sum(bodydims) + offr + N½ .+ (1:1), (ibody-1)*12 .+ (7:9)] = bnd.ainv3 * (VLmat(q3) * Lmat(UnitQuaternion(bnd.p)) * Tmat() + VRᵀmat(q3) * Rmat(UnitQuaternion(bnd.p))) * Rmat(ωbar(ω2, Δt)*Δt/2)*LVᵀmat(q2)
            A[sum(eqcdims) + sum(bodydims) + offr + N½ .+ (3:4), (ibody-1)*12 .+ (7:9)] = dBω([q3.w, q3.x, q3.y, q3.z], ω2, bnd.p) * Rmat(ωbar(ω2, Δt)*Δt/2)*LVᵀmat(q2) # ∇q3B 2x3
        elseif bnd_type <: LinearContactBound
            function d(vars)
                # transforms the velocities of the origin of the link into velocities along all 4 axes of the friction pyramid
                q = UnitQuaternion(vars..., false)
                Bxmat = bnd.Bx
                Bqmat = Bxmat * ∂vrotate∂q(bnd.p, q) * LVᵀmat(q)
                return Bqmat * ω2
            end
            A[sum(eqcdims) + sum(bodydims) + offr + N½ .+ (1:1), (ibody-1)*12 .+ (1:3)] = bnd.ainv3
            A[sum(eqcdims) + sum(bodydims) + offr + N½ .+ (1:1), (ibody-1)*12 .+ (7:9)] = bnd.ainv3 * ∂vrotate∂q(bnd.p,q3) * Rmat(ωbar(ω2, Δt)*Δt/2)*LVᵀmat(q2)
            A[sum(eqcdims) + sum(bodydims) + offr + N½ .+ (3:6), (ibody-1)*12 .+ (7:9)] = FiniteDiff.finite_difference_jacobian(d, [q3.w, q3.x, q3.y, q3.z]) * Rmat(ωbar(ω2, Δt)*Δt/2)*LVᵀmat(q2)
        elseif bnd_type <: ImpactBound
            A[sum(eqcdims) + sum(bodydims) + offr + N½ .+ (1:1), (ibody-1)*12 .+ (1:3)] = bnd.ainv3
            A[sum(eqcdims) + sum(bodydims) + offr + N½ .+ (1:1), (ibody-1)*12 .+ (7:9)] = bnd.ainv3 * (VLmat(q3) * Lmat(UnitQuaternion(bnd.p)) * Tmat() + VRᵀmat(q3) * Rmat(UnitQuaternion(bnd.p))) * LVᵀmat(q2)
        end
        offr += length(ineqc)
    end
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
    attjac = zeros(0,0)
    for i = 1:Nb
        x2, v1, q2, ω1 = unpackdata(data[13*(i-1) .+ (1:13)])
        attjac = cat(attjac, I(6), G(q2), I(3), dims = (1,2))
    end
    ndata = length(data)
    nu = ndata - size(attjac)[1]
    attjac = cat(attjac, I(nu), dims = (1,2))
    return attjac
end

function attitudejacobian_chain(data::AbstractVector, Δt, Nb::Int)
    attjac = zeros(0,0)
    for i = 1:Nb
        x2, v1, q2, ω1 = unpackdata(data[13*(i-1) .+ (1:13)])
        attjac = cat(attjac, I(6), Rmat(ωbar(ω1, Δt)*Δt/2)*LVᵀmat(UnitQuaternion(q2...)), I(3), dims = (1,2))
    end
    ndata = length(data)
    nu = ndata - size(attjac)[1]
    attjac = cat(attjac, I(nu), dims = (1,2))
    return attjac
end

function sensitivities(mech, sol, data)
    setdata!(mech, data)
    setsolution!(mech, sol)
    setentries!(mech)
    datamat = full_data_matrix(mech)
    solmat = full_matrix(mech.system)
    sensi = -1.0 * (solmat \ datamat)
    return sensi
end

function jvp(mech, sol, data, v)
    sensi = sensitivities(mech, sol, data)
    return sensi * v
end

################################################################################
# Finite Difference
################################################################################

function unpackdata(data::AbstractVector)
    x2 = data[1:3]
    v1 = data[4:6]
    q2 = data[7:10]
    ω1 = data[11:13]
    return x2, v1, q2, ω1
end

function setdata!(mechanism::Mechanism, data::AbstractVector)
    off = 0
    for body in mechanism.bodies
        x2, v1, q2, ω1 = unpackdata(data[off+1:end]); off += 13
        body.state.x1 = x2
        body.state.v15 = v1
        body.state.q1 = UnitQuaternion(q2...)
        body.state.ϕ15 = ω1
        body.state.x2[1] = x2
        body.state.q2[1] = UnitQuaternion(q2...)
        setsolution!(body)
        body.state.F2[1] = SVector{3}([0,0,0.])
        body.state.τ2[1] = SVector{3}([0,0,0.])
    end
    for eqc in mechanism.eqconstraints
        dim = controldim(eqc)
        if dim > 0
            u = data[off .+ (1:dim)]; off += dim
            setForce!(mechanism, eqc, u)
        end
    end
    foreach(applyFτ!, mechanism.eqconstraints, mechanism, false)
    return nothing
end

function getdata(mechanism::Mechanism{T}) where T
    data = Vector{T}()
    for body in mechanism.bodies
        x2 = body.state.x1
        v1 = body.state.v15
        qc = body.state.q1
        q2 = [qc.w, qc.x, qc.y, qc.z]
        ω1 = body.state.ϕ15
        push!(data, [x2; v1; q2; ω1]...)
    end
    for eqc in mechanism.eqconstraints
        if controldim(eqc) > 0
            tra = eqc.constraints[findfirst(x -> typeof(x) <: Translational, eqc.constraints)]
            rot = eqc.constraints[findfirst(x -> typeof(x) <: Rotational, eqc.constraints)]
            F = tra.Fτ
            τ = rot.Fτ
            u = [nullspacemat(tra) * F; nullspacemat(rot) * τ]
            push!(data, u...)
        end
    end
    return data
end

function setsolution!(mechanism::Mechanism{T}, sol::AbstractVector) where T
    off = 0
    for (i,eqc) in enumerate(mechanism.eqconstraints)
        nλ = length(eqc)
        λ = sol[off .+ (1:nλ)]; off += nλ
        eqc.λsol[2] = λ
    end
    for (i,body) in enumerate(mechanism.bodies)
        nv = 3
        nω = 3
        v2 = sol[off .+ (1:nv)]; off += nv
        ω2 = sol[off .+ (1:nω)]; off += nω
        body.state.vsol[2] = v2
        body.state.ϕsol[2] = ω2
    end
    for (i,ineqc) in enumerate(mechanism.ineqconstraints)
        N = length(ineqc)
        N½ = Int(N/2)
        s = sol[off .+ (1:N½)]; off += N½
        γ = sol[off .+ (1:N½)]; off += N½
        ineqc.ssol[2] = s
        ineqc.γsol[2] = γ
    end
    return nothing
end

function getsolution(mechanism::Mechanism{T}) where T
    sol = T[]
    for (i,eqc) in enumerate(mechanism.eqconstraints)
        λ = eqc.λsol[2]
        push!(sol, λ...)
    end
    for (i,body) in enumerate(mechanism.bodies)
        v2 = body.state.vsol[2]
        ω2 = body.state.ϕsol[2]
        push!(sol, [v2; ω2]...)
    end
    for (i,ineqc) in enumerate(mechanism.ineqconstraints)
        s = ineqc.ssol[2]
        γ = ineqc.γsol[2]
        push!(sol, [s; γ]...)
    end
    return sol
end

function evaluate_solution!(mechanism::Mechanism, data::AbstractVector; ϵr=1e-8, ϵb=1.0e-8)
    setdata!(mechanism, data)
    mehrotra!(mechanism, opts = InteriorPointOptions(rtol = ϵr, btol = ϵb, undercut=1.2, verbose = true))
    sol = getsolution(mechanism)
    return sol
end

function evaluate_residual!(mechanism::Mechanism, data::AbstractVector, sol::AbstractVector)
    system = mechanism.system
    setdata!(mechanism, data)
    setsolution!(mechanism, sol)
    setentries!(mechanism)
    return full_vector(system)
end

function finitediff_sensitivity(mechanism::Mechanism, data::AbstractVector; ϵr = 1e-8, ϵb=1.0e-8, δ = 1e-5, verbose = false)
    ndata = datadim(mechanism, quat = true)
    nsol = soldim(mechanism)
    jac = zeros(nsol, ndata)

    for i = 1:ndata
        verbose && println("$i / $ndata")
        datap = deepcopy(data)
        datam = deepcopy(data)
        datap[i] += δ
        datam[i] -= δ

        mechanismp = deepcopy(mechanism)
        setdata!(mechanismp, deepcopy(datap))
        mehrotra!(mechanismp, opts = InteriorPointOptions(rtol = ϵr, btol = ϵb, undercut = 1.2, verbose = false))
        solp = getsolution(mechanismp)

        mechanismm = deepcopy(mechanism)
        setdata!(mechanismm, deepcopy(datam))
        mehrotra!(mechanismm, opts = InteriorPointOptions(rtol = ϵr, btol = ϵb, undercut = 1.2, verbose = false))
        solm = getsolution(mechanismm)

        jac[:,i] = (solp - solm) / (2δ)
    end
    return jac
end

function finitediff_data_matrix(mechanism::Mechanism, data::AbstractVector,
        sol::AbstractVector; δ = 1e-8, verbose = false)
    nsol = soldim(mechanism)
    ndata = datadim(mechanism, quat = true)
    jac = zeros(nsol, ndata)

    setdata!(mechanism, deepcopy(data))
    setsolution!(mechanism, deepcopy(sol))

    for i = 1:ndata
        verbose && println("$i / $ndata")
        datap = deepcopy(data)
        datam = deepcopy(data)
        datap[i] += δ
        datam[i] -= δ
        rp = evaluate_residual!(deepcopy(mechanism), datap, deepcopy(sol))
        rm = evaluate_residual!(deepcopy(mechanism), datam, deepcopy(sol))
        jac[:,i] = (rp - rm) / (2δ)
    end
    return jac
end

function finitediff_sol_matrix(mechanism::Mechanism, data::AbstractVector,
        sol::AbstractVector; δ = 1e-8, verbose = false)
    nsol = soldim(mechanism)
    jac = zeros(nsol, nsol)

    setdata!(mechanism, data)
    setsolution!(mechanism, sol)

    for i = 1:nsol
        verbose && println("$i / $nsol")
        solp = deepcopy(sol)
        solm = deepcopy(sol)
        solp[i] += δ
        solm[i] -= δ
        rp = evaluate_residual!(deepcopy(mechanism), data, solp)
        rm = evaluate_residual!(deepcopy(mechanism), data, solm)
        jac[:,i] = (rp - rm) / (2δ)
    end
    return jac
end


################################################################################
# Index and Dimensions
################################################################################

function Fz_indices(Nb::Int)
    return vcat([13*(i-1) .+ [4:6; 11:13] for i = 1:Nb]...)
end

function datadim(mechanism::Mechanism; quat::Bool = false)
    d = 0
    bodies = mechanism.bodies
    for body in bodies
        d += 6 * 2
        quat && (d += 1)
    end
    eqcs = mechanism.eqconstraints
    for eqc in eqcs
        d += controldim(eqc)
    end
    return d
end

function soldim(mechanism::Mechanism)
    d = 0
    d += 6 * length(mechanism.bodies)
    d += sum(length.(mechanism.eqconstraints))
    !isempty(mechanism.ineqconstraints) && (d += sum(length.(mechanism.ineqconstraints)))
    return d
end

function controldim(eqc::EqualityConstraint{T,N,Nc,Cs}; ignore_floating_base::Bool = false) where {T,N,Nc,Cs}
    ignore_floating_base && (N == 0) && return 0
    return 6 - N
end


mech = getmechanism(:hopper)
eqc1 = collect(mech.eqconstraints)[1]
eqc2 = collect(mech.eqconstraints)[2]
length(eqc1)
length(eqc1.constraints[1])
length(eqc1.constraints[2])
length(eqc2)
length(eqc2.constraints[1])
length(eqc2.constraints[2])


controldim(eqc1, ignore_floating_base = false)
controldim(eqc2, ignore_floating_base = false)
