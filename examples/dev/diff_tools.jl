function linearconstraints2(mechanism::Mechanism{T,Nn,Ne,Nb}) where {T,Nn,Ne,Nb}
    Δt = mechanism.Δt
    eqcs = mechanism.eqconstraints

    nc = 0
    for eqc in eqcs
        nc += length(eqc)
    end

    Gl = zeros(T,nc,Nb*12)

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

                ind2 += length(eqc.constraints[i])
                range = oneindc+ind1:oneindc+ind2

                pcol3a12 = offsetrange(parentind,3,12,1)
                pcol3b12 = offsetrange(parentind,3,12,2)
                pcol3c12 = offsetrange(parentind,3,12,3)
                pcol3d12 = offsetrange(parentind,3,12,4)

                ccol3a12 = offsetrange(childind,3,12,1)
                ccol3b12 = offsetrange(childind,3,12,2)
                ccol3c12 = offsetrange(childind,3,12,3)
                ccol3d12 = offsetrange(childind,3,12,4)

                # pXl, pQl = ∂g∂posa(eqc.constraints[i], posargsnext(pstate, Δt)..., posargsnext(cstate, Δt)...) # x3
                # cXl, cQl = ∂g∂posb(eqc.constraints[i], posargsnext(pstate, Δt)..., posargsnext(cstate, Δt)...) # x3

                pXl, pQl = ∂g∂posa(eqc.constraints[i], pbody, cbody, Δt) # x3
                cXl, cQl = ∂g∂posb(eqc.constraints[i], pbody, cbody, Δt) # x3

                mat = constraintmat(eqc.constraints[i])
                # @show mat
                pGlx = mat * pXl
                pGlq = mat * pQl
                cGlx = mat * cXl
                cGlq = mat * cQl
                # @show pGlx
                # @show pGlq
                # @show cGlx
                # @show cGlq

                Gl[range,pcol3a12] = pGlx
                Gl[range,pcol3c12] = pGlq*Rmat(ωbar(pstate.ωc, Δt)*Δt/2)*LVᵀmat(pstate.qc)

                Gl[range,ccol3a12] = cGlx
                Gl[range,ccol3c12] = cGlq*Rmat(ωbar(cstate.ωc, Δt)*Δt/2)*LVᵀmat(cstate.qc)

                # if typeof(eqc.constraints[i]) <: Torque
                #     pXl1, pQl1 = ∂g∂posa1(eqc.constraints[i], pbody, cbody, Δt) # x3
                #     cXl1, cQl1 = ∂g∂posb1(eqc.constraints[i], pbody, cbody, Δt) # x3
                #     @show typeof(pQl1)
                #     @show typeof(cQl1)
                #     pGlq1 = mat * pQl1
                #     cGlq1 = mat * cQl1
                #     Gl[range,pcol3c12] += pGlq1 * LVᵀmat(pstate.qc)
                #     Gl[range,ccol3c12] += cGlq1 * LVᵀmat(cstate.qc)
                # end

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


                # cXl, cQl =  ∂g∂posb(eqc.constraints[i], posargsnext(cstate, Δt)...) # x3
                cXl, cQl =  ∂g∂posb(eqc.constraints[i], mechanism.origin, cbody, Δt) # x3
                # cXl, cQl =  ∂g∂posb(eqc.constraints[i], cstate, Δt) # x3

                mat = constraintmat(eqc.constraints[i])
                cGlx = mat * cXl
                cGlq = mat * cQl

                Gl[range,ccol3a12] = cGlx
                # @show range
                # @show ccol3c12
                # @show size(cGlq)
                # @show size(cGlq*Rmat(ωbar(cstate.ωc, Δt)*Δt/2)*LVᵀmat(cstate.qc))
                Gl[range,ccol3c12] = cGlq*Rmat(ωbar(cstate.ωc, Δt)*Δt/2)*LVᵀmat(cstate.qc)

                # if typeof(eqc.constraints[i]) <: Torque
                #     cXl1, cQl1 = ∂g∂posb1(eqc.constraints[i], mechanism.origin, cbody, Δt)
                #     cQlq1 = mat * cQl1

                #     # @show range
                #     # @show ccol3c12
                #     # @show cQl1
                #     # @show range
                #     # @show ccol3c12
                #     # @show size(cQl1)
                #     # @show size(cQlq1)

                #     Gl[range,ccol3c12] += cQlq1 * LVᵀmat(cstate.qc)
                # end
                ind1 = ind2+1
            end
        end

        oneindc += length(eqc)
    end

    return Gl
end

function linearconstraintmapping2(mechanism::Mechanism{T,Nn,Ne,Nb}) where {T,Nn,Ne,Nb}
    FfzG = zeros(T,Nb*13,Nb*13)

    K = zeros(T,9,9)
    K[1,1] = K[2,4] = K[3,7] = K[4,2] = K[5,5] = K[6,8] = K[7,3] = K[8,6] = K[9,9] = 1
    E = SMatrix{3,3,T,9}(I)

    for eqc in mechanism.eqconstraints
        parentid = eqc.parentid

        if parentid !== nothing
            parentind = parentid - Ne
            body1 = getbody(mechanism, parentid)
            state1 = body1.state
            pcol13 = offsetrange(parentind,13)

            for (i, childid) in enumerate(eqc.childids)
                childind = childid - Ne
                body2 = getbody(mechanism, childid)
                state2 = body2.state
                constraint = eqc.constraints[i]
                ccol13 = offsetrange(childind,13)

                n1 = 1
                n2 = 0
                for j=1:i-1
                    n1 += length(eqc.constraints[j])
                    n2 += length(eqc.constraints[j])
                end
                n2 += length(eqc.constraints[i])
                λ = eqc.λsol[2][n1:n2]

                Aaa = zeros(T,13,13)
                Aab = zeros(T,13,13)
                Aba = zeros(T,13,13)
                Abb = zeros(T,13,13)

                kronproduct = -kron(λ'*Array(constraintmat(constraint)),E)*K

                XX, XQ, QX, QQ = ∂2g∂posaa(constraint, state1.xc, state1.qc, state2.xc, state2.qc)
                Aaa[4:6,1:3] = kronproduct*XX
                Aaa[4:6,7:10] = kronproduct*XQ
                Aaa[11:13,1:3] = kronproduct*QX
                Aaa[11:13,7:10] = kronproduct*QQ

                XX, XQ, QX, QQ = ∂2g∂posab(constraint, state1.xc, state1.qc, state2.xc, state2.qc)
                Aab[4:6,1:3] = kronproduct*XX
                Aab[4:6,7:10] = kronproduct*XQ
                Aab[11:13,1:3] = kronproduct*QX
                Aab[11:13,7:10] = kronproduct*QQ

                XX, XQ, QX, QQ = ∂2g∂posba(constraint, state1.xc, state1.qc, state2.xc, state2.qc)
                Aba[4:6,1:3] = kronproduct*XX
                Aba[4:6,7:10] = kronproduct*XQ
                Aba[11:13,1:3] = kronproduct*QX
                Aba[11:13,7:10] = kronproduct*QQ

                XX, XQ, QX, QQ = ∂2g∂posbb(constraint, state1.xc, state1.qc, state2.xc, state2.qc)
                Abb[4:6,1:3] = kronproduct*XX
                Abb[4:6,7:10] = kronproduct*XQ
                Abb[11:13,1:3] = kronproduct*QX
                Abb[11:13,7:10] = kronproduct*QQ

                FfzG[pcol13,pcol13] += Aaa
                FfzG[pcol13,ccol13] += Aab
                FfzG[ccol13,pcol13] += Aba
                FfzG[ccol13,ccol13] += Abb
            end
        else
            for (i, childid) in enumerate(eqc.childids)
                childind = childid - Ne
                body2 = getbody(mechanism, childid)
                state2 = body2.state
                constraint = eqc.constraints[i]
                ccol13 = offsetrange(childind,13)

                n1 = 1
                n2 = 0
                for i=1:i-1
                    n1 += length(eqc.constraints[i])
                    n2 += length(eqc.constraints[i])
                end
                n2 += length(eqc.constraints[i])
                λ = eqc.λsol[2][n1:n2]

                Abb = zeros(T,13,13)

                kronproduct = -kron(λ'*Array(constraintmat(constraint)),E)*K

                XX, XQ, QX, QQ = ∂2g∂posbb(constraint, state2.xc, state2.qc)
                Abb[4:6,1:3] = kronproduct*XX
                Abb[4:6,7:10] = kronproduct*XQ
                Abb[11:13,1:3] = kronproduct*QX
                Abb[11:13,7:10] = kronproduct*QQ

                FfzG[ccol13,ccol13] += Abb
            end
        end
    end

    return FfzG
end

function linearconstraintmapping3(mechanism::Mechanism{T,Nn,Ne,Nb}) where {T,Nn,Ne,Nb}
    Δt = mechanism.Δt
    FfzG = zeros(T,Nb*13,Nb*13)

    K = zeros(T,9,9)
    K[1,1] = K[2,4] = K[3,7] = K[4,2] = K[5,5] = K[6,8] = K[7,3] = K[8,6] = K[9,9] = 1
    E = SMatrix{3,3,T,9}(I)

    for eqc in mechanism.eqconstraints
        parentid = eqc.parentid

        if parentid !== nothing
            parentind = parentid - Ne
            body1 = getbody(mechanism, parentid)
            state1 = body1.state
            pcol13 = offsetrange(parentind,13)

            for (i, childid) in enumerate(eqc.childids)
                childind = childid - Ne
                body2 = getbody(mechanism, childid)
                state2 = body2.state
                constraint = eqc.constraints[i]
                ccol13 = offsetrange(childind,13)

                n1 = 1
                n2 = 0
                for j=1:i-1
                    n1 += length(eqc.constraints[j])
                    n2 += length(eqc.constraints[j])
                end
                n2 += length(eqc.constraints[i])
                λ = eqc.λsol[2][n1:n2]

                Aaa = zeros(T,13,13)
                Aab = zeros(T,13,13)
                Aba = zeros(T,13,13)
                Abb = zeros(T,13,13)

                kronproduct = -kron(λ'*Array(constraintmat(constraint)),E)*K

                # XX, XQ, QX, QQ = ∂2g∂posaa(constraint, posargsnext(state1, Δt)..., posargsnext(state2, Δt)...)
                XX, XQ, QX, QQ = ∂2g∂posaa(constraint, posargsk(state1)..., posargsk(state2)...)
                Aaa[4:6,1:3] = kronproduct*XX
                Aaa[4:6,7:10] = kronproduct*XQ
                Aaa[11:13,1:3] = kronproduct*QX
                Aaa[11:13,7:10] = kronproduct*QQ

                # XX, XQ, QX, QQ = ∂2g∂posab(constraint, posargsnext(state1, Δt)..., posargsnext(state2, Δt)...)
                XX, XQ, QX, QQ = ∂2g∂posab(constraint, posargsk(state1)..., posargsk(state2)...)
                Aab[4:6,1:3] = kronproduct*XX
                Aab[4:6,7:10] = kronproduct*XQ
                Aab[11:13,1:3] = kronproduct*QX
                Aab[11:13,7:10] = kronproduct*QQ

                # XX, XQ, QX, QQ = ∂2g∂posba(constraint, posargsnext(state1, Δt)..., posargsnext(state2, Δt)...)
                XX, XQ, QX, QQ = ∂2g∂posba(constraint, posargsk(state1)..., posargsk(state2)...)
                Aba[4:6,1:3] = kronproduct*XX
                Aba[4:6,7:10] = kronproduct*XQ
                Aba[11:13,1:3] = kronproduct*QX
                Aba[11:13,7:10] = kronproduct*QQ

                # XX, XQ, QX, QQ = ∂2g∂posbb(constraint, posargsnext(state1, Δt)..., posargsnext(state2, Δt)...)
                XX, XQ, QX, QQ = ∂2g∂posbb(constraint, posargsk(state1)..., posargsk(state2)...)
                Abb[4:6,1:3] = kronproduct*XX
                Abb[4:6,7:10] = kronproduct*XQ
                Abb[11:13,1:3] = kronproduct*QX
                Abb[11:13,7:10] = kronproduct*QQ

                FfzG[pcol13,pcol13] += Aaa
                FfzG[pcol13,ccol13] += Aab
                FfzG[ccol13,pcol13] += Aba
                FfzG[ccol13,ccol13] += Abb
            end
        else
            for (i, childid) in enumerate(eqc.childids)
                childind = childid - Ne
                body2 = getbody(mechanism, childid)
                state2 = body2.state
                constraint = eqc.constraints[i]
                ccol13 = offsetrange(childind,13)

                n1 = 1
                n2 = 0
                for i=1:i-1
                    n1 += length(eqc.constraints[i])
                    n2 += length(eqc.constraints[i])
                end
                n2 += length(eqc.constraints[i])
                λ = eqc.λsol[2][n1:n2]

                Abb = zeros(T,13,13)

                kronproduct = -kron(λ'*Array(constraintmat(constraint)),E)*K

                # XX, XQ, QX, QQ = ∂2g∂posbb(constraint, posargsnext(state2, Δt)...)
                XX, XQ, QX, QQ = ∂2g∂posbb(constraint, posargsk(state2)...)
                Abb[4:6,1:3] = kronproduct*XX
                Abb[4:6,7:10] = kronproduct*XQ
                Abb[11:13,1:3] = kronproduct*QX
                Abb[11:13,7:10] = kronproduct*QQ

                FfzG[ccol13,ccol13] += Abb
            end
        end
    end

    return FfzG
end

function linearconstraintmapping4(mechanism::Mechanism{T,Nn,Ne,Nb}) where {T,Nn,Ne,Nb}
    Δt = mechanism.Δt
    FfzG = zeros(T,Nb*13,Nb*12)

    K = zeros(T,9,9)
    K[1,1] = K[2,4] = K[3,7] = K[4,2] = K[5,5] = K[6,8] = K[7,3] = K[8,6] = K[9,9] = 1
    E = SMatrix{3,3,T,9}(I)

    for eqc in mechanism.eqconstraints
        parentid = eqc.parentid

        if parentid !== nothing
            parentind = parentid - Ne
            body1 = getbody(mechanism, parentid)
            state1 = body1.state
            pcol12 = offsetrange(parentind,12)
            pcol13 = offsetrange(parentind,13)

            for (i, childid) in enumerate(eqc.childids)
                childind = childid - Ne
                body2 = getbody(mechanism, childid)
                state2 = body2.state
                constraint = eqc.constraints[i]
                ccol12 = offsetrange(childind,12)
                ccol13 = offsetrange(childind,13)

                n1 = 1
                n2 = 0
                for j=1:i-1
                    n1 += length(eqc.constraints[j])
                    n2 += length(eqc.constraints[j])
                end
                n2 += length(eqc.constraints[i])
                λ = eqc.λsol[2][n1:n2]

                Aaa = zeros(T,13,12)
                Aab = zeros(T,13,12)
                Aba = zeros(T,13,12)
                Abb = zeros(T,13,12)

                kronproduct = -kron(λ'*Array(constraintmat(constraint)),E)*K

                # XX, XQ, QX, QQ = ∂2g∂posaa(constraint, posargsnext(state1, Δt)..., posargsnext(state2, Δt)...)
                # Aaa[4:6,1:3] = kronproduct*XX
                # Aaa[4:6,7:9] = kronproduct*XQ
                # Aaa[11:13,1:3] = kronproduct*QX
                # Aaa[11:13,7:9] = kronproduct*QQ
                Aaa[[4:6; 11:13], [1:3; 7:9]] -= ∂springforcea∂posa(constraint, body1, body2, Δt)
                Aaa[[4:6; 11:13], [1:3; 7:9]] -= ∂damperforcea∂posa(constraint, body1, body2, Δt)

                # XX, XQ, QX, QQ = ∂2g∂posab(constraint, posargsnext(state1, Δt)..., posargsnext(state2, Δt)...)
                # Aab[4:6,1:3] = kronproduct*XX
                # Aab[4:6,7:9] = kronproduct*XQ
                # Aab[11:13,1:3] = kronproduct*QX
                # Aab[11:13,7:9] = kronproduct*QQ
                Aab[[4:6; 11:13], [1:3; 7:9]] -= ∂springforcea∂posb(constraint, body1, body2, Δt)
                Aab[[4:6; 11:13], [1:3; 7:9]] -= ∂damperforcea∂posb(constraint, body1, body2, Δt)

                # XX, XQ, QX, QQ = ∂2g∂posba(constraint, posargsnext(state1, Δt)..., posargsnext(state2, Δt)...)
                # Aba[4:6,1:3] = kronproduct*XX
                # Aba[4:6,7:9] = kronproduct*XQ
                # Aba[11:13,1:3] = kronproduct*QX
                # Aba[11:13,7:9] = kronproduct*QQ
                Aba[[4:6; 11:13], [1:3; 7:9]] -= ∂springforceb∂posa(constraint, body1, body2, Δt)
                Aba[[4:6; 11:13], [1:3; 7:9]] -= ∂damperforceb∂posa(constraint, body1, body2, Δt)

                # XX, XQ, QX, QQ = ∂2g∂posbb(constraint, posargsnext(state1, Δt)..., posargsnext(state2, Δt)...)
                # Abb[4:6,1:3] = kronproduct*XX
                # Abb[4:6,7:9] = kronproduct*XQ
                # Abb[11:13,1:3] = kronproduct*QX
                # Abb[11:13,7:9] = kronproduct*QQ
                Abb[[4:6; 11:13], [1:3; 7:9]] -= ∂springforceb∂posb(constraint, body1, body2, Δt)
                Abb[[4:6; 11:13], [1:3; 7:9]] -= ∂damperforceb∂posb(constraint, body1, body2, Δt)

                FfzG[pcol13,pcol12] += Aaa
                FfzG[pcol13,ccol12] += Aab
                FfzG[ccol13,pcol12] += Aba
                FfzG[ccol13,ccol12] += Abb
            end
        else
            for (i, childid) in enumerate(eqc.childids)
                body1 = mechanism.origin
                childind = childid - Ne
                body2 = getbody(mechanism, childid)
                state2 = body2.state
                constraint = eqc.constraints[i]
                ccol12 = offsetrange(childind,12)
                ccol13 = offsetrange(childind,13)

                n1 = 1
                n2 = 0
                for i=1:i-1
                    n1 += length(eqc.constraints[i])
                    n2 += length(eqc.constraints[i])
                end
                n2 += length(eqc.constraints[i])
                λ = eqc.λsol[2][n1:n2]

                Abb = zeros(T,13,12)

                kronproduct = -kron(λ'*Array(constraintmat(constraint)),E)*K

                # XX, XQ, QX, QQ = ∂2g∂posbb(constraint, posargsnext(state2, Δt)...)
                # Abb[4:6,1:3] = kronproduct*XX
                # Abb[4:6,7:10] = kronproduct*XQ
                # Abb[11:13,1:3] = kronproduct*QX
                # Abb[11:13,7:10] = kronproduct*QQ
                Abb[[4:6; 11:13], [1:3; 7:9]] -= ∂springforceb∂posb(constraint, body1, body2, Δt)
                Abb[[4:6; 11:13], [1:3; 7:9]] -= ∂damperforceb∂posb(constraint, body1, body2, Δt)

                FfzG[ccol13,ccol12] += Abb
            end
        end
    end

    return FfzG
end


function linearforcemapping2(mechanism::Mechanism{T,Nn,Ne,Nb}) where {T,Nn,Ne,Nb}
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

                # xa, va, qa, ωa = fullargssol(pbody.state)
                # Ma = [I zeros(3,3); zeros(4,3) Rmat(ωbar(ωa, Δt)*Δt/2)*LVᵀmat(qa)]
                xa, qa = posargsk(pbody.state)
                Ma = [I zeros(3,3); zeros(4,3) LVᵀmat(qa)]

                # xb, vb, qb, ωb = fullargssol(cbody.state)
                # Mb = [I zeros(3,3); zeros(4,3) Rmat(ωbar(ωb, Δt)*Δt/2)*LVᵀmat(qb)]
                xb, qb = posargsk(cbody.state)
                Mb = [I zeros(3,3); zeros(4,3) LVᵀmat(qb)]

                cola6 = offsetrange(parentind,6)
                colb6 = offsetrange(childind,6)
                rowav = offsetrange(parentind,3,13,2)
                rowaω = offsetrange(parentind,3,13,4).+1
                rowbv = offsetrange(childind,3,13,2)
                rowbω = offsetrange(childind,3,13,4).+1

                # Fzu[rowav,cola6] = [FaXa FaQa]
                # Fzu[rowaω,cola6] = [τaXa τaQa]
                Fzu[[rowav; rowaω],cola6] = [FaXa FaQa; τaXa τaQa] * Ma
                # Fzu[[rowav; rowaω],cola6] -= ∂springforcea∂posa(joint, pbody, cbody, Δt)
                # Fzu[[rowav; rowaω],cola6] -= ∂damperforcea∂posa(joint, pbody, cbody, Δt)

                # Fzu[rowbv,cola6] = [FbXa FbQa]
                # Fzu[rowbω,cola6] = [τbXa τbQa]
                Fzu[[rowbv; rowbω],cola6] = [FbXa FbQa; τbXa τbQa] * Ma
                # Fzu[[rowbv; rowbω],cola6] -= ∂springforceb∂posa(joint, pbody, cbody, Δt)
                # Fzu[[rowbv; rowbω],cola6] -= ∂damperforceb∂posa(joint, pbody, cbody, Δt)

                # Fzu[rowav,colb6] = [FaXb FaQb]
                # Fzu[rowaω,colb6] = [τaXb τaQb]
                Fzu[[rowav; rowaω],colb6] = [FaXb FaQb; τaXb τaQb] * Mb
                # Fzu[[rowav; rowaω],colb6] -= ∂springforcea∂posb(joint, pbody, cbody, Δt)
                # Fzu[[rowav; rowaω],colb6] -= ∂damperforcea∂posb(joint, pbody, cbody, Δt)

                # Fzu[rowbv,colb6] = [FbXb FbQb]
                # Fzu[rowbω,colb6] = [τbXb τbQb]
                Fzu[[rowbv; rowbω],colb6] = [FbXb FbQb; τbXb τbQb] * Mb
                # Fzu[[rowbv; rowbω],colb6] -= ∂springforceb∂posb(joint, pbody, cbody, Δt)
                # Fzu[[rowbv; rowbω],colb6] -= ∂damperforceb∂posb(joint, pbody, cbody, Δt)
            else
                pbody = mech.origin
                cbody = getbody(mechanism, childid)
                FaXb, FaQb, τaXb, τaQb, FbXb, FbQb, τbXb, τbQb = ∂Fτ∂posb(joint, cbody.state, mechanism.Δt)

                colb6 = offsetrange(childind,6)
                rowbv = offsetrange(childind,3,13,2)
                rowbω = offsetrange(childind,3,13,4).+1

                # x2, v2, q2, ω2 = fullargssol(cbody.state)
                # M = [I zeros(3,3); zeros(4,3) Rmat(ωbar(ω2, Δt)*Δt/2)*LVᵀmat(q2)]
                x2, q2 = posargsk(cbody.state)
                M = [I zeros(3,3); zeros(4,3) LVᵀmat(q2)]
                # Fzu[rowbv,colb6] = [FbXb FbQb]
                # Fzu[rowbω,colb6] = [τbXb τbQb]
                Fzu[[rowbv; rowbω], colb6] = [FbXb FbQb; τbXb τbQb] * M
                # Fzu[[rowbv; rowbω], colb6] -= ∂springforceb∂posb(joint, pbody, cbody, Δt)
                # Fzu[[rowbv; rowbω], colb6] -= ∂damperforceb∂posb(joint, pbody, cbody, Δt)
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
        nu += getcontroldim(eqc)
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

    # Ffz += linearconstraintmapping(mechanism)
    Fz += linearforcemapping2(mechanism)

    n1 = 1
    n2 = 0

    for id in eqcids
        eqc = geteqconstraint(mechanism, id)
        n2 += getcontroldim(eqc)

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

    G = linearconstraints2(mechanism)
    # G, Fλ = linearconstraints(mechanism)

    # F contains the following
    # F = [x3 - v2 Δt;
    #      m(v2 - v1) + g ...
    #      q3 - Δt/2 L(q2) [...];
    #      Jω2
    #      ] ∈ R13
    return Fz, Fu * Bcontrol, G
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
    eqcdims = getdim.(eqcs)
    ineqcdims = getdim.(ineqcs)
    bodydims = 6 * ones(Int, nbodies)
    Fz, Fu, G = data_lineardynamics(mechanism, eqcids)
    data = getdata(mechanism)

    A = zeros(sum(resdims), datadim(mechanism))
    A[1:sum(eqcdims), 1:12Nb] += linearconstraints2(mechanism)
    A[sum(eqcdims) .+ (1:sum(bodydims)), 1:12Nb] += Fz[Fz_indices(length(bodies)),:]
    # A[sum(eqcdims) .+ (1:sum(bodydims)), 1:12Nb] += linearconstraintmapping3(mechanism)[Fz_indices(length(bodies)), :] * attitudejacobian_chain(data, mechanism.Δt, Nb)[1:13Nb,1:12Nb]
    A[sum(eqcdims) .+ (1:sum(bodydims)), 1:12Nb] += linearconstraintmapping3(mechanism)[Fz_indices(length(bodies)), :] * attitudejacobian(data, Nb)[1:13Nb,1:12Nb]
    A[sum(eqcdims) .+ (1:sum(bodydims)), 1:12Nb] += linearconstraintmapping4(mechanism)[Fz_indices(length(bodies)), :]

    offr = 0
    offc = 0
    for body in collect(bodies)
        for ineqc in ineqcs
            if ineqc.parentid == body.id
                cont = ineqc.constraints[1]
                p = cont.p
                Δt = mechanism.Δt
                x3, q3 = posargsnext(body.state, Δt)
                x2, v2, q2, ω2 = fullargssol(body.state)
                M = [I zeros(3,3); zeros(4,3) Rmat(ωbar(ω2, Δt)*Δt/2)*LVᵀmat(q2)]
                A[sum(eqcdims) + offr .+ (1:6), offc .+ [1:3; 7:9]] -= _dN(x3, [q3.w, q3.x, q3.y, q3.z], ineqc.γsol[2][1:1], p) * M
                A[sum(eqcdims) + offr .+ (1:6), offc .+ [1:3; 7:9]] -= _dB(x3, [q3.w, q3.x, q3.y, q3.z], ineqc.γsol[2][2:4], p) * M
            end
        end
        offr += 6
        offc += 12
    end

    nu = isempty(eqcs) ? 0 : sum(getcontroldim.(eqcs))
    A[sum(eqcdims) .+ (1:6Nb), 12Nb .+ (1:nu)] += Fu[Fz_indices(Nb), :]

    offr = 0
    for ineqc in ineqcs
        cont = ineqc.constraints[1]
        body = getbody(mechanism, ineqc.parentid)
        Δt = mechanism.Δt
        N½ = Int(getdim(ineqc)/2)
        x2, v2, q2, ω2 = fullargssol(body.state)
        x2, q2 = posargsk(body.state)
        x3, q3 = posargsnext(body.state, Δt)
        ibody = findfirst(x -> x == body.id, mechanism.bodies.keys)
        A[sum(eqcdims) + sum(bodydims) + offr + N½ .+ (1:1), (ibody-1)*12 .+ (1:3)] = cont.ainv3
        A[sum(eqcdims) + sum(bodydims) + offr + N½ .+ (1:1), (ibody-1)*12 .+ (7:9)] = cont.ainv3 * (VLmat(q3) * Lmat(UnitQuaternion(cont.p)) * Tmat() + VRᵀmat(q3) * Rmat(UnitQuaternion(cont.p))) * Rmat(ωbar(ω2, Δt)*Δt/2)*LVᵀmat(q2)
        A[sum(eqcdims) + sum(bodydims) + offr + N½ .+ (3:4), (ibody-1)*12 .+ (7:9)] = dBω([q3.w, q3.x, q3.y, q3.z], ω2, cont.p) * Rmat(ωbar(ω2, Δt)*Δt/2)*LVᵀmat(q2) # ∇q3B 2x3
        offr += getdim(ineqc)
    end
    return A
end

function Fz_indices(Nb::Int)
    return vcat([13*(i-1) .+ [4:6; 11:13] for i = 1:Nb]...)
end

function getBcontrol(mechanism::Mechanism{T,Nn,Ne,Nb}) where {T,Nn,Ne,Nb}
    bodies = mechanism.bodies
    eqcs = mechanism.eqconstraints
    eqcids = getfield.(eqcs, :id)

    nu = 0
    for id in eqcids
        eqc = geteqconstraint(mechanism, id)
        nu += getcontroldim(eqc)
    end

    Bcontrol = zeros(T, Nb*6, nu)

    n1 = 1
    n2 = 0
    for id in eqcids
        eqc = geteqconstraint(mechanism, id)
        n2 += getcontroldim(eqc)

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
    return Bcontrol
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
        d += getcontroldim(eqc)
    end
    return d
end

function soldim(mechanism::Mechanism)
    d = 0
    d += 6 * length(mechanism.bodies)
    d += sum(getdim.(mechanism.eqconstraints))
    !isempty(mechanism.ineqconstraints) && (d += sum(getdim.(mechanism.ineqconstraints)))
    return d
end

function getdim(eqc::EqualityConstraint{T,N,Nc,Cs}) where {T,N,Nc,Cs}
    return N
end

function getdim(ineqc::InequalityConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs,N½}
    return N
end

function getcontroldim(eqc::EqualityConstraint{T,N,Nc,Cs}; floatingbase::Bool = true) where {T,N,Nc,Cs}
    !floatingbase && (N == 0) && return 0
    cnt = 0
    for joint in eqc.constraints
        cnt += length(joint)
    end
    return 6 - cnt
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
        body.state.xc = x2
        body.state.vc = v1
        body.state.qc = UnitQuaternion(q2...)
        body.state.ωc = ω1
        body.state.xk[1] = x2
        body.state.qk[1] = UnitQuaternion(q2...)
        setsolution!(body)
        body.state.Fk[1] = SVector{3}([0,0,0.])
        body.state.τk[1] = SVector{3}([0,0,0.])
    end
    for eqc in mechanism.eqconstraints
        dim = getcontroldim(eqc)
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
        x2 = body.state.xc
        v1 = body.state.vc
        qc = body.state.qc
        q2 = [qc.w, qc.x, qc.y, qc.z]
        ω1 = body.state.ωc
        push!(data, [x2; v1; q2; ω1]...)
    end
    for eqc in mechanism.eqconstraints
        if getcontroldim(eqc) > 0
            # tra = eqc.constraints[1]
            # rot = eqc.constraints[2]
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
        nλ = getdim(eqc)
        λ = sol[off .+ (1:nλ)]; off += nλ
        eqc.λsol[2] = λ
    end
    for (i,body) in enumerate(mechanism.bodies)
        nv = 3
        nω = 3
        v2 = sol[off .+ (1:nv)]; off += nv
        ω2 = sol[off .+ (1:nω)]; off += nω
        body.state.vsol[2] = v2
        body.state.ωsol[2] = ω2
    end
    for (i,ineqc) in enumerate(mechanism.ineqconstraints)
        N = getdim(ineqc)
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
        ω2 = body.state.ωsol[2]
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
# Integration Scheme
################################################################################
using Test

getx3(x2, v2, Δt) = x2 + v2 * Δt
getq3(q2, ω2, Δt) = UnitQuaternion(q2, false) * ωbar(ω2, Δt) * Δt / 2

function st(x2, v2, Δt)
    x3 = x2 + Δt * v2
    return x3
end

function sr(q2, ω2, Δt; normalize::Bool = false)
    L = Lmat(UnitQuaternion(q2..., normalize))
    q3 = Δt / 2 * L * [sqrt(4 / Δt^2  - ω2' * ω2); ω2]
    return q3
end

function ∂st∂x(x2, v2, Δt)
    return I(3)
end

function ∂st∂v(x2, v2, Δt)
    return Δt * I(3)
end

function ∂sr∂q(q2, ω2, Δt)
    x = Δt / 2 * [sqrt(4 / Δt^2  - ω2' * ω2); ω2]
    return SA[
            x[1]  -x[2]  -x[3]  -x[4];
            x[2]   x[1]   x[4]  -x[3];
            x[3]  -x[4]   x[1]   x[2];
            x[4]   x[3]  -x[2]   x[1];
    ]
end

function ∂sr∂ω(q2, ω2, Δt)
    L = Lmat(UnitQuaternion(q2...))
    return Δt / 2 * L * [- ω2' / sqrt(4 / Δt^2  - ω2' * ω2); I(3)]
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
