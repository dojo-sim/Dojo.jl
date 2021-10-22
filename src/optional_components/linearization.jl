# TODO only works for 1DOF active constraints (eqcids)
# Maximal Coordinates
function linearsystem(mechanism::Mechanism{T,Nn,Ne,Nb}, xd, vd, qd, ωd, Fτd, bodyids, eqcids) where {T,Nn,Ne,Nb}
    statesold = [State{T}() for i=1:Nb]

    # store old state and set new initial state
    for (i,id) in enumerate(bodyids)
        stateold = settempvars!(getbody(mechanism, id), xd[i], vd[i], zeros(T,3), qd[i], ωd[i], zeros(T,3), zeros(T,6))
        statesold[i] = stateold
    end
    for (i,id) in enumerate(eqcids)
        setForce!(mechanism, geteqconstraint(mechanism, id), Fτd[i])
    end

    A, Bu, Bλ, G = lineardynamics(mechanism, eqcids)

    # restore old state
    for (i,id) in enumerate(bodyids)
        body = getbody(mechanism, id)
        body.state = statesold[i]
    end

    return A, Bu, Bλ, G
end
# Minimal Coordinates
function linearsystem(mechanism::Mechanism{T,Nn,Ne,Nb}, xθd, vωd, Fτd, controlledids, controlids) where {T,Nn,Ne,Nb}
    statesold = [State{T}() for i=1:Nb]
    xd = [szeros(T,3) for i=1:Nb]
    vd = [szeros(T,3) for i=1:Nb]
    qd = [one(UnitQuaternion{T}) for i=1:Nb]
    ωd = [szeros(T,3) for i=1:Nb]

    # store old state and set new initial state
    for (id,body) in pairs(mechanism.bodies)
        statesold[id] = deepcopy(body.state)
    end
    for (i,id) in enumerate(controlledids)
        setPosition!(mechanism, geteqconstraint(mechanism, id), [xθd[i]])
    end
    for (id,body) in pairs(mechanism.bodies)
        state = body.state
        xd[id] = state.xc
        vd[id] = state.vc
        qd[id] = state.qc
        ωd[id] = state.ωc
    end
    for (i,id) in enumerate(controlids)
        setForce!(mechanism, geteqconstraint(mechanism, id), [Fτd[i]])
    end

    A, Bu, Bλ, G = lineardynamics(mechanism, controlids)

    # restore old state
    for (id,body) in pairs(mechanism.bodies)
        body.state = statesold[id]
    end

    return A, Bu, Bλ, G, xd, vd, qd, ωd
end


function lineardynamics(mechanism::Mechanism{T,Nn,Ne,Nb}, eqcids) where {T,Nn,Ne,Nb}
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


    # calculate next state
    discretizestate!(mechanism)
    foreach(setsolution!, mechanism.bodies)
    foreach(applyFτ!, eqcs, mechanism, false)
    newton!(mechanism)

    # get state linearization
    Fz = zeros(T,Nb*13,Nb*12)
    Fu = zeros(T,Nb*13,Nb*6)
    Fλ = zeros(T,Nb*13,nc)
    Ffz = zeros(T,Nb*13,Nb*13)
    invFfzquat = zeros(T,Nb*12,Nb*13)

    Bcontrol = zeros(T,Nb*6,nu)

    for (ind,body) in enumerate(bodies)
        col6 = offsetrange(ind,6)
        col12 = offsetrange(ind,12)
        col13 = offsetrange(ind,13)

        Fzi = ∂F∂z(body, Δt)
        Fui = ∂F∂u(body, Δt)
        Ffzi, invFfzquati = ∂F∂fz(body, Δt)

        Fz[col13,col12] = Fzi
        Fu[col13,col6] = Fui
        Ffz[col13,col13] = Ffzi
        invFfzquat[col12,col13] = invFfzquati
    end

    Ffz += linearconstraintmapping(mechanism)
    Fz += linearforcemapping(mechanism)

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

    G, Fλ = linearconstraints(mechanism)

    invFfz = invFfzquat * inv(Ffz)
    A = -invFfz * Fz
    Bu = -invFfz * Fu * Bcontrol
    Bλ = -invFfz * Fλ

    return A, Bu, Bλ, G
end

function linearconstraints(mechanism::Mechanism{T,Nn,Ne,Nb}) where {T,Nn,Ne,Nb}
    Δt = mechanism.Δt
    eqcs = mechanism.eqconstraints

    nc = 0
    for eqc in eqcs
        nc += length(eqc)
    end

    Gl = zeros(T,nc,Nb*12)
    Gr = zeros(T,nc,Nb*13)

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

                pcol3b = offsetrange(parentind,3,13,2)
                pcol3d = offsetrange(parentind,3,13,4)
                pcol3d = first(pcol3d)+1:last(pcol3d)+1

                ccol3a12 = offsetrange(childind,3,12,1)
                ccol3b12 = offsetrange(childind,3,12,2)
                ccol3c12 = offsetrange(childind,3,12,3)
                ccol3d12 = offsetrange(childind,3,12,4)

                ccol3b = offsetrange(childind,3,13,2)
                ccol3d = offsetrange(childind,3,13,4)
                ccol3d = first(ccol3d)+1:last(ccol3d)+1

                pXl, pQl = ∂g∂posa(eqc.constraints[i], posargsnext(pstate, Δt)..., posargsnext(cstate, Δt)...) # x3
                pGr = ∂g∂ʳposa(eqc.constraints[i], posargssol(pstate)..., posargssol(cstate)...) # x2
                pXr, pQr = pGr[:,1:3], pGr[:,4:6]
                cXl, cQl =  ∂g∂posb(eqc.constraints[i], posargsnext(pstate, Δt)..., posargsnext(cstate, Δt)...) # x3
                cGr =  ∂g∂ʳposb(eqc.constraints[i], posargssol(pstate)..., posargssol(cstate)...) # x2
                cXr, cQr = cGr[:,1:3], cGr[:,4:6]

                mat = constraintmat(eqc.constraints[i])
                pGlx = mat * pXl
                pGlq = mat * pQl
                pGrx = mat * pXr
                pGrq = mat * pQr
                cGlx = mat * cXl
                cGlq = mat * cQl
                cGrx = mat * cXr
                cGrq = mat * cQr

                Gl[range,pcol3a12] = pGlx
                Gl[range,pcol3b12] = pGlx*Δt
                Gl[range,pcol3c12] = pGlq*Rmat(ωbar(pstate.ωsol[2],Δt)*Δt/2)*LVᵀmat(pstate.qsol[2])
                Gl[range,pcol3d12] = pGlq*Lmat(pstate.qsol[2])*derivωbar(pstate.ωsol[2],Δt)*Δt/2
                Gr[range,pcol3b] = pGrx
                Gr[range,pcol3d] = pGrq

                Gl[range,ccol3a12] = cGlx
                Gl[range,ccol3b12] = cGlx*Δt
                Gl[range,ccol3c12] = cGlq*Rmat(ωbar(cstate.ωsol[2],Δt)*Δt/2)*LVᵀmat(cstate.qsol[2])
                Gl[range,ccol3d12] = cGlq*Lmat(cstate.qsol[2])*derivωbar(cstate.ωsol[2],Δt)*Δt/2
                Gr[range,ccol3b] = cGrx
                Gr[range,ccol3d] = cGrq

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

                ccol3b = offsetrange(childind,3,13,2)
                ccol3d = offsetrange(childind,3,13,4)
                ccol3d = first(ccol3d)+1:last(ccol3d)+1


                cXl, cQl =  ∂g∂posb(eqc.constraints[i], posargsnext(cstate, Δt)...) # x3
                cGr =  ∂g∂ʳposb(eqc.constraints[i], posargssol(cstate)...) # x2
                cXr, cQr = cGr[:,1:3], cGr[:,4:6]

                mat = constraintmat(eqc.constraints[i])
                cGlx = mat * cXl
                cGlq = mat * cQl
                cGrx = mat * cXr
                cGrq = mat * cQr

                Gl[range,ccol3a12] = cGlx
                Gl[range,ccol3b12] = cGlx*Δt
                Gl[range,ccol3c12] = cGlq*Rmat(ωbar(cstate.ωsol[2],Δt)*Δt/2)*LVᵀmat(cstate.qsol[2])
                Gl[range,ccol3d12] = cGlq*Lmat(cstate.qsol[2])*derivωbar(cstate.ωsol[2],Δt)*Δt/2
                Gr[range,ccol3b] = cGrx
                Gr[range,ccol3d] = cGrq

                ind1 = ind2+1
            end
        end

        oneindc += length(eqc)
    end

    return Gl, -Gr'
end

function linearconstraintmapping(mechanism::Mechanism{T,Nn,Ne,Nb}) where {T,Nn,Ne,Nb}
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

                XX, XQ, QX, QQ = ∂2g∂posaa(constraint, state1.xsol[2], state1.qsol[2], state2.xsol[2], state2.qsol[2])
                Aaa[4:6,1:3] = kronproduct*XX
                Aaa[4:6,7:10] = kronproduct*XQ
                Aaa[11:13,1:3] = kronproduct*QX
                Aaa[11:13,7:10] = kronproduct*QQ

                XX, XQ, QX, QQ = ∂2g∂posab(constraint, state1.xsol[2], state1.qsol[2], state2.xsol[2], state2.qsol[2])
                Aab[4:6,1:3] = kronproduct*XX
                Aab[4:6,7:10] = kronproduct*XQ
                Aab[11:13,1:3] = kronproduct*QX
                Aab[11:13,7:10] = kronproduct*QQ

                XX, XQ, QX, QQ = ∂2g∂posba(constraint, state1.xsol[2], state1.qsol[2], state2.xsol[2], state2.qsol[2])
                Aba[4:6,1:3] = kronproduct*XX
                Aba[4:6,7:10] = kronproduct*XQ
                Aba[11:13,1:3] = kronproduct*QX
                Aba[11:13,7:10] = kronproduct*QQ

                XX, XQ, QX, QQ = ∂2g∂posbb(constraint, state1.xsol[2], state1.qsol[2], state2.xsol[2], state2.qsol[2])
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

                # @show typeof(λ')
                # @show typeof(constraintmat(constraint))
                # @show size(λ')
                # @show size(constraintmat(constraint))
                kronproduct = -kron(λ'*Array(constraintmat(constraint)),E)*K

                XX, XQ, QX, QQ = ∂2g∂posbb(constraint, state2.xsol[2], state2.qsol[2])
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

function linearforcemapping(mechanism::Mechanism{T,Nn,Ne,Nb}) where {T,Nn,Ne,Nb}
    Fzu = zeros(T,Nb*13,Nb*12)

    for eqc in mechanism.eqconstraints
        parentid = eqc.parentid
        for (i,childid) in enumerate(eqc.childids)
            childind = childid - Ne
            if parentid !== nothing
                parentind = parentid - Ne
                FaXa, FaQa, τaXa, τaQa, FbXa, FbQa, τbXa, τbQa = ∂Fτ∂posa(eqc.constraints[i], getbody(mechanism, parentid).state, getbody(mechanism, childid).state)
                FaXb, FaQb, τaXb, τaQb, FbXb, FbQb, τbXb, τbQb = ∂Fτ∂posb(eqc.constraints[i], getbody(mechanism, parentid).state, getbody(mechanism, childid).state)

                cola6 = offsetrange(parentind,6)
                colb6 = offsetrange(childind,6)
                rowav = offsetrange(parentind,3,13,2)
                rowaω = offsetrange(parentind,3,13,4).+1
                rowbv = offsetrange(childind,3,13,2)
                rowbω = offsetrange(childind,3,13,4).+1

                Fzu[rowav,cola6] = [FaXa FaQa]
                Fzu[rowaω,cola6] = [τaXa τaQa]
                Fzu[rowbv,cola6] = [FbXa FbQa]
                Fzu[rowbω,cola6] = [τbXa τbQa]

                Fzu[rowav,colb6] = [FaXb FaQb]
                Fzu[rowaω,colb6] = [τaXb τaQb]
                Fzu[rowbv,colb6] = [FbXb FbQb]
                Fzu[rowbω,colb6] = [τbXb τbQb]
            else
                FaXb, FaQb, τaXb, τaQb, FbXb, FbQb, τbXb, τbQb = ∂Fτ∂posb(eqc.constraints[i], getbody(mechanism, childid).state)

                colb6 = offsetrange(childind,6)
                rowbv = offsetrange(childind,3,13,2)
                rowbω = offsetrange(childind,3,13,4).+1

                Fzu[rowbv,colb6] = [FbXb FbQb]
                Fzu[rowbω,colb6] = [τbXb τbQb]
            end

        end
    end

    return Fzu
end
