function gc(mechanism::Mechanism{T}) where T
    rangeDict = Dict{Int64,UnitRange}()
    ind1 = 1
    ind2 = 0

    for (i,eqc) in enumerate(mechanism.eqconstraints)
        
        ind2 += length(eqc)
        range = ind1:ind2
        rangeDict[i] = range
        ind1 = ind2+1
    end

    gval = zeros(T,ind2)
    
    for (i,eqc) in enumerate(mechanism.eqconstraints)

        gval[rangeDict[i]] = gc(mechanism, eqc)
        
    end 
    return gval
end   

function constraintstep!(mechanism::Mechanism{T,Nn,Ne,Nb},freeids) where {T,Nn,Ne,Nb}
    freebodies = mechanism.bodies[freeids]

    gval=gc(mechanism)
    pinv∂g∂posc=pinv(∂g∂posc(mechanism,freeids))
    stepvec = -pinv∂g∂posc*gval

    for (i,body) in enumerate(freebodies)
            body.state.vsol[1] = stepvec[offsetrange(i, 3, 7, 1)]

            # Limit quaternion step to feasible length
            range = offsetrange(i, 3, 7, 2)
            Δstemp = VLᵀmat(body.state.qk[1]) * stepvec[first(range):(last(range)+1)]
            if norm(Δstemp) > 1
                Δstemp = Δstemp/norm(Δstemp)
            end
            body.state.ωsol[1] = Δstemp
        
    end

    return
end

function initializeConstraints!(mechanism::Mechanism{T,Nn,Ne,Nb}; fixedids = Int64[], freeids = Int64[], ε = 1e-5, newtonIter = 100, lineIter = 10) where {T,Nn,Ne,Nb}
    freebodies = Body[]
    if !isempty(fixedids) && !isempty(freeids)
        error("Specify either free or fixed bodies, not both.")
    elseif !isempty(fixedids)
        freeids = setdiff(getid.(mechanism.bodies),fixedids)
        freebodies = mechanism.bodies[freeids]
    elseif !isempty(freeids)
        freebodies = mechanism.bodies[freeids]
    else
        freeids = getid.(mechanism.bodies)
        freebodies = mechanism.bodies[freeids]
    end
    

    norm0 = norm(gc(mechanism))
    norm1 = norm0
    for n = Base.OneTo(newtonIter)

        for body in freebodies
            body.state.xk[1] = body.state.xc
            body.state.qk[1] = body.state.qc
        end

        constraintstep!(mechanism,freeids) 
    
        # Line search
        for j = Base.OneTo(lineIter)

            for body in freebodies
                
                body.state.xc = body.state.xk[1] + body.state.vsol[1]/(2^(j-1))
                
                w = sqrt(1-norm(body.state.ωsol[1]/(2^(j-1)))^2)
                body.state.qc = body.state.qk[1] * UnitQuaternion(w,body.state.ωsol[1]/(2^(j-1))...,false)
            end
            
            norm1 = norm(gc(mechanism))

            if norm1 < norm0 
                break
            end
        end

        if norm1 < ε
            return
        else
            norm0 = norm1
        end
    end

    display("Constraint initialization did not converge! Tolerance: "*string(norm1))
    return 
end
