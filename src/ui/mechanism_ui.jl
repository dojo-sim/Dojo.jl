function verifyConstraints!(mechanism::Mechanism)
    for eqc in mechanism.eqconstraints
        typeof(eqc) <: EqualityConstraint{T,N,Nc,Cs} where {T,N,Nc,Cs<:Tuple{<:Friction}} && continue
        if norm(gc(mechanism, eqc)) > 1e-3
            @info string("Bad constraint satisfaction at constraint: ", eqc.id, ", |g| = ", norm(gc(mechanism, eqc)))
        end
    end
    return
end


function nameiddict(mechanism::Mechanism)
    dict = Dict{String,Int64}()
    for (id,body) in pairs(mechanism.bodies)
        if body.name != ""
            dict[body.name] = id
        end
    end
    for (id,eqc) in pairs(mechanism.eqconstraints)
        if eqc.name != ""
            dict[eqc.name] = id
        end
    end
    for (id,ineqc) in pairs(mechanism.ineqconstraints)
        if ineqc.name != ""
            dict[ineqc.name] = id
        end
    end

    return dict
end

function minimalCoordinates(mechanism::Mechanism)
    keys = mechanism.eqconstraints.keys
    values = Vector{SVector}()

    for eqc in mechanism.eqconstraints
        push!(values, minimalCoordinates(mechanism, eqc))
    end

    return UnitDict(keys, values)

end

function setPosition!(mechanism::Mechanism, dict::UnitDict)
    for (id,eqc) in pairs(mechanism.eqconstraints)
        setPosition!(mechanism, eqc, dict[id])
    end

    return
end

function setVelocity!(mechanism::Mechanism, dict::UnitDict)
    for (id,eqc) in pairs(mechanism.eqconstraints)
        setVelocity!(mechanism, eqc, dict[id])
    end

    return
end

function setForce!(mechanism::Mechanism, dict::UnitDict)
    for (id,eqc) in pairs(mechanism.eqconstraints)
        setForce!(mechanism, eqc, dict[id])
    end

    return
end

@inline function currentasknot!(mechanism::Mechanism)
    foreach(currentasknot!, mechanism.bodies)
    return
end