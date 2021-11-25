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

function minimalCoordinatesDict(mechanism::Mechanism)
    keys = mechanism.eqconstraints.keys
    values = Vector{SVector}()

    for eqc in mechanism.eqconstraints
        push!(values, minimalCoordinates(mechanism, eqc))
    end

    return UnitDict(keys, values)

end

function minimalCoordinates(mechanism::Mechanism{T}) where {T}
    N = controldim(mechanism)
    x = zeros(T,N)
    off = 0
    for eqc in mechanism.eqconstraints
        n = controldim(eqc)
        x[off .+ (1:n)] += minimalCoordinates(mechanism, eqc)
        off += n
    end
    return x
end

function minimalVelocities(mechanism::Mechanism{T}) where {T}
    N = controldim(mech)
    x = zeros(T,N)
    off = 0
    for eqc in mechanism.eqconstraints
        n = controldim(eqc)
        x[off .+ (1:n)] += minimalVelocities(mechanism, eqc)
        off += n
    end
    return x
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
