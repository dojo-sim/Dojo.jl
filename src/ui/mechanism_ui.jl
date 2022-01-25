# function nameiddict(mechanism::Mechanism)
#     dict = Dict{String,Int64}()
#     for (id,body) in pairs(mechanism.bodies)
#         if body.name != ""
#             dict[body.name] = id
#         end
#     end
#     for (id,eqc) in pairs(mechanism.eqconstraints)
#         if eqc.name != ""
#             dict[eqc.name] = id
#         end
#     end
#     for (id,ineqc) in pairs(mechanism.ineqconstraints)
#         if ineqc.name != ""
#             dict[ineqc.name] = id
#         end
#     end

#     return dict
# end

function minimalCoordinates(mechanism::Mechanism)
    d = Dict()
    for eqc in mechanism.eqconstraints
        push!(d, eqc.id => minimalCoordinates(mechanism, eqc))
    end
    return d
end

function minimalCoordinatesVector(mechanism::Mechanism{T}) where {T}
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

function minimalVelocitiesVector(mechanism::Mechanism{T}) where {T}
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

function setPosition!(mechanism::Mechanism, dict)
    for (id,eqc) in pairs(mechanism.eqconstraints)
        setPosition!(mechanism, eqc, dict[id])
    end

    return
end

function setVelocity!(mechanism::Mechanism, dict)
    for (id,eqc) in pairs(mechanism.eqconstraints)
        setVelocity!(mechanism, eqc, dict[id])
    end

    return
end

function zeroVelocity!(mechanism::Mechanism)
    # velocities
    for (i, body) in enumerate(mechanism.bodies)
        try
            setVelocity!(body, v=zeros(3), ω=zeros(3))
        catch
            nothing
        end
    end
end

function setForce!(mechanism::Mechanism, dict)
    for (id,eqc) in pairs(mechanism.eqconstraints)
        setForce!(eqc, dict[id])
    end

    return
end

@inline function currentasknot!(mechanism::Mechanism)
    for body in mechanism.bodies currentasknot!(body) end
end
