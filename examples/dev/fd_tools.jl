
function fdjac(f, x; δ = 1e-5)
    n = length(f(x))
    m = length(x)
    jac = zeros(n, m)
    for i = 1:m
        xp = deepcopy(x)
        xm = deepcopy(x)
        xp[i] += δ
        xm[i] -= δ
        jac[:,i] = (f(xp) - f(xm)) / (2δ)
    end
    return jac
end

function fd_diagonal∂damper∂ʳvel(j::Joint{T}, x2a::AbstractVector, q2a::UnitQuaternion, x2b::AbstractVector, q2b::UnitQuaternion,
    x1a::AbstractVector, v1a::AbstractVector, q1a::UnitQuaternion, ω1a::AbstractVector, x1b::AbstractVector, v1b::AbstractVector, q1b::UnitQuaternion, ω1b::AbstractVector, Δt) where T
    function f(vω1a)
        v1a = vω1a[1:3]
        ω1a = vω1a[4:6]
        Fτa = damperforcea(j, x2a, q2a, x2b, q2b, x1a, v1a, q1a, ω1a, x1b, v1b, q1b, ω1b, Δt)
        return Fτa
    end
    return fdjac(vω1a -> f(vω1a), Vector([v1a; ω1a]))
end

function fd_offdiagonal∂damper∂ʳvel(j::Joint{T}, x2a::AbstractVector, q2a::UnitQuaternion, x2b::AbstractVector, q2b::UnitQuaternion,
    x1a::AbstractVector, v1a::AbstractVector, q1a::UnitQuaternion, ω1a::AbstractVector, x1b::AbstractVector, v1b::AbstractVector, q1b::UnitQuaternion, ω1b::AbstractVector, Δt) where T
    function f(vω1b)
        v1b = vω1b[1:3]
        ω1b = vω1b[4:6]
        Fτb = damperforceb(j, x2a, q2a, x2b, q2b, x1a, v1a, q1a, ω1a, x1b, v1b, q1b, ω1b, Δt)
        return Fτb
    end
    return fdjac(vω1b -> f(vω1b), Vector([v1b; ω1b]))
end

function fd_offdiagonal∂damper∂ʳvel(j::Joint{T}, x2b::AbstractVector, q2b::UnitQuaternion,
    x1b::AbstractVector, v1b::AbstractVector, q1b::UnitQuaternion, ω1b::AbstractVector, Δt) where T
    function f(vω1b)
        v1b = vω1b[1:3]
        ω1b = vω1b[4:6]
        Fτb = damperforceb(j, x2b, q2b, x1b, v1b, q1b, ω1b, Δt)
        return Fτb
    end
    return fdjac(vω1b -> f(vω1b), Vector([v1b; ω1b]))
end


function fd_diagonal∂spring∂ʳvel(j::Joint{T}, x2a::AbstractVector, q2a::UnitQuaternion, x2b::AbstractVector, q2b::UnitQuaternion,
    x1a::AbstractVector, v1a::AbstractVector, q1a::UnitQuaternion, ω1a::AbstractVector, x1b::AbstractVector, v1b::AbstractVector, q1b::UnitQuaternion, ω1b::AbstractVector, Δt) where T
    function f(vω1a)
        v1a = vω1a[1:3]
        ω1a = vω1a[4:6]
        Fτa = springforcea(j, x2a, q2a, x2b, q2b, x1a, v1a, q1a, ω1a, x1b, v1b, q1b, ω1b, Δt)
        return Fτa
    end
    return fdjac(vω1a -> f(vω1a), Vector([v1a; ω1a]))
end

function fd_offdiagonal∂spring∂ʳvel(j::Joint{T}, x2a::AbstractVector, q2a::UnitQuaternion, x2b::AbstractVector, q2b::UnitQuaternion,
    x1a::AbstractVector, v1a::AbstractVector, q1a::UnitQuaternion, ω1a::AbstractVector, x1b::AbstractVector, v1b::AbstractVector, q1b::UnitQuaternion, ω1b::AbstractVector, Δt) where T
    function f(vω1b)
        v1b = vω1b[1:3]
        ω1b = vω1b[4:6]
        Fτb = springforceb(j, x2a, q2a, x2b, q2b, x1a, v1a, q1a, ω1a, x1b, v1b, q1b, ω1b, Δt)
        return Fτb
    end
    return fdjac(vω1b -> f(vω1b), Vector([v1b; ω1b]))
end

function fd_offdiagonal∂spring∂ʳvel(j::Joint{T}, x2b::AbstractVector, q2b::UnitQuaternion,
    x1b::AbstractVector, v1b::AbstractVector, q1b::UnitQuaternion, ω1b::AbstractVector, Δt) where T
    function f(vω1b)
        v1b = vω1b[1:3]
        ω1b = vω1b[4:6]
        Fτb = springforceb(j, x2b, q2b, x1b, v1b, q1b, ω1b, Δt)
        return Fτb
    end
    return fdjac(vω1b -> f(vω1b), Vector([v1b; ω1b]))
end
