using ConstrainedDynamics
using ConstrainedDynamics: posargsk, fullargssol, Lmat, VLᵀmat, Lᵀmat, ωbar
using Rotations
using Rotations: params
using StaticArrays



function getxqkvector(state)
    x, q = posargsk(state)
    q = params(q)
    
    return [x;q]
end

function getsolestimate(state)
    x, v, q, ω = fullargssol(state)
    q = params(q)

    return [x;v;q;ω]
end

getxk(state) = posargsk(state)[1]
getqk(state) = posargsk(state)[2]



function transfuncpos(vars)
    xa = vars[1:3]
    qa = UnitQuaternion(SA[vars[4:7]...])
    xb = vars[8:10]
    qb = UnitQuaternion(SA[vars[11:14]...])

    pa = vars[15:17]
    pb = vars[18:20]

    return vrotate(xb + vrotate(pb, qb) - (xa + vrotate(pa, qa)), inv(qa))
end

function transfunc3pos(vars)
    g = transfuncpos(vars)

    return g
end

function transfunc2pos(vars)
    g = transfuncpos(vars)
    V1 = vars[21:23]
    V2 = vars[24:26]
    V12 = [V1';V2']

    return V12 * g
end

function transfunc1pos(vars)
    g = transfuncpos(vars)
    V3 = vars[21:23]

    return V3' * g
end


function transfuncvel(vars)
    Δt = 0.01

    xa1 = vars[1:3]
    va1 = vars[4:6]
    qa1 = UnitQuaternion(SA[vars[7:10]...])
    ωa1 = vars[11:13]

    xb1 = vars[14:16]
    vb1 = vars[17:19]
    qb1 = UnitQuaternion(SA[vars[20:23]...])
    ωb1 = vars[24:26]

    xa2 = xa1 + va1 * Δt
    qa2 = qa1 * ωbar(ωa1, Δt) * Δt / 2

    xb2 = xb1 + vb1 * Δt
    qb2 = qb1 * ωbar(ωb1, Δt) * Δt / 2

    pa = vars[27:29]
    pb = vars[30:32]

    return vrotate(xb2 + vrotate(pb, qb2) - (xa2 + vrotate(pa, qa2)), inv(qa2))
end

function transfunc3vel(vars)
    g = transfuncvel(vars)

    return g
end

function transfunc2vel(vars)
    g = transfuncvel(vars)
    V1 = vars[33:35]
    V2 = vars[36:38]
    V12 = [V1';V2']

    return V12 * g
end

function transfunc1vel(vars)
    g = transfuncvel(vars)
    V3 = vars[33:35]

    return V3' * g
end



function rotfuncpos(vars)
    xa = vars[1:3]
    qa = UnitQuaternion(SA[vars[4:7]...], false)
    xb = vars[8:10]
    qb = UnitQuaternion(SA[vars[11:14]...], false)

    offset = UnitQuaternion(SA[vars[15:18]...], false)

    return VLᵀmat(qa) * Lmat(qb) * params(inv(offset))
end

function rotfunc3pos(vars)
    g = rotfuncpos(vars)

    return g
end

function rotfunc2pos(vars)
    g = rotfuncpos(vars)
    V1 = vars[19:21]
    V2 = vars[22:24]
    V12 = [V1';V2']

    return V12 * g
end

function rotfunc1pos(vars)
    g = rotfuncpos(vars)
    V3 = vars[19:21]

    return V3' * g
end


function rotfuncvel(vars)
    Δt = 0.01

    xa1 = vars[1:3]
    va1 = vars[4:6]
    qa1 = UnitQuaternion(SA[vars[7:10]...], false)
    ωa1 = vars[11:13]

    xb1 = vars[14:16]
    vb1 = vars[17:19]
    qb1 = UnitQuaternion(SA[vars[20:23]...], false)
    ωb1 = vars[24:26]

    xa2 = xa1 + va1 * Δt
    qa2 = qa1 * ωbar(ωa1, Δt) * Δt / 2

    xb2 = xb1 + vb1 * Δt
    qb2 = qb1 * ωbar(ωb1, Δt) * Δt / 2

    offset = UnitQuaternion(SA[vars[27:30]...], false)

    return VLᵀmat(qa2) * Lmat(qb2) * params(inv(offset))
end

function rotfunc3vel(vars)
    g = rotfuncvel(vars)

    return g
end

function rotfunc2vel(vars)
    g = rotfuncvel(vars)
    V1 = vars[31:33]
    V2 = vars[34:36]
    V12 = [V1';V2']    

    return V12 * g
end

function rotfunc1vel(vars)
    g = rotfuncvel(vars)
    V3 = vars[31:33]

    return V3' * g
end
