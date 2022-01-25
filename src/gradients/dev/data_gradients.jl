################################################################################
# Data Jacobians
################################################################################
function ∂eqc∂body_data(mechanism::Mechanism, eqc::JointConstraint{T,N},
        body::Body{T}) where {T,N}
    Nd = data_dim(body)
    ∇m = szeros(T,N,1)
    ∇J = szeros(T,N,6)
    ∇z1 = szeros(T,N,6)
    ∇z2 = ∂g∂z(mechanism, eqc, body) * ∂i∂z(body, mechanism.Δt, attjac=true)
    ∇g = [∇m ∇J ∇z1 ∇z2]
    return ∇g
end

function ∂eqc∂eqc_data(mechanism::Mechanism, eqc::JointConstraint{T,N}) where {T,N}
    Nd = data_dim(eqc)
    return szeros(T,N,Nd)
end


function ∂body∂body_data(mechanism::Mechanism, body::Body{T}) where T
    Nd = data_dim(body)
    N = 6
    x1, v15, q1, ϕ15 = fullargs1(body.state)
    x2, v25, q2, ϕ25 = fullargssol(body.state)
    x3, q3 = posargs3(body.state, Δt)
    # Mass
    ezg = SA{T}[0; 0; -mechanism.g]
    ∇m = [- 1 / Δt * (x2 - x1) + Δt/2 * ezg + 1 / Δt * (x3 - x2) + Δt/2 * ezg;
          szeros(T,3,1)]
    # Inertia
    ∇J = -4 / Δt * LVᵀmat(q2)' * LVᵀmat(q1) * ∂inertia(VLᵀmat(q1) * vector(q2))
    ∇J += -4 / Δt * LVᵀmat(q2)' * Tmat() * RᵀVᵀmat(q3) * ∂inertia(VLᵀmat(q2) * vector(q3))
    ∇J = [szeros(T,3,6); ∇J]

    # configuration 1: z1 = x1, q1
    ∇x1 = 1 / Δt * body.m * SMatrix{3,3,T,9}(Diagonal(sones(T,3)))
    ∇q1 = -4 / Δt * LVᵀmat(q2)' * ∂qLVᵀmat(body.J * VLᵀmat(q1) * vector(q2))
    ∇q1 += -4 / Δt * LVᵀmat(q2)' * LVᵀmat(q1) * body.J * ∂qVLᵀmat(vector(q2))
    ∇q1 *= LVᵀmat(q1) # attjac
    ∇z1 = [∇x1 szeros(T,3,3);
           szeros(T,3,3) ∇q1]

    # configuration 2: z2 = x2, q2
    # TODO
    ∇z2 = szeros(T,N,6)

    return [∇m ∇J ∇z1 ∇z2]
end

function ∂body∂eqc_data(mechanism::Mechanism{T}, eqc::JointConstraint{T},
        body::Body{T}) where {T}
    Nd = data_dim(eqc)
    N = 6
    x1, v15, q1, ϕ15 = fullargs1(body.state)
    x2, v25, q2, ϕ25 = fullargssol(body.state)
    x3, q3 = posargs3(body.state, Δt)
    ∇u = Diagonal(SVector{6,T}(-1,-1,-1,-2,-2,-2)) * ∂Fτ∂u(mechanism, eqc, body)
    ∇spring = springforce(mechanism, eqc, body, unitary=true)
    ∇damper = damperforce(mechanism, eqc, body, unitary=true)
    # TODO
    nu = controldim(eqc)
    ∇spring_offset = szeros(T,N,nu)
    return [∇u ∇spring ∇damper ∇spring_offset]
end

function ∂body∂ineqc_data(mechanism::Mechanism, ineqc::ContactConstraint{T,N,Nc,Cs,N½},
        body::Body{T}) where {T,N,Nc,Cs<:Tuple{NonlinearContact{T,N}},N½}
    Nd = data_dim(ineqc)
    bound = ineqc.constraints[1]
    offset = bound.offset
    x3, q3 = posargs3(body.state, mechanism.Δt)
    γ = ineqc.γsol[2]

    ∇cf = szeros(T,3,1)

    X = forcemapping(bound)
    # this what we differentiate: Qᵀγ = - skew(p - vrotate(offset, inv(q3))) * VRmat(q3) * LᵀVᵀmat(q3) * X' * γ
    ∇p = - ∂pskew(VRmat(q3) * LᵀVᵀmat(q3) * X' * γ)
    ∇off = - ∂pskew(VRmat(q3) * LᵀVᵀmat(q3) * X' * γ) * -∂vrotate∂p(offset, inv(q3))

    ∇X = szeros(T,3,Nd)
    ∇Q = [∇cf ∇p ∇off]
    return [∇X; ∇Q]
end


function ∂ineqc∂ineqc_data(mechanism::Mechanism, ineqc::ContactConstraint{T,N,Nc,Cs,N½},
        body::Body{T}) where {T,N,Nc,Cs<:Tuple{NonlinearContact{T,N}},N½}
    Nd = data_dim(ineqc)
    bound = ineqc.constraints[1]
    p = bound.p
    offset = bound.offset
    x2, v25, q2, ϕ25 = fullargssol(body.state)
    x3, q3 = posargs3(body.state, mechanism.Δt)
    s = ineqc.ssol[2]
    γ = ineqc.γsol[2]

    ∇cf = SA[0,γ[1],0,0]
    ∇off = [-bound.ainv3; szeros(T,1,3); -bound.Bx * skew(vrotate(ϕ25, q3))]
    ∇p = [bound.ainv3 * ∂vrotate∂p(bound.p, q3); szeros(T,1,3); bound.Bx * skew(vrotate(ϕ25, q3)) * ∂vrotate∂p(bound.p, q3)]

    ∇compμ = szeros(T,N½,Nd)
    ∇g = [∇cf ∇p ∇off]
    return [∇compμ; ∇g]
end

function ∂ineqc∂body_data(mechanism::Mechanism, ineqc::ContactConstraint{T,N,Nc,Cs,N½},
        body::Body{T}) where {T,N,Nc,Cs,N½}
    Nd = data_dim(body)
    ∇compμ = szeros(T,N½,Nd)
    ∇m = szeros(T,N½,1)
    ∇J = szeros(T,N½,6)
    ∇z1 = szeros(T,N½,6)
    ∇z3 = ∂g∂z(mechanism, ineqc, body)
    ∇z2 = ∇z3 * ∂i∂z(body, mechanism.Δt, attjac=true) # 4x7 * 7x6 = 4x6
    ∇g = [∇m ∇J ∇z1 ∇z2]
    return [∇compμ; ∇g]
end

################################################################################
# System Data Jacobians
################################################################################
function ∂ineqc_data!(data_system::System, mechanism::Mechanism{T}) where {T}
    # ∂body∂ineqcdata
    for ineqc in mechanism.ineqconstraints
        pbody = getbody(mechanism, ineqc.parentid)
        data_system.matrix_entries[pbody.id, ineqc.id].value += ∂body∂ineqc_data(mechanism, ineqc, pbody)
    end
    # ∂ineqc∂ineqcdata
    for ineqc in mechanism.ineqconstraints
        pbody = getbody(mechanism, ineqc.parentid)
        data_system.matrix_entries[ineqc.id, ineqc.id].value += ∂ineqc∂ineqc_data(mechanism, ineqc, pbody)
    end
    return nothing
end

function ∂eqc_data!(data_system::System, mechanism::Mechanism{T}) where {T}
    # ∂eqc∂eqcdata
    for eqc in mechanism.eqconstraints
        data_system.matrix_entries[eqc.id, eqc.id].value += ∂eqc∂eqc_data(mechanism, eqc)
    end
    # ∂body∂eqcdata
    # TODO adapt this to handle cycles
    for body in mechanism.bodies
        for eqc in mechanism.eqconstraints
            if (body.id == eqc.parentid) || (body.id ∈ eqc.childids)
                data_system.matrix_entries[body.id, eqc.id].value += ∂body∂eqc_data(mechanism, eqc, body)
            end
        end
    end
    return nothing
end

function ∂body_data!(data_system::System, mechanism::Mechanism{T}) where {T}
    # ∂eqc∂bodydata
    # TODO adapt this to handle cycles
    for body in mechanism.bodies
        for eqc in mechanism.eqconstraints
            if (body.id == eqc.parentid) || (body.id ∈ eqc.childids)
                data_system.matrix_entries[eqc.id, body.id].value += ∂eqc∂body_data(mechanism, eqc, body)
            end
        end
    end
    # ∂body∂bodydata
    for body in mechanism.bodies
        data_system.matrix_entries[body.id, body.id].value += ∂body∂body_data(mechanism, body)
    end
    # ∂ineqc∂bodydata
    for ineqc in mechanism.ineqconstraints
        pbody = getbody(mechanism, ineqc.parentid)
        data_system.matrix_entries[ineqc.id, pbody.id].value += ∂ineqc∂body_data(mechanism, ineqc, pbody)
    end
    return nothing
end
