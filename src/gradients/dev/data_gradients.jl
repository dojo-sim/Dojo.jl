################################################################################
# Data Jacobians
################################################################################
function ∂eqc∂body_data(mechanism::Mechanism, eqc::EqualityConstraint{T,N},
        body::Body{T}) where {T,N}
    Nd = data_dim(body)
    ∇m = szeros(T,N,1)
    ∇J = szeros(T,N,6)
    ∇v15 = szeros(T,N,3)
    ∇ϕ15 = szeros(T,N,3)
    ∇z2 = -∂g∂z(mechanism, eqc, body) * ∂i∂z(body, mechanism.Δt, attjac=true)
    ∇g = [∇m ∇J ∇v15 ∇ϕ15 ∇z2]
    return ∇g
end

function ∂eqc∂eqc_data(mechanism::Mechanism, eqc::EqualityConstraint{T,N}) where {T,N}
    Nd = data_dim(eqc)
    return szeros(T,N,Nd)
end


function ∂body∂body_data(mechanism::Mechanism, body::Body{T}) where T
    Δt = mechanism.Δt
    Nd = data_dim(body)
    N = 6
    x1, v15, q1, ϕ15 = fullargs1(body.state)
    x2, v25, q2, ϕ25 = fullargssol(body.state)
    x3, q3 = posargs3(body.state, Δt)
    # Mass
    ezg = SA{T}[0; 0; -mechanism.g]
    ∇m = [1 / Δt * (x2 - x1) - Δt/2 * ezg - 1 / Δt * (x3 - x2) - Δt/2 * ezg;
          szeros(T,3,1)]
    # Inertia
    ∇J = 4 / Δt * LVᵀmat(q2)' * LVᵀmat(q1) * ∂inertia(VLᵀmat(q1) * vector(q2))
    ∇J += 4 / Δt * LVᵀmat(q2)' * Tmat() * RᵀVᵀmat(q3) * ∂inertia(VLᵀmat(q2) * vector(q3))
    ∇J = [szeros(T,3,6); ∇J]

    # initial conditions: v15, ϕ15
    ∇v15 = body.m * SMatrix{3,3,T,9}(Diagonal(sones(T,3)))
    ∇q1 = -4 / Δt * LVᵀmat(q2)' * ∂qLVᵀmat(body.J * VLᵀmat(q1) * vector(q2))
    ∇q1 += -4 / Δt * LVᵀmat(q2)' * LVᵀmat(q1) * body.J * ∂qVLᵀmat(vector(q2))
    ∇ϕ15 = ∇q1 * ∂integrator∂ϕ(q2, -ϕ15, Δt)
    ∇15 = [∇v15 szeros(T,3,3);
           szeros(T,3,3) ∇ϕ15]

    # current configuration: z2 = x2, q2
    ∇tra_x2 = - 2 / Δt * body.m * SMatrix{3,3,T,9}(Diagonal(sones(T,3)))
    ∇tra_q2 = szeros(T,3,3)
    ∇rot_x2 = szeros(T,3,3)
    ∇rot_q2 = -4 / Δt * VLᵀmat(q2) * LVᵀmat(q1) * body.J * VLᵀmat(q1)
    ∇rot_q2 += -4 / Δt * VLᵀmat(q2) * Tmat() * RᵀVᵀmat(q3) * body.J * ∂qVLᵀmat(vector(q3))
    ∇rot_q2 += -4 / Δt * ∂qVLᵀmat(LVᵀmat(q1) * body.J * VLᵀmat(q1) * vector(q2) + Tmat() * RᵀVᵀmat(q3) * body.J * VLᵀmat(q2) * vector(q3))
    ∇rot_q2 *= LVᵀmat(q2)
    ∇z2 = [∇tra_x2 ∇tra_q2;
           ∇rot_x2 ∇rot_q2]
    return [∇m ∇J ∇15 ∇z2]
end

LVᵀmat(q2)' - VLᵀmat(q2)

function ∂body∂eqc_data(mechanism::Mechanism{T}, eqc::EqualityConstraint{T},
        body::Body{T}) where {T}
    Δt = mechanism.Δt
    Nd = data_dim(eqc)
    N = 6
    x1, v15, q1, ϕ15 = fullargs1(body.state)
    x2, v25, q2, ϕ25 = fullargssol(body.state)
    x3, q3 = posargs3(body.state, Δt)
    ∇u = Diagonal(SVector{6,T}(1,1,1,2,2,2)) * ∂Fτ∂u(mechanism, eqc, body)
    ∇spring = springforce(mechanism, eqc, body, unitary=true)
    ∇damper = damperforce(mechanism, eqc, body, unitary=true)
    return [∇u ∇spring ∇damper]
end

function ∂body∂ineqc_data(mechanism::Mechanism, ineqc::InequalityConstraint{T,N,Nc,Cs,N½},
        body::Body{T}) where {T,N,Nc,Cs<:Tuple{ContactBound{T,N}},N½}
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


function ∂ineqc∂ineqc_data(mechanism::Mechanism, ineqc::InequalityConstraint{T,N,Nc,Cs,N½},
        body::Body{T}) where {T,N,Nc,Cs<:Tuple{ContactBound{T,N}},N½}
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

function ∂ineqc∂body_data(mechanism::Mechanism, ineqc::InequalityConstraint{T,N,Nc,Cs,N½},
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
function create_data_system(eqcs::Vector{<:EqualityConstraint}, bodies::Vector{<:Body},
        ineqcs::Vector{<:InequalityConstraint})
    nodes = [eqcs; bodies; ineqcs]
    A = adjacencyMatrix(eqcs, bodies, ineqcs)
    dimrow = length.(nodes)
    dimcol = data_dim.(nodes)
    data_system = System(A, dimrow, dimcol)
    return data_system
end

function ∂data!(data_system::System, mech::Mechanism)
    ∂ineqc_data!(data_system, mech::Mechanism)
    ∂body_data!(data_system, mech::Mechanism)
    ∂eqc_data!(data_system, mech::Mechanism)
    return nothing
end

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


function ∂body∂z(body::Body{T}, Δt::T; attjac::Bool=true) where T
    state = body.state
    q2 = state.q2[1]
    # ϕ25 = state.ϕsol[2]
    Z3 = szeros(T,3,3)
    Z34 = szeros(T,3,4)
    ZT = attjac ? szeros(T,6,6) : szeros(T,6,7)
    ZR = szeros(T,7,6)

    x1, q1 = posargs1(state)
    x2, q2 = posargs2(state)
    x3, q3 = posargs3(state, Δt)

    AposT = [-I Z3]
    # AvelT = [Z3 -I*body.m] # solving for impulses

    AposR = [-∂integrator∂q(q2, ϕ25, Δt, attjac = attjac) szeros(4,3)]

    rot_q1(q) = -4 / Δt * LVᵀmat(q2)' * Lmat(UnitQuaternion(q..., false)) * Vᵀmat() * body.J * Vmat() * Lmat(UnitQuaternion(q..., false))' * vector(q2)
    rot_q2(q) = -4 / Δt * LVᵀmat(UnitQuaternion(q..., false))' * Tmat() * Rmat(getq3(UnitQuaternion(q..., false), state.ϕsol[2], Δt))' * Vᵀmat() * body.J * Vmat() * Lmat(UnitQuaternion(q..., false))' * vector(getq3(UnitQuaternion(q..., false), state.ϕsol[2], Δt)) + -4 / Δt * LVᵀmat(UnitQuaternion(q..., false))' * Lmat(getq3(UnitQuaternion(q..., false), -state.ϕ15, Δt)) * Vᵀmat() * body.J * Vmat() * Lmat(getq3(UnitQuaternion(q..., false), -state.ϕ15, Δt))' * q

    # dynR_ϕ15 = -1.0 * FiniteDiff.finite_difference_jacobian(rot_q1, vector(q1)) * ∂integrator∂ϕ(q2, -state.ϕ15, Δt)
    dynR_q2 = FiniteDiff.finite_difference_jacobian(rot_q2, vector(q2))
    AvelR = attjac ? [dynR_q2 * LVᵀmat(q2) dynR_ϕ15] : [dynR_q2 dynR_ϕ15]

    return [[AposT;AvelT] ZT;
             ZR [AposR;AvelR]]
end
