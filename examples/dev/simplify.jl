using Dojo
using MeshCat

vis = Visualizer()
open(vis)


mech = getmechanism(:humanoid, contact=true, Δt=0.05, g=-9.81, spring=30.0, damper=5.0)
initialize!(mech, :humanoid, rot=[0.1,0,0], tran=[0,0,1.5])

function ctrl!(mechanism, k)
    nu = controldim(mechanism)
    u = szeros(nu)
    set_control!(mechanism, u)
    return
end

storage = simulate!(mech, 2.3, ctrl!, record=true, verbose=false)
visualize(mech, storage, vis=vis)


l = 1.0
m = 1.0
joint_axis = [1.0; 0; 0]
width, depth = 0.1, 0.1
p2 = [0; 0; l/2] # joint connection point

# Links
origin = Origin{Float64}()
body1 = Box(width, depth, l, m)

# Constraints
joint_between_origin_and_body1 = EqualityConstraint(Revolute(origin, body1,
    joint_axis; p2=p2, spring = 0, damper = 0,
    rot_joint_limits = [SVector{1}([0.25 * π]), SVector{1}([π])]
    ))
bodies = [body1]
eqcs = [joint_between_origin_and_body1]

length(joint_between_origin_and_body1)
length(joint_between_origin_and_body1.constraints[1])
length(joint_between_origin_and_body1.constraints[2])

################################################################################
# Analytical Jacobian
################################################################################

function create_data_system(eqcs::Vector{<:EqualityConstraint}, bodies::Vector{<:Body},
        ineqcs::Vector{<:InequalityConstraint})
    nodes = [eqcs; bodies; ineqcs]
    A = adjacencyMatrix(eqcs, bodies, ineqcs)
    dimrow = length.(nodes)
    dimcol = datadim.(nodes)
    data_system = System(A, dimrow, dimcol)
    return data_system
end

# Body
datadim(body::Body) = 21 # 1+6+7+7 [m,flat(J),x1,q1,x2,q2]
# Eqconstraints
datadim(eqc::EqualityConstraint) = 2 + sum(datadim.(eqc.constraints)) # [spring, damper, utra, urot, tra_spring_offset, rot_spring_offset]
datadim(joint::Rotational{T,Nλ,Nb,N,Nb½,N̄λ}) where {T,Nλ,Nb,N,Nb½,N̄λ} = 2N̄λ # [u, spring, damper, spring_offset]
datadim(joint::Translational{T,Nλ,Nb,N,Nb½,N̄λ}) where {T,Nλ,Nb,N,Nb½,N̄λ} = 2N̄λ # [u, spring, damper, spring_offset]
# Ineqconstraints
datadim(ineqc::InequalityConstraint) = sum(datadim.(ineqc.constraints))
datadim(bound::ContactBound) = 7 # [cf, p, offset]
datadim(bound::LinearContactBound) = 7 # [cf, p, offset]
datadim(bound::ImpactBound) = 6 # [p, offset]

function ineqc_dataineqc(mechanism::Mechanism, ineqc::InequalityConstraint{T,N,Nc,Cs,N½},
        body::Body{T}) where {T,N,Nc,Cs<:Tuple{ContactBound{T,N}},N½}
    Nd = datadim(ineqc)
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

function ineqc_databody(mechanism::Mechanism, ineqc::InequalityConstraint{T,N,Nc,Cs,N½},
        body::Body{T}) where {T,N,Nc,Cs,N½}
    Nd = datadim(body)
    Δt = mechanism.Δt
    bound = ineqc.constraints[1]
    state = body.state
    v25 = state.vsol[2]
    ϕ25 = state.ϕsol[2]
    x2, q2 = posargs2(state)

    ∇compμ = szeros(T,N½,Nd)
    ∇m = szeros(T,N½,1)
    ∇J = szeros(T,N½,6)
    ∇z1 = szeros(T,N½,7)
    ∇v = ∂g∂v(mechanism, ineqc, body)
    ∇z2 = ∇v * ∂i∂z(q2, ϕ25, Δt, attjac=true)' # 4x6 * 6x7 = 4x7
    ∇g = [∇m ∇J ∇z1 ∇z2]
    return [∇compμ; ∇g]
end

function eqc_databody(mechanism::Mechanism, eqc::EqualityConstraint{T,N},
        body::Body{T}) where {T,N}
    Nd = datadim(body)
    ∇m = szeros(T,N,1)
    ∇J = szeros(T,N,6)
    ∇z1 = szeros(T,N,7)
    ∇z2 = ∂g∂z(mechanism, eqc, body) * ∂i∂z(body, mechanism.Δt, attjac=false)
    ∇g = [∇m ∇J ∇z1 ∇z2]
    return ∇g
end

function eqc_dataeqc(mechanism::Mechanism, eqc::EqualityConstraint{T,N}) where {T,N}
    Nd = datadim(eqc)
    return szeros(T,N,Nd)
end

function body_databody(mechanism::Mechanism, body::Body{T}) where T
    Nd = datadim(body)
    N = 6
    x1, v15, q1, ϕ15 = fullargs1(body.state)
    x2, v25, q2, ϕ25 = fullargssol(body.state)
    x3, q3 = posargs3(body.state, Δt)
    # Mass
    ezg = SA{T}[0; 0; -mechanism.g]
    ∇m = [- 1 / Δt * (x2 - x1) + Δt/2 * ezg + 1 / Δt * (x3 - x2) + Δt/2 * ezg;
          szeros(T,3,1)]
    # Inertia
    # TODO
    ∇J = szeros(T,N,6)
    # TODO
    ∇z1 = szeros(T,N,7)
    # TODO
    ∇z2 = szeros(T,N,7)
    return [∇m ∇J ∇z1 ∇z2]
end

function body_dataeqc(mechanism::Mechanism{T}, eqc::EqualityConstraint{T},
        body::Body{T}) where {T}
    Nd = datadim(eqc)
    N = 6
    x1, v15, q1, ϕ15 = fullargs1(body.state)
    x2, v25, q2, ϕ25 = fullargssol(body.state)
    x3, q3 = posargs3(body.state, Δt)
    nu = controldim(eqc)
    @show nu
    # TODO
    ∇u = szeros(T,N,nu)
    # TODO
    ∇spring = szeros(T,N,1)
    # TODO
    ∇damper = szeros(T,N,1)
    # TODO
    ∇spring_offset = szeros(T,N,nu)
    return [∇u ∇spring ∇damper ∇spring_offset]
end

function body_dataineqc(mechanism::Mechanism, ineqc::InequalityConstraint{T,N,Nc,Cs,N½},
        body::Body{T}) where {T,N,Nc,Cs<:Tuple{ContactBound{T,N}},N½}
    Nd = datadim(ineqc)
    x1, v15, q1, ϕ15 = fullargs1(body.state)
    x2, v25, q2, ϕ25 = fullargssol(body.state)
    x3, q3 = posargs3(body.state, Δt)

    # TODO
    ∇cf = szeros(T,6,1)
    # TODO
    ∇p = szeros(T,6,3)
    # TODO
    ∇off = szeros(T,6,3)
    ∇g = [∇cf ∇p ∇off]
    return ∇g
end

function dataineqc!(data_system::System, mechanism::Mechanism{T}) where {T}
    # ∂body∂ineqcdata
    for ineqc in mechanism.ineqconstraints
        pbody = getbody(mechanism, ineqc.parentid)
        data_system.matrix_entries[pbody.id, ineqc.id].value += body_dataineqc(mechanism, ineqc, pbody)
    end
    # ∂ineqc∂ineqcdata
    for ineqc in mechanism.ineqconstraints
        pbody = getbody(mechanism, ineqc.parentid)
        data_system.matrix_entries[ineqc.id, ineqc.id].value += ineqc_dataineqc(mechanism, ineqc, pbody)
    end
    return nothing
end

function dataeqc!(data_system::System, mechanism::Mechanism{T}) where {T}
    # ∂eqc∂eqcdata
    for eqc in mechanism.eqconstraints
        data_system.matrix_entries[eqc.id, eqc.id].value += eqc_dataeqc(mechanism, eqc)
    end
    # ∂body∂eqcdata
    # TODO adapt this to handle cycles
    for body in mechanism.bodies
        for eqc in mechanism.eqconstraints
            if (body.id == eqc.parentid) || (body.id ∈ eqc.childids)
                data_system.matrix_entries[body.id, eqc.id].value += body_dataeqc(mechanism, eqc, body)
            end
        end
    end
    return nothing
end

function databody!(data_system::System, mechanism::Mechanism{T}) where {T}
    # ∂eqc∂bodydata
    # TODO adapt this to handle cycles
    for body in mechanism.bodies
        for eqc in mechanism.eqconstraints
            if (body.id == eqc.parentid) || (body.id ∈ eqc.childids)
                data_system.matrix_entries[eqc.id, body.id].value += eqc_databody(mechanism, eqc, body)
            end
        end
    end
    # ∂body∂bodydata
    for body in mechanism.bodies
        data_system.matrix_entries[body.id, body.id].value += body_databody(mechanism, body)
    end
    # ∂ineqc∂bodydata
    for ineqc in mechanism.ineqconstraints
        pbody = getbody(mechanism, ineqc.parentid)
        data_system.matrix_entries[ineqc.id, pbody.id].value += ineqc_databody(mechanism, ineqc, pbody)
    end
    return nothing
end

# vis = Visualizer()
# open(vis)

mech = getsnake(Nb=3, damper=0.0, spring=0.0);
function ctrl!(mech, k)
    nu = controldim(mech)
    u = mech.Δt * 0.00 * sones(nu)
    set_control!(mech, u)
    return nothing
end

initialize!(mech, :snake, x=[0,0,1.0])
storage = simulate!(mech, 1.35, ctrl!, record=true, verbose=false)
visualize(mech, storage, vis=vis)
ineqc0 = mech.ineqconstraints.values[1]
eqc0 = mech.eqconstraints.values[2]
body0 = mech.bodies.values[1]
∇v0 = ∂g∂v(mech, ineqc0, body0)
∇0 = ineqc_databody(mech, ineqc0, body0)
∇0 = ineqc_dataineqc(mech, ineqc0, body0)

data_system = create_data_system(mech.eqconstraints.values,
    mech.bodies.values, mech.ineqconstraints.values);

∂g∂z(mech, eqc0, body0)
∂i∂z(body0, mech.Δt, attjac=false)
∂g∂z(mech, eqc0, body0) * ∂i∂z(body0, mech.Δt, attjac=false)

eqc_databody(mech, eqc0, body0)


dataineqc!(data_system, mech)
plot(Gray.(1e10 .* abs.(full_matrix(data_system))))

databody!(data_system, mech)
plot(Gray.(1e10 .* abs.(full_matrix(data_system))))

dataeqc!(data_system, mech)
plot(Gray.(1e10 .* abs.(full_matrix(data_system))))


full_matrix(data_system)

eqcs = mech.eqconstraints.values
bodies = mech.bodies.values
ineqcs = mech.ineqconstraints.values
A = adjacencyMatrix(eqcs, bodies, ineqcs)
plot(Gray.(A))

data_bound!(data_system, mech)
















mech.origin.id
getfield.(mech.eqconstraints.values, :id)
getfield.(mech.bodies.values, :id)
getfield.(mech.ineqconstraints.values, :id)



full_vector(system) = vcat(getfield.(system.vector_entries,:value)...)
mech.system
eqcs = mech.eqconstraints.values
bodies = mech.bodies.values
ineqcs = mech.ineqconstraints.values
A, dims = adjacencyMatrix(eqcs, bodies, ineqcs)


full_vector(system) = vcat(getfield.(system.vector_entries,:value)...)
mech.system.matrix_entries.rowval

full_matrix(mech.system)
