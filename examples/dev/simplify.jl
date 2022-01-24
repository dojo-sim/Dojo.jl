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
datadim(body::Body) = 19 # 1+6+6+6 [m,flat(J),x1,q1,x2,q2] with attjac
# Eqconstraints
datadim(eqc::EqualityConstraint) = 2 + sum(datadim.(eqc.constraints)) # [spring, damper, utra, urot, tra_spring_offset, rot_spring_offset]
datadim(joint::Rotational{T,Nλ,Nb,N,Nb½,N̄λ}) where {T,Nλ,Nb,N,Nb½,N̄λ} = 2N̄λ # [u, spring, damper, spring_offset]
datadim(joint::Translational{T,Nλ,Nb,N,Nb½,N̄λ}) where {T,Nλ,Nb,N,Nb½,N̄λ} = 2N̄λ # [u, spring, damper, spring_offset]
# Ineqconstraints
datadim(ineqc::InequalityConstraint) = sum(datadim.(ineqc.constraints))
datadim(bound::ContactBound) = 7 # [cf, p, offset]
datadim(bound::LinearContactBound) = 7 # [cf, p, offset]
datadim(bound::ImpactBound) = 6 # [p, offset]

function lift_inertia(j::SVector{6,T}) where T
    J = SMatrix{3,3,T,9}(
        [j[1] j[2] j[3];
         j[2] j[4] j[5];
         j[3] j[5] j[6]])
end
function flatten_inertia(J::SMatrix{3,3,T,9}) where T
    j = SVector{6,T}([J[1,1], J[1,2], J[1,3], J[2,2], J[2,3], J[3,3]])
end
function ∂inertia(p)
    SA[
        p[1]  p[2]  p[3]  0     0     0;
        0     p[1]  0     p[2]  p[3]  0;
        0     0     p[1]  0     p[2]  p[3];
    ]
end

function eqc_databody(mechanism::Mechanism, eqc::EqualityConstraint{T,N},
        body::Body{T}) where {T,N}
    Nd = datadim(body)
    ∇m = szeros(T,N,1)
    ∇J = szeros(T,N,6)
    ∇z1 = szeros(T,N,6)
    ∇z2 = ∂g∂z(mechanism, eqc, body) * ∂i∂z(body, mechanism.Δt, attjac=true)
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
    ∇J = -4 / Δt * LVᵀmat(q2)' * LVᵀmat(q1) * ∂inertia(VLᵀmat(q1) * vector(q2))
    ∇J += -4 / Δt * LVᵀmat(q2)' * Tmat() * RᵀVᵀmat(q3) * ∂inertia(VLᵀmat(q2) * vector(q3))
    ∇J = [szeros(T,3,6); ∇J]

    ∇x1 = 1 / Δt * body.m * SMatrix{3,3,T,9}(Diagonal(sones(T,3)))
    ∇q1 = -4 / Δt * LVᵀmat(q2)' * ∂qLVᵀmat(body.J * VLᵀmat(q1) * vector(q2))
    ∇q1 += -4 / Δt * LVᵀmat(q2)' * LVᵀmat(q1) * body.J * ∂qVLᵀmat(vector(q2))
    ∇q1 *= LVᵀmat(q1) # attjac
    ∇z1 = [∇x1 szeros(T,3,3);
           szeros(T,3,3) ∇q1]

    # TODO
    ∇z2 = szeros(T,N,6)

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
    ∇compμ = szeros(T,N½,Nd)
    ∇m = szeros(T,N½,1)
    ∇J = szeros(T,N½,6)
    ∇z1 = szeros(T,N½,6)
    ∇z3 = ∂g∂z(mechanism, ineqc, body)
    ∇z2 = ∇z3 * ∂i∂z(body, mechanism.Δt, attjac=true) # 4x7 * 7x6 = 4x6
    ∇g = [∇m ∇J ∇z1 ∇z2]
    return [∇compμ; ∇g]
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

mech = getsnake(Nb=3, damper=0.0, spring=0.0, contact_type=:contact);
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

∇0 = eqc_dataeqc(mech, eqc0)
∇0 = eqc_databody(mech, eqc0, body0)
∇0 = ineqc_databody(mech, ineqc0, body0)
∇0 = ineqc_dataineqc(mech, ineqc0, body0)
∇0 = body_dataeqc(mech, eqc0, body0)
∇0 = body_databody(mech, body0)
∇0 = body_dataineqc(mech, ineqc0, body0)

data_system = create_data_system(mech.eqconstraints.values,
    mech.bodies.values, mech.ineqconstraints.values);

∂g∂z(mech, eqc0, body0)
∂i∂z(body0, mech.Δt, attjac=false)
∂g∂z(mech, eqc0, body0) * ∂i∂z(body0, mech.Δt, attjac=false)

eqc_databody(mech, eqc0, body0)
sum(ones(10,5), dims=2)

dataineqc!(data_system, mech)
plot(Gray.(1e10 .* abs.(full_matrix(data_system))))

databody!(data_system, mech)
plot(Gray.(1e10 .* abs.(full_matrix(data_system))))

dataeqc!(data_system, mech)
plot(Gray.(1e10 .* abs.(full_matrix(data_system))))
plot(log.(10, abs.(sum(full_matrix(data_system), dims=1)[1,:])))


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

















ineqc0
bound0 = ineqc0.constraints[1]
λ = nothing
s = SVector{4}([1.0; rand(3)/10])
γ = SVector{4}([1.0; rand(3)/10])
g(bound0, s, γ, x3, q3, v25, ϕ25)
jacv0 = ∂g∂v(bound0, x3, q3, x2, v25, q2, ϕ25, λ, Δt)
jacv1 = FiniteDiff.finite_difference_jacobian(
    v -> g(bound0, s, γ, getx3(x2, SVector{3}(v[1:3]), Δt),
        getq3(q2, SVector{3}(v[4:6]), Δt), SVector{3}(v[1:3]), SVector{3}(v[4:6])), [v25; ϕ25])

jacz0 = ∂g∂z(bound0, x3, q3, x2, v25, q2, ϕ25, λ, Δt)
jacz1 = FiniteDiff.finite_difference_jacobian(
    z -> g(bound0, s, γ, SVector{3}(z[1:3]),
        UnitQuaternion(z[4:7],false), SVector{3}(v[1:3]), SVector{3}(v[4:6])), [x3; vector(q3)])

norm(jacv0 - jacv1)
jacv0 - jacv1
norm(jacz0 - jacz1)
jacz0 - jacz1



# using Symbolics
# @variables j[1:6]
# @variables p[1:3]
# J = SMatrix{3,3,Num,9}([j[1] j[2] j[3];
#      j[2] j[4] j[5];
#      j[3] j[5] j[6]])
# Symbolics.jacobian(J * p, j)
# lift_inertia(flatten_inertia(J)) - J
# flatten_inertia(J)
