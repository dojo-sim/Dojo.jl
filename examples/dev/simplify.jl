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
datadim(eqc::EqualityConstraint) = sum(datadim.(eqc.constraints))
datadim(joint::Rotational{T,Nλ,Nb,N,Nb½,N̄λ}) where {T,Nλ,Nb,N,Nb½,N̄λ} = 2 + N̄λ # [spring, damper, spring_offset]
datadim(joint::Translational{T,Nλ,Nb,N,Nb½,N̄λ}) where {T,Nλ,Nb,N,Nb½,N̄λ} = 2 + N̄λ # [spring, damper, spring_offset]
# Ineqconstraints
datadim(ineqc::InequalityConstraint) = sum(datadim.(ineqc.constraints))
datadim(bound::ContactBound) = 7 # [cf, p, offset]
datadim(bound::LinearContactBound) = 7 # [cf, p, offset]
datadim(bound::ImpactBound) = 6 # [p, offset]

mech = getsnake(Nb=3);
eqcs = mech.eqconstraints.values
bodies = mech.bodies.values
ineqcs = mech.ineqconstraints.values
A = adjacencyMatrix(eqcs, bodies, ineqcs)
plot(Gray.(A))

data_system = create_data_system(eqcs, bodies, ineqcs);
full_matrix(data_system)
data_system.matrix_entries

∂g∂data(bound::Bound, body::Body, id, λ, Δt) = ∂g∂data(bound, posargs3(body.state, Δt)..., fullargssol(body.state)..., λ, Δt)
function ∂g∂data()

end




storage = simulate!(mech, 2.3, ctrl!, record=true, verbose=false)
visualize(mech, storage, vis=vis)
mech.system.matrix_entries[5,6]


data_system = System(A, dims);














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
