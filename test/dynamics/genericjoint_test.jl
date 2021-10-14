using ConstrainedDynamics
using ConstrainedDynamics: vrotate, Vmat, GenericJoint, AbstractBody, orthogonalrows, deepcopy
using LinearAlgebra


function Fixed_abstract(body1::AbstractBody{T}, body2; p1 = szeros(T, 3), p2 = szeros(T, 3), qoffset = one(UnitQuaternion{T})) where T  
    vertices = (p1, p2)
    function gfixed(joint::GenericJoint, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
       G= [
           vrotate(xb + vrotate(vertices[2], qb) - (xa + vrotate(vertices[1], qa)), inv(qa));
           Vmat(qa \ qb / qoffset)
       ]

       return G
   end
   function gfixed(joint::GenericJoint, xb::AbstractVector, qb::UnitQuaternion)
       G= [
           xb + vrotate(vertices[2], qb) - vertices[1];
           Vmat(qb / qoffset)
       ]

       return G
   end

   return GenericJoint{6}(body1, body2, gfixed) 
end
function Revolute_abstract(body1::AbstractBody{T}, body2, axis; p1 = szeros(T, 3), p2 = szeros(T, 3), qoffset = one(UnitQuaternion{T})) where T  
    vertices = (p1, p2)
    V1, V2, V3 = orthogonalrows(axis)  
    V12 = [V1;V2]
    
    function grevolute(joint::GenericJoint, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
       G= [
           vrotate(xb + vrotate(vertices[2], qb) - (xa + vrotate(vertices[1], qa)), inv(qa));
           V12*Vmat(qa \ qb / qoffset)
       ]

       return G
   end
   function grevolute(joint::GenericJoint, xb::AbstractVector, qb::UnitQuaternion)
       G= [
           xb + vrotate(vertices[2], qb) - vertices[1];
           V12*Vmat(qb / qoffset)
       ]

       return G
   end

   return GenericJoint{5}(body1, body2,grevolute)
end



joint_axis = [1.0;0.0;0.0]
l = 1.0
w, d = .1, .1

vert11 = [0.;0.;l / 2]
vert12 = -vert11
vert21 = [0.;0.;l / 2]

origin = Origin{Float64}()
origin_abstract = deepcopy(origin)
link1 = Box(w, d, l, l)
link1_abstract = deepcopy(link1)
link2 = deepcopy(link1)
link2_abstract = deepcopy(link1)

joint0to1 = EqualityConstraint(Fixed(origin, link1; p2=vert11))
joint1to2 = EqualityConstraint(Revolute(link1, link2,joint_axis; p1=vert12, p2=vert21))
joint0to1_abstract = EqualityConstraint(Fixed_abstract(origin_abstract, link1_abstract; p2=vert11))
joint1to2_abstract = EqualityConstraint(Revolute_abstract(link1_abstract, link2_abstract,joint_axis; p1=vert12, p2=vert21))

links = [link1;link2]
links_abstract = [link1_abstract;link2_abstract]
constraints = [joint0to1;joint1to2]
constraints_abstract = [joint0to1_abstract;joint1to2_abstract]

mech = Mechanism(origin, links, constraints)
mech_abstract= Mechanism(origin_abstract, links_abstract, constraints_abstract)

setPosition!(origin,link1,p2 = vert11)
setPosition!(origin_abstract,link1_abstract,p2 = vert11)

function compare(s1,s2)
    @test all(isapprox.(norm.(s1.x[1]-s2.x[1]), 0.0; atol = 1e-8)) &&
        all(isapprox.(norm.(s1.q[1]./s2.q[1]), 1.0; atol = 1e-8)) &&
        all(isapprox.(norm.(s1.x[2]-s2.x[2]), 0.0; atol = 1e-8)) &&
        all(isapprox.(norm.(s1.q[2]./s2.q[2]), 1.0; atol = 1e-8)) &&
        all(isapprox.(norm.(s1.ω[1]-s2.ω[1]), 0.0; atol = 1e-8)) &&
        all(isapprox.(norm.(s1.v[1]-s2.v[1]), 0.0; atol = 1e-8)) &&
        all(isapprox.(norm.(s1.ω[2]-s2.ω[2]), 0.0; atol = 1e-8)) &&
        all(isapprox.(norm.(s1.v[2]-s2.v[2]), 0.0; atol = 1e-8))
end

for i=1:10
    randang = pi/4*rand()
    setPosition!(link1,link2,p1 = vert12,p2 = vert21,Δq = UnitQuaternion(RotX(randang)))
    setPosition!(link1_abstract,link2_abstract,p1 = vert12,p2 = vert21,Δq = UnitQuaternion(RotX(randang)))
    storage = simulate!(mech, 1, record = true)
    storage_abstract = simulate!(mech_abstract, 1, record = true)

    compare(storage,storage_abstract)
end