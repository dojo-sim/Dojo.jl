using Dojo
import DifferentiableCollisions as dc
using LinearAlgebra
using Printf
using StaticArrays
import FiniteDiff as FD2
import MeshCat as mc
import Random
using Colors
using SparseArrays
using Combinatorics


function Gbar_2_body(state1, state2)
    # attitude jacobian for two rigid bodies
    Gbar1 = blockdiag(sparse(I(3)),sparse(Dojo.Lmat(state1.q2)*Dojo.Vᵀmat()))
    Gbar2 = blockdiag(sparse(I(3)),sparse(Dojo.Lmat(state2.q2)*Dojo.Vᵀmat()))
    Gbar = Matrix(blockdiag(Gbar1,Gbar2))
end

# modify the torque to account for the extra 2x from quaternion stuff
const τ_mod = Diagonal(kron(ones(2),[ones(3);0.5*ones(3)]))

function Dojo.impulse_map(relative::Symbol, model::Dojo.Contact, pbody::Dojo.Node, cbody::Dojo.Node)
    pbody.shape.primitive.r = pbody.state.x2
    pbody.shape.primitive.q = Dojo.vector(pbody.state.q2)
    cbody.shape.primitive.r = cbody.state.x2
    cbody.shape.primitive.q = Dojo.vector(cbody.state.q2)

    D_α = reshape(dc.proximity_gradient(model.collision.primitive1,model.collision.primitive2)[2], 1, 14)
    E = h * τ_mod * (D_α*Gbar_2_body(pbody.state, cbody.state))'

    if relative == :parent 
        return E[1:6]
    elseif relative == :child 
        return E[7:12]
    end
end

# time step
const h = 0.05


Random.seed!(1)

origin = Origin()
body1 = Cylinder(0.5,3.0,1;orientation_offset=Dojo.RotY(pi/2))
body1.inertia = diagm([body1.inertia[3,3];body1.inertia[2,2];body1.inertia[1,1]])
body2 = Sphere(1.0,1)
body3 = Pyramid(1,1.5,1)
body4 = Box(1,1,1,1)
bodies = [body1;body2;body3;body4]
joint1 = JointConstraint(Floating(origin, body1))
joint2 = JointConstraint(Floating(origin, body2))
joint3 = JointConstraint(Floating(origin, body3))
joint4 = JointConstraint(Floating(origin, body4))
joints = [joint1;joint2;joint3;joint4]
dojo_contacts = ContactConstraint(ImpactContact(bodies))

mech = Mechanism(origin, bodies, joints, dojo_contacts; timestep=h)

# create the indexing named tuple
N_bodies = length(bodies) 

# initial conditions of everything
rs = [5*(@SVector randn(3)) for i = 1:N_bodies]
rs[end] = SA[0,0,0.0]
qs = [SA[1,0,0,0.0] for i = 1:N_bodies]

# initial velocities
vs =  [SA[1,1,1.0] for i = 1:N_bodies]
for i = 1:N_bodies
    vs[i] = -.5*rs[i]
end
ωs = [deg2rad.(20*(@SVector randn(3))) for i = 1:N_bodies]

for (i,joint) in enumerate(mech.joints)
    set_minimal_coordinates!(mech, joint, [rs[i];0;0;0])
    set_minimal_velocities!(mech, joint, [vs[i];ωs[i]])
end

storage = simulate!(mech,79*0.05;record=true)

visualize(mech, storage; visualize_floor=false, framerate=Inf)
