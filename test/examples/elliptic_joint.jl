using ConstrainedDynamics
using ConstrainedDynamics: GenericJoint,Vmat,params
using StaticArrays


# Parameters
length1 = 1.0
width, depth = 0.1, 0.1

# Links
origin = Origin{Float64}()
link1 = Box(width, depth, length1, length1)


@inline function g(joint::GenericJoint, xb::AbstractVector, qb::UnitQuaternion)
    a=5
    b=2

    if xb[3]==0
        xb[2] > 0 ? α = -pi/2 : α = pi/2
    elseif xb[3] > 0
        α = -atan(b^2/a^2*xb[2]/xb[3])
    else
        α = -atan(b^2/a^2*xb[2]/xb[3]) - pi
    end

    qα = UnitQuaternion(cos(α/2),sin(α/2),0,0,false)

    eqc1=Vmat(qb\qα)
    eqc2=SA[xb[1]; (xb[2]^2/a^2)+(xb[3]^2/b^2)-1]
    G= [eqc1;eqc2]
    return G
end

# Constraints
joint_between_origin_and_link1 = EqualityConstraint(GenericJoint{5}(origin,link1,g))


links = [link1]
constraints = [joint_between_origin_and_link1]


mech = Mechanism(origin, links, constraints)
setPosition!(link1,x=[0;0.0;2],q = UnitQuaternion(RotX(0)))
setVelocity!(link1,v=[0;2.0;0])
