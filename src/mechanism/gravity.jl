# gravity
get_gravity(g::T) where T <: Real = SVector{3,T}([0.0; 0.0; g])
get_gravity(g::Vector{T}) where T = SVector{3,T}(g)
get_gravity(g::SVector) = g

# compensation
function gravity_compensation(mechanism::Mechanism)
    # only works with revolute joints for now
    nu = control_dimension(mechanism)
    u = zeros(nu)
    off  = 0
    for joint in mechanism.joints
        nu = control_dimension(joint)
        if joint.parent_id != 0
            body = get_body(mechanism, joint.parent_id)
            rot = joint.rotational
            A = Matrix(nullspace_mask(rot))
            input = spring_impulses(mechanism, joint, body)
            F = input[1:3]
            τ = input[4:6]
            u[off .+ (1:nu)] = -A * τ
        else
            @warn "need to treat the joint to origin"
        end
        off += nu
    end
    return u
end