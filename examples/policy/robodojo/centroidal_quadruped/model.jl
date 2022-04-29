"""
    centroidal quadruped
    q = (p, r, f1, f2, f3, f4)
        p - body position
        r - body orientation (modified Rodriques parameters)
        f1 - foot 1 position
        f2 - foot 2 position
        f3 - foot 3 position
        f4 - foot 4 position
"""
mutable struct CentroidalQuadruped114{T} <: RoboDojo.Model{T}
    # dimensions
	nq::Int # generalized coordinates
    nu::Int # controls
    nw::Int # parameters
    nc::Int # contact points


	# parameters
	body_height::T
	foot_x::T
	foot_y::T

	# parameters
    mass_body::T
    inertia_body::Matrix{T}
    mass_foot::T

	# joint friction
	friction_joint::Vector{T}
	spring_stiffness_joint::Vector{T}

	# environment
    friction_foot_world::Vector{T}
    gravity::T
end

function L_mult(x)
    [x[1] -transpose(x[2:4]);
     x[2:4] x[1] * I(3) + skew(x[2:4])]
end

# right quaternion multiply as matrix
function R_mult(x)
    [x[1] -transpose(x[2:4]); x[2:4] x[1] * I(3) - skew(x[2:4])]
end

# rotation matrix
function quat_rot_mat(q)
    H = [zeros(1, 3); I(3)]
    transpose(H) * L_mult(q) * transpose(R_mult(q)) * H
end

function quat_from_mrp(p)
    """Quaternion (scalar first) from MRP"""
    return (1.0 / (1.0 + dot(p, p))) * [(1 - dot(p, p)); 2.0 * p]
end

function mrp_rot_mat(x)
    quat_rot_mat(quat_from_mrp(x))
end

# Kinematics
function kinematics(model::CentroidalQuadruped114, q)
	q[6 .+ (1:12)]
end

lagrangian(model::CentroidalQuadruped114, q, q̇) = 0.0

function mass_matrix(model::CentroidalQuadruped114, q)
    cat(
        model.mass_body * Diagonal(ones(3)),     # body position
        model.inertia_body,                      # body orienation
        model.mass_foot * Diagonal(ones(3 * 4)), # feet position
        dims=(1, 2)
        )
end

function dynamics_bias(model::CentroidalQuadruped114, q, q̇)
	offsets = [
		[+model.foot_x, +model.foot_y, -model.body_height],
		[+model.foot_x, -model.foot_y, -model.body_height],
		[-model.foot_x, +model.foot_y, -model.body_height],
		[-model.foot_x, -model.foot_y, -model.body_height],
		]

	spring = model.spring_stiffness_joint
    [
        model.mass_body * [0,0,model.gravity] + spring .* (+sum(offsets) + 4 * q[1:3] - q[7:9] - q[10:12] - q[13:15] - q[16:18]);     # body position
        skew(q̇[4:6]) * model.inertia_body * q̇[4:6]; # body orienation
        model.mass_foot * [0,0,model.gravity] + spring .* (-offsets[1] + q[7:9] - q[1:3]);
        model.mass_foot * [0,0,model.gravity] + spring .* (-offsets[2] + q[10:12] - q[1:3]);
        model.mass_foot * [0,0,model.gravity] + spring .* (-offsets[3] + q[13:15] - q[1:3]);
        model.mass_foot * [0,0,model.gravity] + spring .* (-offsets[4] + q[16:18] - q[1:3]);
    ]
end

function signed_distance(model::CentroidalQuadruped114, q)

    position_foot1 = q[6 .+ (1:3)]
    position_foot2 = q[9 .+ (1:3)]
    position_foot3 = q[12 .+ (1:3)]
	position_foot4 = q[15 .+ (1:3)]

	# return [q[3]; position_foot1[3]; position_foot2[3]; position_foot3[3]; position_foot4[3]]
	return [position_foot1[3]; position_foot2[3]; position_foot3[3]; position_foot4[3]]
end

function input_jacobian(model::CentroidalQuadruped114, q)
    position_body = q[1:3]
    orientation_body = q[3 .+ (1:3)]
	R = mrp_rot_mat(orientation_body)
	# R = euler_rotation_matrix(orientation_body)

	# kinematics in world frame
	r1 = q[6 .+ (1:3)] - position_body
	r2 = q[9 .+ (1:3)] - position_body
	r3 = q[12 .+ (1:3)] - position_body
	r4 = q[15 .+ (1:3)] - position_body

	z3 = zeros(3, 3)

	transpose([
        I(3) z3   I(3) I(3) I(3) I(3);
		z3   I(3) transpose(R) * skew(r1) transpose(R) * skew(r2) transpose(R) * skew(r3) transpose(R) * skew(r4);
        z3   z3   -I(3)    z3    z3   z3;
        z3   z3   z3    -I(3)    z3   z3;
        z3   z3   z3       z3 -I(3)   z3;
        z3   z3   z3       z3    z3 -I(3)
    ])
end

function disturbance_jacobian(model::CentroidalQuadruped114, q)
    @SMatrix [1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
              0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
			  0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
end

function contact_jacobian(model::CentroidalQuadruped114, q)
    z3 = zeros(3, 3)

    [
        # I(3) z3   z3   z3   z3   z3;
        z3   z3 I(3)   z3   z3   z3;
		z3   z3   z3 I(3)   z3   z3;
        z3   z3   z3   z3 I(3)   z3;
        z3   z3   z3   z3   z3 I(3);
    ]
end

function RD.nominal_configuration(model::CentroidalQuadruped114)
    [
        0.0; 0.0; model.body_height;
        0.0; 0.0; 0.0;
        model.foot_x; model.foot_y; 0.0;
        model.foot_x;-model.foot_y; 0.0;
       -model.foot_x; model.foot_y; 0.0;
       -model.foot_x;-model.foot_y; 0.0;
    ]
end

function nominal_state(model::CentroidalQuadruped114)
	[nominal_configuration(model); zeros(model.nq)]
end

# friction coefficients
function RD.friction_coefficients(model::CentroidalQuadruped114)
    model.friction_foot_world
end

# dimensions
nq = 3 + 3 + 3 * 4       # generalized coordinates
nu = 6 + 3 * 4           # controls
nw = 3                   # parameters
nc = 4                   # contact points

# inertial properties
scale = 10.0
mass_body = 13.5 / scale
i_xx = 0.0178533 / scale
i_xy = 0.0 / scale
i_xz = 0.0 / scale
i_yz = 0.0 / scale
i_yy = 0.0377999 / scale
i_zz = 0.0456542 / scale
inertia_body = Array(Diagonal([i_xx, i_yy, i_zz]))
mass_foot = 0.2

# parameters
body_height = 0.3
foot_x = 0.17
foot_y = 0.15

friction_joint = 0.15 * ones(3)
spring_stiffness_joint = [5,5,15.0]
friction_foot_world = 0.3*ones(nc)     # coefficient of friction
gravity = 9.81                   # gravity

centroidal_quadruped = CentroidalQuadruped114(nq, nu, nw, nc,
				body_height,
				foot_x,
				foot_y,
				mass_body,
				inertia_body,
				mass_foot,
				friction_joint,
				spring_stiffness_joint,
				friction_foot_world,
				gravity,
                )


floating_base_dim(::CentroidalQuadruped114) = 6
name(::CentroidalQuadruped114) = :centroidal_quadruped
