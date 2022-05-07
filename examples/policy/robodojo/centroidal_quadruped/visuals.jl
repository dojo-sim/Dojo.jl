function build_robot!(vis::Visualizer, model::CentroidalQuadruped114;
	body_height=0.025,
	body_length=0.17,
	body_width=0.15,
	foot_radius=0.025,
	color_opacity=1.0,
	name::Symbol=:robot)

	body_mat = MeshPhongMaterial(color=RGBA(0.0, 0.0, 0.0, color_opacity))
	foot_mat = MeshPhongMaterial(color=RGBA(1.0, 165.0 / 255.0, 0.0, color_opacity))

	default_background!(vis)

	setobject!(vis[name][:body],
    	GeometryBasics.Rect(Vec(-body_length, -body_width, -body_height),Vec(2.0 * body_length, 2.0 * body_width, body_height)),
    	body_mat)

	foot1 = setobject!(vis[name][:foot1], GeometryBasics.Sphere(Point3f0(0),
        convert(Float32, foot_radius)),
        foot_mat)

	foot2 = setobject!(vis[name][:foot2], GeometryBasics.Sphere(Point3f0(0),
		convert(Float32, foot_radius)),
		foot_mat)

	foot3 = setobject!(vis[name][:foot3], GeometryBasics.Sphere(Point3f0(0),
		convert(Float32, foot_radius)),
		foot_mat)

	foot4 = setobject!(vis[name][:foot4], GeometryBasics.Sphere(Point3f0(0),
		convert(Float32, foot_radius)),
		foot_mat)
end

function set_robot!(vis::Visualizer, model::CentroidalQuadruped114, q::AbstractVector;
	p_shift=[0.0; 0.0; 0.025], name::Symbol=:robot)

	R = mrp_rotation_matrix(q[4:6])

	p_body  = q[0  .+ (1:3)] + p_shift
	p_foot1 = q[6  .+ (1:3)] + p_shift
	p_foot2 = q[9  .+ (1:3)] + p_shift
	p_foot3 = q[12 .+ (1:3)] + p_shift
	p_foot4 = q[15 .+ (1:3)] + p_shift

	settransform!(vis[name][:body], compose(Translation(p_body), LinearMap(mrp_rotation_matrix(q[4:6]))))
	settransform!(vis[name][:foot1], Translation(p_foot1))
	settransform!(vis[name][:foot2], Translation(p_foot2))
	settransform!(vis[name][:foot3], Translation(p_foot3))
	settransform!(vis[name][:foot4], Translation(p_foot4))

end

function visualize!(vis, model::CentroidalQuadruped114, q;
	Δt=0.1,
	body_height=0.025,
	body_length=0.17,
	body_width=0.15,
	foot_radius=0.025,
	color_opacity=1.0,
	fixed_camera=false,
	anim=MeshCat.Animation(convert(Int, floor(1.0 / Δt))),
	name::Symbol=:robot)

	build_robot!(vis, model,
		body_height=body_height,
		body_length=body_length,
		body_width=body_width,
		foot_radius=foot_radius,
		color_opacity=color_opacity,
		name=name)



	for (t, qt) in enumerate(q)
		MeshCat.atframe(anim, t) do
			set_robot!(vis, model, qt, name=name)
		end
	end

	MeshCat.setanimation!(vis, anim)

	return anim
end
