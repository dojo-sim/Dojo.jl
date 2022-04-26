function visualize_drone!(vis, model::Drone, Q, U;
    Δt=0.1,
    r=0.05,
    xT=[zero(q[end]) for q in Q],
    color=RGBA(0,0,0,1.0),
    color_goal=RGBA(0,1,0,1.0),
    color_thrust=RGBA(1.0, 0.0, 0.0, 1.0))

    N = length(Q)
    T = length(Q[1])
    @assert all([length(Q[i]) == length(U[i]) for i = 1:N])

    setvisible!(vis["/Background"], true)
	setprop!(vis["/Background"], "top_color", RGBA(1.0, 1.0, 1.0, 1.0))
	setprop!(vis["/Background"], "bottom_color", RGBA(1.0, 1.0, 1.0, 1.0))
	setvisible!(vis["/Axes"], false)
	setvisible!(vis["/Grid"], false)

    body = MeshCat.Cylinder(Point3f0(-0.5 * model.body_length, 0.0, 0.0), Point3f0(0.5 * model.body_length, 0.0, 0.0), convert(Float32, r))
    thruster1 = MeshCat.Cylinder(Point3f0(0.0, 0.0, -0.25 * model.body_length), Point3f0(0.0, 0.0, 0.15 * model.body_length), convert(Float32, 0.65 * r))
    thruster2 = MeshCat.Cylinder(Point3f0(0.0, 0.0, -0.25 * model.body_length), Point3f0(0.0, 0.0, 0.15 * model.body_length), convert(Float32, 0.65 * r))
    thrust1 = MeshCat.Cylinder(Point3f0(0.0, 0.0, -0.3 * model.body_length), Point3f0(0.0, 0.0, 0.0), convert(Float32, 0.35 * r))
    thrust2 = MeshCat.Cylinder(Point3f0(0.0, 0.0, -0.3 * model.body_length), Point3f0(0.0, 0.0, 0.0), convert(Float32, 0.35 * r))

    sphere_goal = GeometryBasics.Sphere(Point(0.0,0.0,0.0), 1.5 * r)

    for i = 1:N
        setobject!(vis["goal_$i"], sphere_goal, MeshPhongMaterial(color=color_goal))
        settransform!(vis["goal_$i"], Translation([xT[i][1]; 0.1; xT[i][2]]))
        setobject!(vis["drone_$i"]["body"], body, MeshPhongMaterial(color=color))
        setobject!(vis["drone_$i"]["thrust1"], thrust1, MeshPhongMaterial(color=color_thrust))
        setobject!(vis["drone_$i"]["thrust2"], thrust2, MeshPhongMaterial(color=color_thrust))
        setobject!(vis["drone_$i"]["thruster1"], thruster1, MeshPhongMaterial(color=color))
        setobject!(vis["drone_$i"]["thruster2"], thruster2, MeshPhongMaterial(color=color))
    end

    anim = MeshCat.Animation(convert(Int, floor(1.0 / Δt)))

    for t = 1:T
        MeshCat.atframe(anim,t) do
            for i = 1:N
                settransform!(vis["drone_$i"], compose(Translation([Q[i][t][1]; 0.0; Q[i][t][2]]), LinearMap(RotY(-Q[i][t][3]))))
                settransform!(vis["drone_$i"]["thrust1"], compose(Translation([0.5 * model.body_length; 0.0; 0.0]), LinearMap(RotY(-U[i][t][2]))))
                settransform!(vis["drone_$i"]["thrust2"], compose(Translation([-0.5 * model.body_length; 0.0; 0.0]), LinearMap(RotY(-U[i][t][4]))))
                settransform!(vis["drone_$i"]["thruster1"], compose(Translation([0.5 * model.body_length; 0.0; 0.0]), LinearMap(RotY(-U[i][t][2]))))
                settransform!(vis["drone_$i"]["thruster2"], compose(Translation([-0.5 * model.body_length; 0.0; 0.0]), LinearMap(RotY(-U[i][t][4]))))
            end
        end
    end

    # set camera
    settransform!(vis["/Cameras/default"],
        compose(Translation(0.0, -50.0, -1.0), LinearMap(RotZ(- pi / 2))))
    setprop!(vis["/Cameras/default/rotated/<object>"], "zoom", 20)

    MeshCat.setanimation!(vis,anim)
end
