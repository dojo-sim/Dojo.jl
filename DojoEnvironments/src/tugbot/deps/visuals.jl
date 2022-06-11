function tugbot_visualizer!(vis::Visualizer, mech::Mechanism;)
    material = MeshPhongMaterial(color=Colors.RGBA(0.0, 0.0, 0.0, 1.0))
    drone = MeshFileObject(joinpath(dirname(pathof(DojoEnvironments)), "tugbot/deps/drone.obj"))
    setobject!(vis["robot"]["bodies"]["drone__id_1"][:mesh], drone)
    settransform!(vis["robot"]["bodies"]["drone__id_1"][:mesh], MeshCat.LinearMap(RotX(π/2)))
end
