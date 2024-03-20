@setup_workload begin
    @compile_workload begin 
        # One full simulation
        origin = Origin()
        body = Cylinder(0.1, 1, 1)
        joint = JointConstraint(Revolute(origin, body, [1;0;0]; child_vertex=[0;0;1/2]))
        mechanism = Mechanism(origin, [body], [joint])
        controller!(mechanism, k) = set_input!(joint, [0])
        set_minimal_coordinates!(mechanism, joint, [pi/4])
        set_minimal_velocities!(mechanism, joint, [0.2])
        storage = simulate!(mechanism, mechanism.timestep * 2, controller!, record=true)

        # Common shapes
        mesh = Mesh("", 1, rand(3,3))
        box = Box(1,1,1,1.0)
        capsule = Capsule(1,1,1.0)
        cylinder = Cylinder(1, 1, 1.0)
        sphere = Sphere(1,1.0)
        pyramid = Pyramid(1,1,1.0)

        # Common joints
        joint_axis = [1;0;0]
        JointConstraint(Floating(origin, mesh))
        JointConstraint(Fixed(origin, box))
        JointConstraint(Prismatic(origin, capsule, joint_axis))
        JointConstraint(Planar(origin, cylinder, joint_axis))
        JointConstraint(FixedOrientation(origin, sphere))
        JointConstraint(Revolute(origin, pyramid, joint_axis))
        JointConstraint(Cylindrical(origin, mesh, joint_axis))
        JointConstraint(PlanarAxis(origin, box, joint_axis))
        JointConstraint(FreeRevolute(origin, capsule, joint_axis))
        JointConstraint(Orbital(origin, cylinder, joint_axis))
        JointConstraint(PrismaticOrbital(origin, sphere, joint_axis))
        JointConstraint(PlanarOrbital(origin, pyramid, joint_axis))
        JointConstraint(FreeOrbital(origin, mesh, joint_axis))
        JointConstraint(Spherical(origin, box))
        JointConstraint(CylindricalFree(origin, capsule, joint_axis))
        JointConstraint(PlanarFree(origin, cylinder, joint_axis))
    end
end