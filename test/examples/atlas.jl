using ConstrainedDynamics


path = "examples/examples_files/atlas_simple.urdf"
mech = Mechanism(path, floating=false, g = -.5)
