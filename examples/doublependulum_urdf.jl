using ConstrainedDynamics
using ConstrainedDynamicsVis


path = "examples/examples_files/doublependulum.urdf"
mech = Mechanism(path)

storage = simulate!(mech, 10., record = true)
visualize(mech, storage)
