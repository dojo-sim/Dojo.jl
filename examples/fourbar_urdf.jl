using ConstrainedDynamics
using ConstrainedDynamicsVis


path = "examples/examples_files/fourbar.urdf"
mech = Mechanism(path)

initializeConstraints!(mech, newtonIter=1000)

steps = Base.OneTo(1000)
storage = Storage{Float64}(steps,4)

storage = simulate!(mech, storage, record = true)
visualize(mech, storage)
