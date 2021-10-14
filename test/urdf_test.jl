using ConstrainedDynamics

Mechanism("urdf_test.urdf", floating=false)
@test true

Mechanism("urdf_test.urdf", floating=true)
@test true
