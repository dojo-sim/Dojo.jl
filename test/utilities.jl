@test length(Dojo.module_dir()) > 0
@test Dojo.scn(1.0) == " 1.0e+0"
@test Dojo.scn(0.0) == " 0.0e+0"
@test Dojo.scn(Inf) == " Inf"
@test Dojo.scn(123.1, 
    digits=0) == " 1e+2"


# show functions
test_origin = Origin{Float64}()
display(test_origin)
@test true

test_body = Box(1.0,2.0,3.0,4.0)
display(test_body)
display(test_body.shape)
display(test_body.state)
@test true

test_joint = JointConstraint(Fixed(test_origin,test_body))
display(test_joint)
display(test_joint.translational)
display(test_joint.rotational)
@test true

test_mechanism = Mechanism(test_origin,[test_body],[test_joint])
display(test_mechanism)
# display(test_mechanism.system) # TODO Removed from tests because '_show_with_braille_patterns' fails during tests for Julia 1.6, but it works in GBS for 1.6. 
@test true

test_storage = Storage(3,4)
display(test_storage)
@test true

test_options = SolverOptions()
display(test_options)
@test true
