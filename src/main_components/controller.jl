"""
$(TYPEDEF)

Users can create subtypes of `Controller` to implement controllers. 
If such a controller object is passed to `simulate!`, at each time step, `controller.control!(mechanism, controller, timestep)` will be called. 
This way, the user has access to all attributes of the `mechanism`, the `controller`, and the time step. Note that the control function must be called `control!`.
"""
abstract type Controller end