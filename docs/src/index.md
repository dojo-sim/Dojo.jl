# Get Started 

__[Dojo](https://github.com/dojo-sim/Dojo.jl) is a differentiable simulator for robotics__, prioritizing accurate physics and useful gradients. The simulator is written in pure Julia in order to be both performant and easy to use.

## Features
* __Maximal-Coordinates Representation__: Fast and efficient conversion between [maximal](background_representations/maximal_representation.md) and [minimal](background_representations/minimal_representation.md) representations
* __Smooth Gradients__: Simulation with [hard contact](background_contact/impact.md) and useful [gradients](background_representations/gradients.md) through contact events
* __Open Source__: Code is available on [GitHub](https://github.com/dojo-sim/Dojo.jl) and distributed under the MIT License
* __Python Interface__: [dojopy](https://github.com/dojo-sim/dojopy)

## Installation
Dojo can be installed using the Julia package manager for Julia `v1.6` and higher. Inside the Julia REPL, type `]` to enter the Pkg REPL mode then run

`pkg> add Dojo`

## Related talks


## Credits

The following people are involved in the development of Dojo:

__Primary Development__
* [Simon Le Cleac'h](https://simon-lc.github.io/) (main development, contact modeling, interior-point solver, gradients)
* [Taylor Howell](https://thowell.github.io/) (main development, contact modeling, interior-point solver, gradients)
* [Jan Bruedigam](https://github.com/janbruedigam) (main development, maximal representation and graph-based solver)


* [Zico Kolter](https://zicokolter.com/)
* [Mac Schwager](https://web.stanford.edu/~schwager/)
* [Zachary Manchester](https://www.ri.cmu.edu/ri-faculty/zachary-manchester/) (principal investigator)

__Additional Contributions__
* [Suvansh Sanjeev](https://suvan.sh/) (PyTorch interface)
* [Benjamin Bokser](http://www.benbokser.com/) (REx Hopper)

Development by the [Robotic Exploration Lab](https://roboticexplorationlab.org/).
 
If this project is useful for your work please consider
* [Citing](citing.md) the relevant papers
* Leaving a star on the [GitHub repository](https://github.com/dojo-sim/Dojo.jl)

## Licence
Dojo.jl is licensed under the MIT License. For more details click [here](https://github.com/dojo-sim/Dojo.jl/blob/main/LICENSE.md).
