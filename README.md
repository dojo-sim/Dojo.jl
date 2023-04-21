[![CI](https://github.com/dojo-sim/Dojo.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/dojo-sim/Dojo.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/dojo-sim/Dojo.jl/branch/main/graph/badge.svg?token=NMS3JQZ2OE)](https://codecov.io/gh/dojo-sim/Dojo.jl)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://dojo-sim.github.io/Dojo.jl/dev)

# Dojo
A differentiable physics engine for robotics
- arXiv preprint: https://arxiv.org/abs/2203.00806
- Python interface: https://github.com/dojo-sim/dojopy
- site: https://sites.google.com/view/dojo-sim
- video presentation: https://youtu.be/TRtOESXJxJQ

[![IMAGE ALT TEXT](https://i.ytimg.com/vi/TRtOESXJxJQ/hq720.jpg?sqp=-oaymwEcCOgCEMoBSFXyq4qpAw4IARUAAIhCGAFwAcABBg==&rs=AOn4CLD1RdCHZ0Z1zSkv1N-PD0Ds79lDiA)](https://youtu.be/TRtOESXJxJQ "Dojo: A Differentiable Simulator for Robotics")

# Update April 2023
- We are no longer actively developing Dojo, but pull requests are always welcome.
- We have updated or removed examples to account for changes since the initial version of Dojo.
- Additional developments on differentiable simulation:
  - Differentiable collision detection (Kevin Tracy): [capsules](https://arxiv.org/abs/2207.00202), [convex primitives](https://arxiv.org/abs/2207.00669) 
  - Single-level contact dynamics + collision detection (Simon Le Cleac'h): [Silico](https://arxiv.org/pdf/2212.06764.pdf)

# Examples

## Simulation
<p float="left">
	<img src="docs/src/assets/animations/atlas_drop.gif" width="120"/>
	<img src="docs/src/assets//animations/astronaut.gif" width="210"/>
	<img src="docs/src/assets/animations/dzhanibekov.gif" width="180"/>
	<img src="docs/src/assets/animations/tippetop.gif" width="180"/>
</p>

## Learning and Control
<p float="left">
	<img src="docs/src/assets/animations/quadruped.gif" width="350"/>
	<img src="docs/src/assets/animations/ant_ars.gif" width="350"/>
	<img src="docs/src/assets/animations/quadrotor.gif" width="150"/>
</p>

## System Identification
<p float="left">
	<img src="docs/src/assets/animations/box_learning.gif" width="200"/>
	<img src="docs/src/assets/animations/cone_learning.gif" width="200"/>
	<img src="docs/src/assets/animations/box_toss.gif" width="300"/>
</p>


## Interfacing Other Packages
| [ReinforcementLearning.jl](https://github.com/JuliaReinforcementLearning/ReinforcementLearning.jl): DQN | [ControlSystems.jl](https://github.com/JuliaControl/ControlSystems.jl): LQR |
| - | -|
| <img src="docs/src/assets/animations/cartpole_rl.gif" width="250"/> | <img src="docs/src/assets/animations/cartpole_lqr.gif" width="250"/> |




## Installation

`Dojo` can be added via the Julia package manager (type `]`):
```julia
pkg> add Dojo
```
For convenience mechanisms and environments, add `DojoEnvironments` additionally:
```julia
pkg> add DojoEnvironments
```

## Citing
```
@article{howelllecleach2022,
	title={Dojo: A Differentiable Simulator for Robotics},
	author={Howell, Taylor and Le Cleac'h, Simon and Bruedigam, Jan and Kolter, Zico and Schwager, Mac and Manchester, Zachary},
	journal={arXiv preprint arXiv:2203.00806},
	url={https://arxiv.org/abs/2203.00806},
	year={2022}
}
```

## How To Contribute
Please submit a pull request or open an issue.
See the [docs](https://dojo-sim.github.io/Dojo.jl/dev/contributing.html) for contribution ideas.
