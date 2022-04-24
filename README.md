[![CI](https://github.com/dojo-sim/Dojo.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/dojo-sim/Dojo.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/dojo-sim/Dojo.jl/branch/main/graph/badge.svg?token=NMS3JQZ2OE)](https://codecov.io/gh/dojo-sim/Dojo.jl)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://dojo-sim.github.io/Dojo.jl/dev)

# Dojo
A differentiable simulator for robotics
- arXiv preprint: https://arxiv.org/abs/2203.00806
- Python interface: https://github.com/dojo-sim/dojopy
- site: https://sites.google.com/view/dojo-sim
- presentation video: https://youtu.be/TRtOESXJxJQ

# Examples

## Simulation

### Atlas drop
<img src="examples/animations/atlas_drop.gif" width="100"/>

### REx Hopper drop
<img src="examples/animations/rexhopper.gif" width="200"/>

### Astronaut
<img src="examples/animations/astronaut.gif" width="200"/>

### Friction-cone comparison
- blue = Dojo nonlinear friction cone
- orange = Dojo linear friction cone
- black = MuJoCo nonlinear friction cone
- magenta = MuJoCo linear friction model
<img src="examples/animations/cone_compare_mujoco.gif" height="150"/>

### Dzhanibekov effect
<img src="examples/animations/dzhanibekov.gif" width="150"/>

### Tippe top
<img src="examples/animations/tippetop.gif" width="150"/>

### Pendulum swing-up
<img src="examples/animations/pendulum.gif" width="150"/>

## Trajectory Optimization

### Box
<img src="examples/animations/box_right.gif" width="200"/>

<img src="examples/animations/box_up.gif" width="95"/>

### Hopper
<img src="examples/animations/hopper_max.gif" width="100"/>

### Quadruped
<img src="examples/animations/quadruped_min.gif" width="200"/>

### Atlas
<img src="examples/animations/atlas_ilqr.gif" width="200"/>

### Cart-pole
<img src="examples/animations/cartpole_max.gif" width="200"/>

## Reinforcement Learning

### Half Cheetah
<img src="examples/animations/halfcheetah_ars.gif" width="600"/>

### Ant
<img src="examples/animations/ant_ars_no_grid.gif" width="300"/>

## Real-To-Sim

### Learning
Learning geometry:
<img src="examples/animations/box_learning.gif" width="200"/>
Learning friction coefficient:
<img src="examples/animations/cone_learning.gif" width="200"/>

### Toss
<img src="examples/animations/box_toss.gif" width="300"/>

## Installation

`Dojo` can be added via the Julia package manager (type `]`):
```julia
pkg> add Dojo
```

The latest version can be added by calling:
```julia
pkg> add Dojo#main
```

## Citing
```
@article{howelllecleach2022,
	title={Dojo: A Differentiable Simulator for Robotics},
	author={A. Howell, Taylor and Le Cleac'h, Simon and Kolter, Zico and Schwager, Mac and Manchester, Zachary},
	journal={arXiv preprint arXiv:2203.00806},
	url={https://arxiv.org/abs/2203.00806},
	year={2022}
}
```

## How To Contribute
Please submit a pull request, open an issue, or reach out to: thowell@stanford.edu (Taylor) or simonlc@stanford.edu (Simon)
