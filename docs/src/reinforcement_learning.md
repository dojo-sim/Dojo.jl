# Reinforcement Learning

Policy optimization is performed using the reinforcement-learning algorithm [augmented random search (ARS)](https://arxiv.org/abs/1803.07055) to optimize static linear policies for locomotion. A number of [Gym-like environments](https://gym.openai.com/) are created with Dojo for this application.

## [Half Cheetah](https://github.com/dojo-sim/Dojo.jl/blob/main/examples/reinforcement_learning/halfcheetah_ars.jl)

```@raw html
<img src="./assets/animations/halfcheetah_ars.gif" width="600"/>
```

The cheetah-like robot has rewards on forward velocity and costs on control usage. 

## [Ant](https://github.com/dojo-sim/Dojo.jl/blob/main/examples/reinforcement_learning/ant_ars.jl) 

```@raw html
<img src="./assets/animations/ant_ars_no_grid.gif" width="300"/>
```

The insect-like robot has rewards on forward velocity and survival and costs on control usage and contact forces.

## Gradient-Free
The original ARS is a derivative-free algorithm and only utilizes the relative costs for each sample policy in order to assemble a new search direction for the policy parameters. While this approach is general, it is potentially very sample inefficient.

## Gradient-Based
We formulate a gradient-based version, augmented gradient search (AGS), that utilizes gradients from Dojo. We find that by using this information the policy optimization can often be an order of magntitude more sample efficient in terms of queries to the simulator.