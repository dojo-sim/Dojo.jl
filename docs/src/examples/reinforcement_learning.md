# Reinforcement Learning

We have implemented a few learning examples.

## Ant

```@raw html
<img src="../assets/animations/ant_ars_no_grid.gif" width="300"/>
```
Policy optimization is performed using the reinforcement-learning algorithm [augmented random search (ARS)](https://arxiv.org/abs/1803.07055) to optimize static linear policies for locomotion.
The insect-like robot has rewards on forward velocity and survival and costs on control usage and contact forces.

## Quadruped
```@raw html
<img src="../assets/animations/quadruped_walking.gif" width="300"/>
```
A very basic random-sampling algorithm is used to find parameters for the periodic gait of a quadruped. 

## Cartpole
```@raw html
<img src="../assets/animations/cartpole_rl.gif" width="300"/>
```
We have modified the cartpole example in the `ReinforcementLearning` package to use `Dojo`'s dynamics. This allows us to combine advanced learning algorithms with accurate dynamics simulation.