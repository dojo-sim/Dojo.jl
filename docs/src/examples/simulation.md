# Simulation 

Dojo can simulate a number of interesting physical behaviors. 
We include notebooks (generated upon installation) for the examples below.

## Atlas Drop

```@raw html
<img src="../assets/animations/atlas_drop.gif" width="100"/>
```
The humanoid [Atlas](https://www.bostondynamics.com/atlas) is dropped onto a flat surface. 
Dojo is able to simulate hard contact and prevent interpenetration of the robot's feet with the floor.
In comparison, when the same system is simulated in [MuJoCo](https://mujoco.org), **centimeters** of interpenetration occur.


## Friction Cone Comparison

```@raw html
<img src="../assets/animations/cone_compare_mujoco.gif" width="300"/>
```
Blocks are simulated with initial velocity before impacting and sliding along a flat surface. We compare Dojo's nonlinear cone (blue) with a linearized approximation (orange) and MuJoCo's default linear cone (magenta). The linearized cones exhibit drift due to the approximation, whereas **Dojo's nonlinear cone produces the expected sliding behavior**.

## Dzhanibekov Effect

```@raw html
<img src="../assets/animations/dzhanibekov.gif" width="150"/>
```
Dojo simulates the [unstable rotational motion](https://en.wikipedia.org/wiki/Tennis_racket_theorem) of a rigid body about its second primary moment of inertia. Using [non-Euclidean optimization for quaternions](https://roboticexplorationlab.org/papers/planning_with_attitude.pdf) enables continuous simulation of rotating objects without singularity issues.

## Tippetop

```@raw html
<img src="../assets/animations/tippetop.gif" width="300"/>
```

A spinning object osciallates between up and down configurations as a result of its mass distribution.



