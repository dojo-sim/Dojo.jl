# Solver Options

The solver we implemented has several options all accessible via [`SolverOptions`](@ref). Here is a list describing what they are, their typical values and their impact on the solver's behavior.


Increasing $b_{tol}$ results in smoothed contact dynamics.

| option                  | min $\cdot$ default $\cdot$ max values | effect | require tuning |
| ----------------------- | -------------------------------------- | ------ | -------------- |
| `rtol`                  | $10^{-6}$ $\cdot$ $10^{-4}$ $\cdot$ $10^{-2}$ | larger leads to faster solve (usually takes the same value as `btol`) | rarely |
| `btol`                  | $10^{-6}$ $\cdot$ $10^{-4}$ $\cdot$ $10^{-2}$ | larger results in smoothed contact dynamics and faster solve | rarely |
| `ls_scale`              | $0.3$ $\cdot$ $0.5$ $\cdot$ $0.8$             | larger potentially increase step size at the cost of more residual evaluations | never |
| `max_ls`                | $5$ $\cdot$ $10$ $\cdot$ $15$                 | larger allows for taking smaller steps | never |
| `undercut`              | $2$ $\cdot$ $+\infty$ $\cdot$ $+\infty$       | larger is more robust but can generate stiffer gradients | rarely |
| `no_progress_max`       | $3$ $\cdot$ $3$ $\cdot$ $5$                   | smaller will increase the undercut faster | never |
| `no_progress_undercut`  | $3$ $\cdot$ $10$ $\cdot$ $100$                | larger will increase the undercut faster | never |
| `verbose`               | true $\cdot$ false $\cdot$ false              | printing the status of the solver | often |
