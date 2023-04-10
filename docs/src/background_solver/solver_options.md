# Solver Options

The solver has several options accessible via [`SolverOptions`](@ref). Below is a list describing their effect on the solver's behavior, typical values, and if they need to be tuned by the user.


| option                  | default |range                                  | effect | tuning         | 
| ----------------------- | --------|-------------------------------------- | ------ | -------------- |
| `rtol`                  | ``10^{-4}`` |``[10^{-6}, 10^{-2}]``| larger leads to faster solve (usually takes the same value as `btol`)          | rarely |
| `btol`                  | ``10^{-4}`` |``[10^{-6}, 10^{-2}]``| larger results in smoothed contact dynamics and faster solve                   | rarely |
| `ls_scale`              | ``0.5``     |``[0.3, 0.8]``        | larger potentially increase step size at the cost of more residual evaluations | never |
| `max_ls`                | ``10``      |``[1, 25]``           | larger allows for taking smaller steps                                         | never |
| `undercut`              | ``+\infty`` |``[2, +\infty]``      | larger is more robust but can generate stiffer gradients                       | rarely |
| `no_progress_max`       | ``3``       |``[3, 5]``            | smaller will increase the undercut faster                                      | never |
| `no_progress_undercut`  | ``10``      |``[3, 100]``          | larger will increase the undercut faster                                       | never |
| `verbose`               | ``\text{false}``   |``\{\text{true}, \text{false}\}``   | printing the status of the solver                                              | often |
