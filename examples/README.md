# Dojo.jl examples

This directory contains examples using Dojo.
The `.jl` files in each subdirectory can be executed directly or they can processed using [Literate.jl](https://github.com/fredrikekre/Literate.jl) to obtain notebooks. 

Building of the `Dojo` package generates notebooks and they can be run locally by performing the following steps:

1. [install Dojo.jl](https://github.com/dojo-sim/Dojo.jl)
2. [install IJulia](https://github.com/JuliaLang/IJulia.jl) (`add` it to the default project)
3. in the Julia REPL, run (do once)
   ```
   using Pkg
   Pkg.build("Dojo")
   ```
4. interact with notebooks
   ```
   using IJulia, Dojo
   notebook(dir=joinpath(dirname(pathof(Dojo)), "..", "examples"))
   ````