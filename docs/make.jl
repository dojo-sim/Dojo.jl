push!(LOAD_PATH, "../src/")

using Documenter, Dojo

makedocs(
    modules = [Dojo],
    format = Documenter.HTML(prettyurls = false),
    sitename = "Dojo",
    pages = [
        ##############################################
        ## MAKE SURE TO SYNC WITH docs/src/index.md ##
        ##############################################
        "Basics" => [
            "index.md",
            "install.md",
            "get_started.md",
            "notations.md",
           ],

        "Creating a Mechanism" => [
            "def_mechanism.md",
            "load_mechanism.md",
            "mechanism_interfaces.md",
           ],

        "Creating a Simulation" => [
            "def_simulation.md",
            "def_controller.md",
           ],

        "Creating Environments" => [
            "def_environment.md",
            "load_environment.md",
        ],

        "Gradients from Simulator" => [
            "gradients.md",
           ],

        "Examples" => [
            "simulation.md",
            "control.md",
            "trajectory_optimization.md",
            "reinforcement_learning.md",
            "system_id.md",
           ],

        "State Representations" => [
            "maximal_rep.md",
            "minimal_rep.md",
        ],

        "Contact Models" => [
            "impact.md",
            "friction.md",
        ],

        "Interior-Point Solver" => [
            "interior_point.md",
            "solver_options.md",
        ],

        "faq.md",
        "api.md",
        "citing.md"
    ]
)

deploydocs(
    repo = "github.com/dojo-sim/Dojo.jl.git",
)
