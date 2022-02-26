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
            "define_mechanism.md",
            "load_mechanism.md",
            "mechanism_interfaces.md",
           ],

        "Creating a Simulation" => [
            "define_simulation.md",
            "define_controller.md",
           ],

        "Creating Environments" => [
            "define_environment.md",
            "load_environment.md",
        ],

        "Gradients from Simulator" => [
            "gradients.md",
           ],

        "Environments" => [
            "atlas_env.md",
            "quadruped_env.md",
            "rexhopper_env.md",
            "classic_env.md",
            "gym_env.md",
        ],

        "Examples" => [
            "simulation.md",
            "control.md",
            "trajectory_optimization.md",
            "reinforcement_learning.md",
            "system_identification.md",
           ],

        "State Representations" => [
            "maximal_representation.md",
            "minimal_representation.md",
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
        "contributing.md",
        "citing.md"
    ]
)

deploydocs(
    repo = "github.com/dojo-sim/Dojo.jl.git",
)
