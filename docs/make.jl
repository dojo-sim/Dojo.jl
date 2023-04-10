push!(LOAD_PATH, "../src/")

using Documenter, Dojo

# copy animations to docs/src/assets/animations
# path_assets_animations = joinpath(@__DIR__, "src/assets/animations")
# !isdir(path_assets_animations) && mkdir(path_assets_animations)
# path_animations = joinpath(@__DIR__, "../examples/animations")
# files = readdir(path_animations)
# filter!(x -> endswith(x, ".gif"), files)
# for file in files 
#     cp(joinpath(path_animations, file), joinpath(path_assets_animations, file), force=true)
# end

makedocs(
    modules = [Dojo],
    format = Documenter.HTML(prettyurls=false),
    sitename = "Dojo",
    pages = [
        ##############################################
        ## MAKE SURE TO SYNC WITH docs/src/index.md ##
        ##############################################
        "index.md",

        "Examples" => [
            "examples/simulation.md",
            "examples/reinforcement_learning.md",
            "examples/system_identification.md",
            "examples/trajectory_optimization.md",
           ],

        "Creating a Mechanism" => [
            "creating_mechanism/overview.md",
            "creating_mechanism/mechanism_directly.md",
            "creating_mechanism/mechanism_existing.md",
            "creating_mechanism/environment_existing.md",
            "creating_mechanism/tippetop.md",
            "creating_mechanism/quadruped.md",
           ],

        "Creating a Simulation" => [
            "creating_simulation/define_simulation.md",
            "creating_simulation/define_controller.md",
            "creating_simulation/simulation_with_gradients.md",
           ],

        "Environments" => [
            "environments/overview.md",
        ],

        "Background: Contact Models" => [
            "background_contact/contact_models.md",
            "background_contact/impact.md",
            "background_contact/nonlinear_friction.md",
            "background_contact/linearized_friction.md",
            "background_contact/collisions.md",
        ],

        "Background: Representations" => [
            "background_representations/maximal_representation.md",
            "background_representations/minimal_representation.md",
            "background_representations/gradients.md",
        ],

        "Background: Solver" => [
            "background_solver/interior_point.md",
            "background_solver/solver_options.md",
        ],

        "api.md",
        "contributing.md",
        "citing.md"
    ]
)

deploydocs(
    repo = "github.com/dojo-sim/Dojo.jl.git",
)
