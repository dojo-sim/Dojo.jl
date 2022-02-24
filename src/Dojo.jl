module Dojo

using LinearAlgebra
using StaticArrays
using ForwardDiff
using FiniteDiff
using StaticArrays: SUnitRange
using Rotations
using Rotations: RotationError, params, lmult, rmult, tmat, vmat, hmat, skew, pure_quaternion
using Colors: RGBA, RGB
using LightXML
using Parameters
using SparseArrays

using Plots
using Random
using MeshCat
using Meshing
using GeometryBasics
using LightGraphs
using DocStringExtensions
using JLD2

using Statistics
import Distributions: Uniform, Normal

export Origin,
    Body,
    JointConstraint,
    ContactConstraint,
    Mechanism,
    Storage,
    UnitQuaternion,
    Rotational,
    Translational,

    Box,
    Cylinder,
    Sphere,
    Pyramid,
    Mesh,

    Fixed,
    Prismatic,
    Planar,
    FixedOrientation,
    Revolute,
    Cylindrical,
    PlanarAxis,
    FreeRevolute,
    Orbital,
    PrismaticOrbital,
    PlanarOrbital,
    FreeOrbital,
    Spherical,
    CylindricalFree,
    PlanarFree,

    NonlinearContact,
    UnitQuaternion,

    set_maximal_configuration!,
    set_minimal_coordinates!,
    set_maximal_velocity!,
    set_input!,
    add_input!,
    get_node,
    get_body,
    get_joint_constraint,
    get_contact_constraint,
    simulate!,
    disassemble,

    minimal_coordinates,
    minimal_velocities,
    minimal_dimension,
    maximal_dimension,
    maximal_to_minimal,
    minimal_to_maximal,
    maximal_to_minimal_jacobian,
    minimal_to_maximal_jacobian,

    RotX,
    RotY,
    RotZ,
    RGBA,

    szeros,
    sones,
    srand,

    getdim,
    control_dimension,

    get_mechanism,
    initialize!,
    get_data,
    set_data!,
    get_solution,
    attitude_jacobian,
    finitediff_sol_matrix,
    full_matrix,
    full_data_matrix,
    finitediff_data_matrix,
    finitediff_sensitivity,

    momentum,
    kinetic_energy,
    potential_energy,
    mechanical_energy,

    mehrotra!,
    SolverOptions,

    Environment,
    get_minimal_state

# Utilities
include(joinpath("util", "util.jl"))
include(joinpath("util", "custom_static.jl"))
include(joinpath("util", "quaternion.jl"))
include(joinpath("util", "mrp.jl"))

# Graph system
include(joinpath("graph", "entry.jl"))
include(joinpath("graph", "system.jl"))
include(joinpath("graph", "setup_functions.jl"))
include(joinpath("graph", "ldu.jl"))

# Mechanism
include(joinpath("mechanism", "shapes.jl"))
include(joinpath("mechanism", "node.jl"))
include(joinpath("mechanism", "state.jl"))
include(joinpath("mechanism", "body.jl"))
include(joinpath("mechanism", "origin.jl"))
include(joinpath("joints", "constraint.jl"))
include(joinpath("contacts", "constructor.jl"))
include(joinpath("mechanism", "id.jl"))
include(joinpath("mechanism", "mechanism.jl"))
include(joinpath("mechanism", "maximal.jl"))
include(joinpath("mechanism", "minimal.jl"))
include(joinpath("mechanism", "system.jl"))
include(joinpath("mechanism", "methods.jl"))
include(joinpath("mechanism", "set.jl"))

# Simulation
include(joinpath("simulation", "step.jl"))
include(joinpath("simulation", "storage.jl"))
include(joinpath("simulation", "simulate.jl"))

# Mechanics
include(joinpath("mechanics", "momentum.jl"))
include(joinpath("mechanics", "energy.jl"))

# Joints
include(joinpath("joints", "joint.jl"))
include(joinpath("joints", "translational", "constructor.jl"))
include(joinpath("joints", "translational", "impulses.jl"))
include(joinpath("joints", "translational", "input.jl"))
include(joinpath("joints", "translational", "springs.jl"))
include(joinpath("joints", "translational", "dampers.jl"))
include(joinpath("joints", "translational", "minimal.jl"))
include(joinpath("joints", "rotational", "constructor.jl"))
include(joinpath("joints", "rotational", "impulses.jl"))
include(joinpath("joints", "rotational", "input.jl"))
include(joinpath("joints", "rotational", "springs.jl"))
include(joinpath("joints", "rotational", "dampers.jl"))
include(joinpath("joints", "rotational", "minimal.jl"))
include(joinpath("joints", "limits.jl"))
include(joinpath("joints", "prototypes.jl"))
include(joinpath("joints", "minimal.jl"))
include(joinpath("joints", "impulses.jl"))

# Contacts
include(joinpath("contacts", "constraints.jl"))
include(joinpath("contacts", "contact.jl"))
include(joinpath("contacts", "cone.jl"))
include(joinpath("contacts", "impact.jl"))
include(joinpath("contacts", "linear.jl"))
include(joinpath("contacts", "nonlinear.jl"))
include(joinpath("contacts", "utils.jl"))

# Solver
include(joinpath("solver", "methods.jl"))
include(joinpath("solver", "mehrotra.jl"))
include(joinpath("solver", "linesearch.jl"))

# Integrator
include(joinpath("integrators", "integrator.jl"))
include(joinpath("integrators", "constraint.jl"))

# User interface
include(joinpath("ui", "mechanism_ui.jl"))
include(joinpath("ui", "mechanism_traversal.jl"))
include(joinpath("ui", "urdf.jl"))
include(joinpath("ui", "convert_shape.jl"))
include(joinpath("ui", "visualize.jl"))
include(joinpath("ui", "visualizer_utils.jl"))
include(joinpath("ui", "colors.jl"))

# Gradients
include(joinpath("gradients", "finite_difference.jl"))
include(joinpath("gradients", "data.jl"))
include(joinpath("gradients", "data_gradients.jl"))
include(joinpath("gradients", "utils.jl"))

# Environments
include(joinpath("..", "env", "mechanisms.jl"))
include(joinpath("..", "env", "environment.jl"))

# Utilities
include(joinpath("..", "examples", "reinforcement_learning", "ars.jl"))

end
