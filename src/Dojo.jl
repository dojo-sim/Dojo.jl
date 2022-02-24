module Dojo

using LinearAlgebra
using Random
using StaticArrays
using SparseArrays
using FiniteDiff
using StaticArrays: SUnitRange
using Rotations
using Rotations: RotationError, params, lmult, rmult, tmat, vmat, hmat, skew, pure_quaternion
using Parameters

using Colors: RGBA, RGB
using FFMPEG
using LightXML
using MeshCat
using Meshing
using GeometryBasics
using LightGraphs

using JLD2
using DocStringExtensions

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

    set_maximal_coordinates!,
    set_minimal_coordinates!,
    set_maximal_velocities!,
    set_minimal_velocities!,
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
include(joinpath("utilities", "methods.jl"))
include(joinpath("utilities", "custom_static.jl"))

# Orientation
include(joinpath("orientation", "quaternion.jl"))
include(joinpath("orientation", "mrp.jl"))
include(joinpath("orientation", "axis_angle.jl"))
include(joinpath("orientation", "rotate.jl"))

# Graph system
include(joinpath("graph", "entry.jl"))
include(joinpath("graph", "system.jl"))
include(joinpath("graph", "setup_functions.jl"))
include(joinpath("graph", "ldu.jl"))

# Node 
include(joinpath("mechanism", "node.jl"))
include(joinpath("mechanism", "id.jl"))

# Bodies
include(joinpath("bodies", "shapes.jl"))
include(joinpath("bodies", "state.jl"))
include(joinpath("bodies", "constructor.jl"))
include(joinpath("bodies", "origin.jl"))
include(joinpath("bodies", "set.jl"))

# Mechanism
include(joinpath("joints", "constraint.jl"))
include(joinpath("contacts", "constructor.jl"))
include(joinpath("contacts", "contact.jl"))

include(joinpath("mechanism", "constructor.jl"))
include(joinpath("mechanism", "gravity.jl"))
include(joinpath("mechanism", "state.jl"))
include(joinpath("mechanism", "system.jl"))
include(joinpath("mechanism", "methods.jl"))
include(joinpath("mechanism", "set.jl"))
include(joinpath("mechanism", "get.jl"))
include(joinpath("mechanism", "urdf.jl"))
include(joinpath("mechanism", "traversal.jl"))

# Simulation
include(joinpath("simulation", "step.jl"))
include(joinpath("simulation", "storage.jl"))
include(joinpath("simulation", "simulate.jl"))

# Mechanics
include(joinpath("mechanics", "momentum.jl"))
include(joinpath("mechanics", "energy.jl"))

# Joints
include(joinpath("joints", "orthogonal.jl"))
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
include(joinpath("contacts", "cone.jl"))
include(joinpath("contacts", "impact.jl"))
include(joinpath("contacts", "linear.jl"))
include(joinpath("contacts", "nonlinear.jl"))
include(joinpath("contacts", "utils.jl"))

# Solver
include(joinpath("solver", "methods.jl"))
include(joinpath("solver", "options.jl"))
include(joinpath("solver", "mehrotra.jl"))
include(joinpath("solver", "line_search.jl"))

# Integrator
include(joinpath("integrators", "integrator.jl"))
include(joinpath("integrators", "constraint.jl"))

# Visualizer
include(joinpath("visuals", "visualizer.jl"))
include(joinpath("visuals", "utilities.jl"))
include(joinpath("visuals", "colors.jl"))

# Data 
include(joinpath("mechanism", "data.jl"))

# Gradients
include(joinpath("gradients", "finite_difference.jl"))
include(joinpath("gradients", "state.jl"))
include(joinpath("gradients", "data.jl"))
include(joinpath("gradients", "utilities.jl"))

# Environments
include(joinpath("..", "env", "mechanisms.jl"))
include(joinpath("..", "env", "environment.jl"))

# Optimizers
include(joinpath("..", "examples", "reinforcement_learning", "ars.jl"))

end
