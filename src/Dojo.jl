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

    set_position,
    set_velocity!,
    set_input!,
    add_force!,
    getid,
    get_node,
    get_body,
    get_joint_constraint,
    get_contact_constraint,
    simulate!,
    initializeConstraints!,
    disassemble,
    minimal_coordinates,
    minimal_velocities,

    RotX,
    RotY,
    RotZ,
    RGBA,

    szeros,
    sones,
    srand,

    getmechanism,
    initialize!,
    getdim,
    control_dimension,

    getmechanism,
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
    getMinState

# Utilities
include(joinpath("util", "util.jl"))
include(joinpath("util", "custom_static.jl"))
include(joinpath("util", "quaternion.jl"))

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
include(joinpath("contacts", "constraint.jl"))
include(joinpath("mechanism", "mechanism.jl"))
include(joinpath("mechanism", "system.jl"))
include(joinpath("mechanism", "methods.jl"))

# Simulation
include(joinpath("simulation", "step.jl"))
include(joinpath("simulation", "storage.jl"))
include(joinpath("simulation", "simulate.jl"))

# Mechanics
include(joinpath("mechanics", "momentum.jl"))
include(joinpath("mechanics", "energy.jl"))

# Joints
include(joinpath("joints", "joint.jl"))
include(joinpath("joints", "translational", "constraint.jl"))
include(joinpath("joints", "translational", "constraint_limits.jl"))
include(joinpath("joints", "translational", "input.jl"))
include(joinpath("joints", "translational", "force.jl"))
include(joinpath("joints", "translational", "minimal.jl"))
include(joinpath("joints", "rotational", "constraint.jl"))
include(joinpath("joints", "rotational", "constraint_limits.jl"))
include(joinpath("joints", "rotational", "input.jl"))
include(joinpath("joints", "rotational", "torque.jl"))
include(joinpath("joints", "rotational", "minimal.jl"))
include(joinpath("joints", "prototypes.jl"))

# Contacts 
include(joinpath("contacts", "contact.jl"))
include(joinpath("contacts", "cone.jl"))
include(joinpath("contacts", "impact.jl"))
include(joinpath("contacts", "linear.jl"))
include(joinpath("contacts", "nonlinear.jl"))
include(joinpath("contacts", "constructor.jl"))
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
include(joinpath("ui", "initialize.jl"))
include(joinpath("ui", "urdf.jl"))
include(joinpath("ui", "convert_shape.jl"))
include(joinpath("ui", "visualize.jl"))
include(joinpath("ui", "visualizer_utils.jl"))
include(joinpath("ui", "colors.jl"))

# Gradients
include(joinpath("gradients", "index.jl"))
include(joinpath("gradients", "analytical.jl"))
include(joinpath("gradients", "finite_difference.jl"))

# Environments
include(joinpath("..", "env", "mechanisms.jl"))
include(joinpath("..", "env", "environment.jl"))

# Utilities
include(joinpath("..", "examples", "trajectory_optimization", "utils.jl"))
include(joinpath("..", "examples", "reinforcement_learning", "ars.jl"))

end
