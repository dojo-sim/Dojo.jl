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
using GeometryBasics
using LightGraphs
using DocStringExtensions


export Origin,
    Body,
    EqualityConstraint,
    InequalityConstraint,
    Mechanism,
    Controller,
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

    ContactBound,
    UnitQuaternion,

    setPosition!,
    setVelocity!,
    setForce!,
    addForce!,
    getid,
    getcomponent,
    getbody,
    geteqconstraint,
    getineqconstraint,
    simulate!,
    initializeConstraints!,
    disassemble,
    minimalCoordinates,
    minimalVelocities,

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
    controldim,

    getmechanism,
    initialize!,
    getdata,
    setdata!,
    getsolution,
    attitudejacobian,
    finitediff_sol_matrix,
    full_matrix,
    full_data_matrix,
    finitediff_data_matrix,
    finitediff_sensitivity

# Utilities
include(joinpath(module_dir(), "src", "util", "util.jl"))
include(joinpath(module_dir(), "src", "util", "custom_static.jl"))
include(joinpath(module_dir(), "src", "util", "customdict.jl"))
include(joinpath(module_dir(), "src", "util", "quaternion.jl"))

# Graph system
include(joinpath(module_dir(), "src", "graph", "entry.jl"))
include(joinpath(module_dir(), "src", "graph", "system.jl"))
include(joinpath(module_dir(), "src", "graph", "setup_functions.jl"))
include(joinpath(module_dir(), "src", "graph", "ldu.jl"))

# Simulation
include(joinpath(module_dir(), "src", "simulation", "controller.jl"))
include(joinpath(module_dir(), "src", "simulation", "storage.jl"))

# Mechanism
include(joinpath(module_dir(), "src", "mechanism", "shapes.jl"))
include(joinpath(module_dir(), "src", "mechanism", "component.jl"))
include(joinpath(module_dir(), "src", "mechanism", "state.jl"))
include(joinpath(module_dir(), "src", "mechanism", "body.jl"))
include(joinpath(module_dir(), "src", "mechanism", "abstract_constraint.jl"))
include(joinpath(module_dir(), "src", "mechanism", "equality_constraint.jl"))
include(joinpath(module_dir(), "src", "mechanism", "inequality_constraint.jl"))
include(joinpath(module_dir(), "src", "mechanism", "mechanism.jl"))
include(joinpath(module_dir(), "src", "mechanism", "system.jl"))
include(joinpath(module_dir(), "src", "mechanism", "methods.jl"))
include(joinpath(module_dir(), "src", "mechanism", "energy.jl"))

# Joints
include(joinpath(module_dir(), "src", "joints", "abstract_joint.jl"))
include(joinpath(module_dir(), "src", "joints", "joint.jl"))
include(joinpath(module_dir(), "src", "joints", "translational.jl"))
include(joinpath(module_dir(), "src", "joints", "rotational.jl"))
include(joinpath(module_dir(), "src", "joints", "generic_joint.jl"))
include(joinpath(module_dir(), "src", "joints", "prototypes.jl"))
include(joinpath(module_dir(), "src", "joints", "force.jl"))
include(joinpath(module_dir(), "src", "joints", "torque.jl"))

# Inequality constraints
include(joinpath(module_dir(), "src", "bounds", "bound.jl"))
include(joinpath(module_dir(), "src", "bounds", "cone.jl"))
include(joinpath(module_dir(), "src", "bounds", "contact.jl"))
include(joinpath(module_dir(), "src", "bounds", "impact.jl"))
include(joinpath(module_dir(), "src", "bounds", "linear_contact.jl"))

# Solver
include(joinpath(module_dir(), "src", "solver", "methods.jl"))
include(joinpath(module_dir(), "src", "solver", "mehrotra.jl"))
include(joinpath(module_dir(), "src", "solver", "linesearch.jl"))

# Variational integrator
include(joinpath(module_dir(), "src", "discretization", "integrator.jl"))
include(joinpath(module_dir(), "src", "discretization", "body.jl"))

# User interface
include(joinpath(module_dir(), "src", "ui", "mechanism_ui.jl"))
include(joinpath(module_dir(), "src", "ui", "simulate.jl"))
include(joinpath(module_dir(), "src", "ui", "initialize.jl"))
include(joinpath(module_dir(), "src", "ui", "urdf.jl"))
include(joinpath(module_dir(), "src", "ui", "convert_shape.jl"))
include(joinpath(module_dir(), "src", "ui", "visualize.jl"))

# Differentiation
include(joinpath(module_dir(), "src", "diff", "diff_tools.jl"))
include(joinpath(module_dir(), "src", "diff", "finitediff_tools.jl"))

# Environments
include(joinpath(module_dir(), "env", "mechanisms.jl"))
include(joinpath(module_dir(), "src", "environment", "environment.jl")
include(joinpath(module_dir(), "src", "environment", "pendulum.jl")
