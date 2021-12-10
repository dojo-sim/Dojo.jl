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
using GeometryBasics
using LightGraphs
using DocStringExtensions


export Origin,
    Body,
    EqualityConstraint,
    InequalityConstraint,
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
    finitediff_sensitivity,

    momentum, 
    kineticEnergy,
    potentialEnergy,
    mechanicalEnergy,

    mehrotra!,
    InteriorPointOptions

# Utilities
include(joinpath("util", "util.jl"))
include(joinpath("util", "custom_static.jl"))
include(joinpath("util", "customdict.jl"))
include(joinpath("util", "quaternion.jl"))

# Graph system
include(joinpath("graph", "entry.jl"))
include(joinpath("graph", "system.jl"))
include(joinpath("graph", "setup_functions.jl"))
include(joinpath("graph", "ldu.jl"))

# Mechanism
include(joinpath("mechanism", "shapes.jl"))
include(joinpath("mechanism", "component.jl"))
include(joinpath("mechanism", "state.jl"))
include(joinpath("mechanism", "body.jl"))
include(joinpath("mechanism", "abstract_constraint.jl"))
include(joinpath("mechanism", "equality_constraint.jl"))
include(joinpath("mechanism", "inequality_constraint.jl"))
include(joinpath("mechanism", "mechanism.jl"))
include(joinpath("mechanism", "system.jl"))
include(joinpath("mechanism", "methods.jl"))

# Simulation
include(joinpath("simulation", "step.jl"))
include(joinpath("simulation", "storage.jl"))
include(joinpath("simulation", "simulate.jl"))

# Energy
include(joinpath("mechanism", "energy.jl"))

# Joints
include(joinpath("joints", "abstract_joint.jl"))
include(joinpath("joints", "joint.jl"))
include(joinpath("joints", "translational.jl"))
include(joinpath("joints", "rotational.jl"))
include(joinpath("joints", "generic_joint.jl"))
include(joinpath("joints", "prototypes.jl"))
include(joinpath("joints", "force.jl"))
include(joinpath("joints", "torque.jl"))

# Inequality constraints
include(joinpath("bounds", "bound.jl"))
include(joinpath("bounds", "cone.jl"))
include(joinpath("bounds", "contact.jl"))
include(joinpath("bounds", "impact.jl"))
include(joinpath("bounds", "linear_contact.jl"))

# Solver
include(joinpath("solver", "methods.jl"))
include(joinpath("solver", "mehrotra.jl"))
include(joinpath("solver", "linesearch.jl"))

# Variational integrator
include(joinpath("discretization", "integrator.jl"))
include(joinpath("discretization", "body.jl"))

# User interface
include(joinpath("ui", "mechanism_ui.jl"))
include(joinpath("ui", "initialize.jl"))
include(joinpath("ui", "urdf.jl"))
include(joinpath("ui", "convert_shape.jl"))
include(joinpath("ui", "visualize.jl"))

# Differentiation
include(joinpath("diff", "diff_tools.jl"))
include(joinpath("diff", "finitediff_tools.jl"))

# Environments
include(joinpath("..", "env", "mechanisms.jl"))
include(joinpath("..", "env", "environment.jl"))

# Utilities
include(joinpath("..", "examples", "trajectory_optimization", "utils.jl"))

end
