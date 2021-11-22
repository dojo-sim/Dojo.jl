using LinearAlgebra
using Parameters
using StaticArrays
using ForwardDiff
using FiniteDiff
using StaticArrays: SUnitRange
using Rotations
using Rotations: RotationError, pure_quaternion, params, lmult, rmult, tmat, vmat, hmat, skew
using Colors: RGBA, RGB
using LightXML
using SparseArrays
using Parameters
using Test

using Plots
using Random
using MeshCat
using GeometryBasics
using LightGraphs
using DocStringExtensions
using FiniteDiff

export Origin,
    Body,
    EqualityConstraint,
    InequalityConstraint,
    LinearInequalityConstraint,
    Mechanism,
    Controller,
    Storage,
    UnitQuaternion,

    Box,
    Cylinder,
    Sphere,
    Pyramid,
    Mesh,

    Floating,
    Prismatic,
    Spherical,
    Cylindrical,
    Revolute,
    Planar,
    PlanarFree,
    Fixed,
    FixedOrientation,
    CylindricalFree,

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
    srand

include(joinpath(module_dir(), "src", "util", "util.jl"))
include(joinpath(module_dir(), "src", "util", "custom_static.jl"))
include(joinpath(module_dir(), "src", "util", "customdict.jl"))
include(joinpath(module_dir(), "src", "util", "quaternion.jl"))

include(joinpath(module_dir(), "src", "graph", "entry.jl"))
include(joinpath(module_dir(), "src", "graph", "system.jl"))
include(joinpath(module_dir(), "src", "graph", "setup_functions.jl"))
include(joinpath(module_dir(), "src", "graph", "ldu.jl"))


include(joinpath(module_dir(), "src", "optional_components", "shapes.jl"))
include(joinpath(module_dir(), "src", "optional_components", "storage.jl"))

include(joinpath(module_dir(), "src", "main_components", "component.jl"))
include(joinpath(module_dir(), "src", "main_components", "state.jl"))
include(joinpath(module_dir(), "src", "main_components", "body.jl"))
include(joinpath(module_dir(), "src", "main_components", "abstractconstraint.jl"))
include(joinpath(module_dir(), "src", "main_components", "equalityconstraint.jl"))
include(joinpath(module_dir(), "src", "main_components", "inequalityconstraint.jl"))
include(joinpath(module_dir(), "src", "main_components", "controller.jl"))
include(joinpath(module_dir(), "src", "main_components", "mechanism_struct.jl"))
include(joinpath(module_dir(), "src", "main_components", "system.jl"))
include(joinpath(module_dir(), "src", "main_components", "mechanism_functions.jl"))

include(joinpath(module_dir(), "src", "joints", "abstract_joint.jl"))

include(joinpath(module_dir(), "src", "bounds", "bound.jl"))
include(joinpath(module_dir(), "src", "bounds", "contact_bounds.jl"))
include(joinpath(module_dir(), "src", "bounds", "impact_bounds.jl"))
include(joinpath(module_dir(), "src", "bounds", "linear_contact_bounds.jl"))

include(joinpath(module_dir(), "src", "joints", "joint.jl"))
include(joinpath(module_dir(), "src", "joints", "translational.jl"))
include(joinpath(module_dir(), "src", "joints", "rotational.jl"))
include(joinpath(module_dir(), "src", "joints", "genericjoint.jl"))
include(joinpath(module_dir(), "src", "joints", "prototypes.jl"))

include(joinpath(module_dir(), "src", "solver", "solverfunctions.jl"))
include(joinpath(module_dir(), "src", "solver", "initconstraints.jl"))
include(joinpath(module_dir(), "src", "solver", "mehrotra.jl"))
include(joinpath(module_dir(), "src", "solver", "linesearch.jl"))

include(joinpath(module_dir(), "src", "discretization", "Linear.jl"))

include(joinpath(module_dir(), "src", "ui", "mechanism_ui.jl"))
include(joinpath(module_dir(), "src", "ui", "simulate.jl"))
include(joinpath(module_dir(), "src", "ui", "initialize.jl"))
include(joinpath(module_dir(), "src", "ui", "urdf.jl"))

include(joinpath(module_dir(), "examples", "diff_tools.jl"))
include(joinpath(module_dir(), "examples", "finitediff_tools.jl"))
include(joinpath(module_dir(), "examples", "mechanism_zoo.jl"))
include(joinpath(module_dir(), "src", "joints", "force.jl"))
include(joinpath(module_dir(), "src", "joints", "torque.jl"))
include(joinpath(module_dir(), "src", "vis", "convertshape.jl"))
include(joinpath(module_dir(), "src", "vis", "visualize.jl"))
include(joinpath(module_dir(), "src", "optional_components", "energy.jl"))
