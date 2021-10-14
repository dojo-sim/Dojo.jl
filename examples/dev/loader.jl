using LinearAlgebra
using Parameters
using StaticArrays
using ForwardDiff
using StaticArrays: SUnitRange
using Rotations
using Rotations: RotationError, pure_quaternion, params, lmult, rmult, tmat, vmat, hmat, skew
using Colors: RGBA, RGB
using LightXML
using GraphBasedSystems
using GraphBasedSystems: Entry
using SparseArrays
using Symbolics
using FiniteDiff

using Plots
using Random
using MeshCat
using GeometryBasics

using DocStringExtensions

export Origin,
    Body,
    EqualityConstraint,
    InequalityConstraint,
    LinearInequalityConstraint,
    FullInequalityConstraint,
    Friction,
    Friction,
    FullFriction,
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

    Impact,
    ConeBound,
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
    getfriction,
    getineqconstraint,
    simulate!,
    initializeConstraints!,
    disassemble,
    minimalCoordinates,
    minimalVelocities,
    linearsystem,

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

include(joinpath(module_dir(), "src", "optional_components", "shapes.jl"))
include(joinpath(module_dir(), "src", "optional_components", "storage.jl"))

include(joinpath(module_dir(), "src", "main_components", "component.jl"))
include(joinpath(module_dir(), "src", "main_components", "state.jl"))
include(joinpath(module_dir(), "src", "main_components", "body.jl"))
include(joinpath(module_dir(), "src", "main_components", "abstractconstraint.jl"))
include(joinpath(module_dir(), "src", "main_components", "equalityconstraint.jl"))
include(joinpath(module_dir(), "src", "main_components", "inequalityconstraint.jl"))
include(joinpath(module_dir(), "src", "main_components", "friction.jl"))
include(joinpath(module_dir(), "src", "main_components", "controller.jl"))
include(joinpath(module_dir(), "src", "main_components", "mechanism_struct.jl"))
include(joinpath(module_dir(), "src", "main_components", "system.jl"))
include(joinpath(module_dir(), "src", "main_components", "mechanism_functions.jl"))

include(joinpath(module_dir(), "src", "joints", "abstract_joint.jl"))

include(joinpath(module_dir(), "src", "bounds", "bound.jl"))
include(joinpath(module_dir(), "src", "bounds", "impact.jl"))
include(joinpath(module_dir(), "src", "bounds", "cone_bounds.jl"))
include(joinpath(module_dir(), "src", "bounds", "contact_bounds.jl"))
include(joinpath(module_dir(), "src", "bounds", "friction_bounds.jl"))

include(joinpath(module_dir(), "src", "joints", "joint.jl"))
include(joinpath(module_dir(), "src", "joints", "translational.jl"))
include(joinpath(module_dir(), "src", "joints", "rotational.jl"))
include(joinpath(module_dir(), "src", "joints", "genericjoint.jl"))
include(joinpath(module_dir(), "src", "joints", "prototypes.jl"))
# include(joinpath(module_dir(), "src", "joints", "friction.jl"))

include(joinpath(module_dir(), "src", "solver", "solverfunctions.jl"))
include(joinpath(module_dir(), "src", "solver", "initconstraints.jl"))
include(joinpath(module_dir(), "src", "solver", "newton.jl"))
include(joinpath(module_dir(), "src", "solver", "mehrotra.jl"))
include(joinpath(module_dir(), "src", "solver", "linesearch.jl"))
include(joinpath(module_dir(), "src", "optional_components", "linearization.jl"))

include(joinpath(module_dir(), "src", "discretization", "Linear.jl"))
# include(joinpath(module_dir(), "src", "discretization", "Quadratic.jl"))

include(joinpath(module_dir(), "src", "ui", "mechanism_ui.jl"))
include(joinpath(module_dir(), "src", "ui", "simulate.jl"))
include(joinpath(module_dir(), "src", "ui", "initialize.jl"))
include(joinpath(module_dir(), "src", "ui", "urdf.jl"))

include(joinpath(module_dir(), "examples", "dev", "diff_tools_control_contact.jl"))
include(joinpath(module_dir(), "examples", "dev", "mechanism_zoo.jl"))
include(joinpath(module_dir(), "src", "vis", "convertshape.jl"))
include(joinpath(module_dir(), "src", "vis", "visualize.jl"))
