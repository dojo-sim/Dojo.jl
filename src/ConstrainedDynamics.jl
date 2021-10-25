module ConstrainedDynamics

using LinearAlgebra
using StaticArrays
using ForwardDiff
using StaticArrays: SUnitRange
using Rotations
using Rotations: RotationError, pure_quaternion, params, lmult, rmult, tmat, vmat, hmat, skew
using Colors: RGBA, RGB
using LightXML
# using GraphBasedSystems
# using GraphBasedSystems: Entry
using Parameters
using SparseArrays
# using Symbolics
# using FiniteDiff

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
    srand,

    getmechanism,
    initialize!,
    getdim,
    getcontroldim

include(joinpath("util", "util.jl"))
include(joinpath("util", "custom_static.jl"))
include(joinpath("util", "customdict.jl"))
include(joinpath("util", "quaternion.jl"))

include(joinpath(module_dir(), "src", "graph", "entry.jl"))
include(joinpath(module_dir(), "src", "graph", "system.jl"))
include(joinpath(module_dir(), "src", "graph", "setup_functions.jl"))
include(joinpath(module_dir(), "src", "graph", "ldu.jl"))

include(joinpath("optional_components", "shapes.jl"))
include(joinpath("optional_components", "storage.jl"))

include(joinpath("main_components", "component.jl"))
include(joinpath("main_components", "state.jl"))
include(joinpath("main_components", "body.jl"))
include(joinpath("main_components", "abstractconstraint.jl"))
include(joinpath("main_components", "equalityconstraint.jl"))
include(joinpath("main_components", "inequalityconstraint.jl"))
include(joinpath("main_components", "friction.jl"))
include(joinpath("main_components", "controller.jl"))
include(joinpath("main_components", "mechanism_struct.jl"))
include(joinpath("main_components", "system.jl"))
include(joinpath("main_components", "mechanism_functions.jl"))

include(joinpath("joints", "abstract_joint.jl"))

include(joinpath("bounds", "bound.jl"))
include(joinpath("bounds", "impact.jl"))
include(joinpath("bounds", "cone_bounds.jl"))
include(joinpath("bounds", "contact_bounds.jl"))
include(joinpath("bounds", "friction_bounds.jl"))

include(joinpath("joints", "joint.jl"))
include(joinpath("joints", "translational.jl"))
include(joinpath("joints", "rotational.jl"))
include(joinpath("joints", "genericjoint.jl"))
include(joinpath("joints", "prototypes.jl"))
# include(joinpath("joints", "friction.jl"))

include(joinpath("solver", "solverfunctions.jl"))
include(joinpath("solver", "initconstraints.jl"))
include(joinpath("solver", "newton.jl"))
include(joinpath("solver", "mehrotra.jl"))
include(joinpath("solver", "linesearch.jl"))
include(joinpath("optional_components", "linearization.jl"))

include(joinpath("discretization", "Linear.jl"))
# include(joinpath("discretization", "Quadratic.jl"))

include(joinpath("ui", "mechanism_ui.jl"))
include(joinpath("ui", "simulate.jl"))
include(joinpath("ui", "initialize.jl"))
include(joinpath("ui", "urdf.jl"))

include(joinpath("..", "examples", "dev", "mechanism_zoo.jl"))
include(joinpath("..", "examples", "dev", "diff_tools.jl"))
include(joinpath("joints", "fjoint.jl"))
include(joinpath("joints", "force.jl"))
include(joinpath("joints", "torque.jl"))
include(joinpath("vis", "convertshape.jl"))
include(joinpath("vis", "visualize.jl"))

include(joinpath("graph", "entry.jl"))
include(joinpath("graph", "system.jl"))
include(joinpath("graph", "setup_functions.jl"))
include(joinpath("graph", "ldu.jl"))

end
