module Dojo

using LinearAlgebra
using StaticArrays
using ForwardDiff
using FiniteDiff
using StaticArrays: SUnitRange
using Rotations
using Rotations: RotationError, pure_quaternion, params, lmult, rmult, tmat, vmat, hmat, skew
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
    getcontroldim,

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


include(joinpath("util", "util.jl"))
include(joinpath("util", "custom_static.jl"))
include(joinpath("util", "customdict.jl"))
include(joinpath("util", "quaternion.jl"))

include(joinpath("graph", "entry.jl"))
include(joinpath("graph", "system.jl"))
include(joinpath("graph", "setup_functions.jl"))
include(joinpath("graph", "ldu.jl"))

include(joinpath("optional_components", "shapes.jl"))
include(joinpath("optional_components", "storage.jl"))

include(joinpath("main_components", "component.jl"))
include(joinpath("main_components", "state.jl"))
include(joinpath("main_components", "body.jl"))
include(joinpath("main_components", "abstractconstraint.jl"))
include(joinpath("main_components", "equalityconstraint.jl"))
include(joinpath("main_components", "inequalityconstraint.jl"))
include(joinpath("main_components", "controller.jl"))
include(joinpath("main_components", "mechanism_struct.jl"))
include(joinpath("main_components", "system.jl"))
include(joinpath("main_components", "mechanism_functions.jl"))

include(joinpath("joints", "abstract_joint.jl"))

include(joinpath("bounds", "bound.jl"))
include(joinpath("bounds", "contact_bounds.jl"))

include(joinpath("joints", "joint.jl"))
include(joinpath("joints", "translational.jl"))
include(joinpath("joints", "rotational.jl"))
include(joinpath("joints", "genericjoint.jl"))
include(joinpath("joints", "prototypes.jl"))

include(joinpath("solver", "solverfunctions.jl"))
include(joinpath("solver", "initconstraints.jl"))
include(joinpath("solver", "mehrotra.jl"))
include(joinpath("solver", "linesearch.jl"))

include(joinpath("discretization", "Linear.jl"))

include(joinpath("ui", "mechanism_ui.jl"))
include(joinpath("ui", "simulate.jl"))
include(joinpath("ui", "initialize.jl"))
include(joinpath("ui", "urdf.jl"))

include(joinpath("..", "examples", "mechanism_zoo.jl"))
include(joinpath("..", "examples", "diff_tools.jl"))
include(joinpath("joints", "force.jl"))
include(joinpath("joints", "torque.jl"))
include(joinpath("vis", "convertshape.jl"))
include(joinpath("vis", "visualize.jl"))

include(joinpath("optional_components", "energy.jl"))

end
