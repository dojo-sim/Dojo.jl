module Dojo

# constants
global const REG = 1.0e-10::Float64

#TODO: remove
using FiniteDiff

using LinearAlgebra
using Random
using StaticArrays
using SparseArrays
using StaticArrays: SUnitRange
using Quaternions
using Statistics

using Colors
using Colors: RGBA, RGB
using FFMPEG
using LightXML
using MeshCat
import MeshCat: render
using Meshing
using GeometryBasics
using GraphBasedSystems
using CoordinateTransformations

using DocStringExtensions

# Utilities
include(joinpath("utilities", "methods.jl"))
include(joinpath("utilities", "custom_static.jl"))
include(joinpath("utilities", "normalize.jl"))

# Orientation
include(joinpath("orientation", "quaternion.jl"))
include(joinpath("orientation", "mrp.jl"))
include(joinpath("orientation", "axis_angle.jl"))
include(joinpath("orientation", "mapping.jl"))
include(joinpath("orientation", "rotate.jl"))

# Graph objects
include(joinpath("mechanism", "node.jl"))
include(joinpath("mechanism", "edge.jl"))
include(joinpath("mechanism", "id.jl"))

# Bodies
include(joinpath("bodies", "shapes.jl"))
include(joinpath("bodies", "state.jl"))
include(joinpath("bodies", "constructor.jl"))
include(joinpath("bodies", "origin.jl"))
include(joinpath("bodies", "set.jl"))

# Mechanism
include(joinpath("joints", "constraints.jl"))
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
include(joinpath("contacts", "collisions", "collision.jl"))
include(joinpath("contacts", "collisions", "point_to_segment.jl"))
include(joinpath("contacts", "collisions", "point_to_box_v2.jl"))
include(joinpath("contacts", "collisions", "sphere_halfspace.jl"))
include(joinpath("contacts", "collisions", "sphere_sphere.jl"))
include(joinpath("contacts", "collisions", "sphere_capsule.jl"))
include(joinpath("contacts", "collisions", "sphere_box.jl"))
include(joinpath("contacts", "collisions", "string.jl"))
include(joinpath("contacts", "velocity.jl"))
include(joinpath("contacts", "impact.jl"))
include(joinpath("contacts", "linear.jl"))
include(joinpath("contacts", "nonlinear.jl"))
include(joinpath("contacts", "utilities.jl"))

# Solver
include(joinpath("solver", "linear_system.jl"))
include(joinpath("solver", "centering.jl"))
include(joinpath("solver", "complementarity.jl"))
include(joinpath("solver", "violations.jl"))
include(joinpath("solver", "options.jl"))
include(joinpath("solver", "initialization.jl"))
include(joinpath("solver", "correction.jl"))
include(joinpath("solver", "mehrotra.jl"))
include(joinpath("solver", "line_search.jl"))
include(joinpath("solver", "initialize_constraints.jl"))

# Integrator
include(joinpath("integrators", "integrator.jl"))
include(joinpath("integrators", "constraint.jl"))

# Visualizer
include(joinpath("visuals", "visualizer.jl"))
include(joinpath("visuals", "set.jl"))
include(joinpath("visuals", "convert.jl"))
include(joinpath("visuals", "colors.jl"))

# Data
include(joinpath("mechanism", "data.jl"))

# Gradients
include(joinpath("gradients", "contact.jl"))
include(joinpath("gradients", "finite_difference.jl"))
include(joinpath("gradients", "state.jl"))
include(joinpath("gradients", "data.jl"))
include(joinpath("gradients", "utilities.jl"))


# Bodies
export
    Body,
    Origin,
    Box,
    Capsule,
    Cylinder,
    Sphere,
    Pyramid,
    Mesh,
    CombinedShapes,
    get_body,
    get_node,
    set_external_force!,
    add_external_force!

# Joints
export
    Rotational,
    Translational,
    JointConstraint,
    Floating,
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
    get_joint,
    nullspace_mask

# Contacts
export
    ContactConstraint,
    ImpactContact,
    LinearContact,
    NonlinearContact,
    get_contact,
    get_sdf,
    contact_location,
    damper_impulses,
    contact_constraint

# Collision
export
    SphereHalfSpaceCollision,
    SphereSphereCollision,
    SphereCapsuleCollision,
    SphereBoxCollision,
    StringCollision,
    distance,
    contact_point,
    contact_normal,
    contact_tangent

# Inputs
export
    set_input!,
    add_input!,
    input_dimension

# Mechanism
export
    Mechanism,
    initialize!,
    set_floating_base,
    zero_coordinates!,
    zero_velocities!

# Maximal
export
    set_maximal_state!,
    set_maximal_configurations!,
    set_maximal_velocities!,
    get_maximal_state,
    get_maximal_gradients!,
    maximal_dimension

# Minimal
export
    set_minimal_state!,
    set_minimal_coordinates!,
    set_minimal_velocities!,
    get_minimal_state,
    get_minimal_gradients!,
    minimal_coordinates,
    minimal_velocities,
    minimal_dimension

# Maximal <-> Minimal
export
    maximal_to_minimal,
    minimal_to_maximal

# Simulation
export
    simulate!,
    step!,
    generate_storage

# Orientation
export
    Quaternion,
    attitude_jacobian

# Data
export
    get_data,
    set_data!,
    get_solution

# Gradients
export
    maximal_to_minimal_jacobian,
    minimal_to_maximal_jacobian

# Mechanics
export
    kinetic_energy,
    potential_energy,
    mechanical_energy,
    momentum

# Solver
export
    mehrotra!,
    SolverOptions

# Initialization
export 
    initialize_constraints!

# Linear System "Ax = b"
export
    full_matrix

# Visuals
export
    Visualizer,
    visualize,
    render,
    set_background!,
    set_floor!,
    set_surface!,
    set_light!,
    set_camera!,
    RGBA,
    orange,
    cyan,
    magenta

# Static
export
    szeros,
    sones,
    srand

# Utilities
export
    Storage

end
