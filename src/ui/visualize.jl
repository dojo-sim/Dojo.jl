function transform(x, q, shape)
    scale_transform = MeshCat.LinearMap(diagm(shape.scale))
    x_transform = MeshCat.Translation(x + vrotate(shape.xoffset, q))
    q_transform = MeshCat.LinearMap(q * shape.axis_offset)
    return MeshCat.compose(x_transform, q_transform, scale_transform)
end

MeshCat.js_scaling(s::AbstractVector) = s
MeshCat.js_position(p::AbstractVector) = p


function set_node!(x, q, id, shape, shapevisualizer, showshape) where {T,N}
    if showshape
        # TODO currently setting props directly because MeshCat/Rotations doesn't convert scaled rotation properly.
        # If this changes, do similarily to origin
        setprop!(shapevisualizer, "scale", MeshCat.js_scaling(shape.scale))
        setprop!(shapevisualizer, "position", MeshCat.js_position(x + vrotate(shape.xoffset, q)))
        setprop!(shapevisualizer, "quaternion", MeshCat.js_quaternion(q * shape.axis_offset))
    end
    return
end

function animate_node!(storage::Storage{T,N}, id, shape, animation, shapevisualizer, showshape) where {T,N}
    for i=1:N
        x = storage.x[id][i]
        q = storage.q[id][i]
        atframe(animation, i) do
            set_node!(x, q, id, shape, shapevisualizer, showshape)
        end
    end
    return
end

function MeshCat.setobject!(subvisshape, visshape, shapes::Shapes; transparent=false)
    for (i, s) in enumerate(shapes.shape)
        v = subvisshape["node_$i"]
        setobject!(v, visshape[i], s, transparent=transparent)
        scale_transform = MeshCat.LinearMap(diagm(s.scale))
        x_transform = MeshCat.Translation(s.xoffset)
        q_transform = MeshCat.LinearMap(s.axis_offset)
        t = MeshCat.compose(x_transform, q_transform, scale_transform)
        settransform!(v, t)
    end
end

function MeshCat.setobject!(subvisshape, visshape, shape::Shape; transparent=false)
    setobject!(subvisshape, visshape, MeshPhongMaterial(color=(transparent ? RGBA(0.75, 0.75, 0.75, 0.5) : shape.color)))
end

function MeshCat.setobject!(subvisshape, visshape::Vector, shape::Capsule; transparent=false)
    setobject!(subvisshape["cylinder"], visshape[1], MeshPhongMaterial(color=(transparent ? RGBA(0.75, 0.75, 0.75, 0.5) : shape.color)))
    setobject!(subvisshape["cap1"], visshape[2], MeshPhongMaterial(color=(transparent ? RGBA(0.75, 0.75, 0.75, 0.5) : shape.color)))
    setobject!(subvisshape["cap2"], visshape[3], MeshPhongMaterial(color=(transparent ? RGBA(0.75, 0.75, 0.75, 0.5) : shape.color)))
end

function MeshCat.setobject!(subvisshape, visshape, shape::Mesh; transparent=false)
    if visshape.mtl_library == ""
        visshape = MeshFileGeometry(visshape.contents, visshape.format)
        setobject!(subvisshape, visshape, MeshPhongMaterial(color=shape.color))
    else
        setobject!(subvisshape, visshape)
    end
end

function build_robot(mechanism::Mechanism; vis::Visualizer=Visualizer(),
        show_contact=false, name::Symbol=:robot, color=nothing) where {T,N}

    bodies = mechanism.bodies
    origin = mechanism.origin
    set_background!(vis)
    set_light!(vis)
    set_floor!(vis)

    # Bodies and Contacts
    for (id,body) in enumerate(bodies)
        if color !== nothing
            shape = deepcopy(body.shape)
            set_color!(shape, color)
        else
            shape = body.shape
        end
        visshape = convert_shape(shape)
        subvisshape = nothing
        if visshape !== nothing
            subvisshape = vis[name][:bodies][Symbol(body.name, "_id_$id")]
            setobject!(subvisshape,visshape,shape,transparent=show_contact)
        end

        if show_contact
            for (jd, contact) in enumerate(mechanism.contacts)
                if contact.parent_id == body.id
                    radius = abs(contact.model.offset[3])
                    (radius == 0.0) && (radius = 0.01)
                    contact_shape = Sphere(radius,
                        xoffset=(contact.model.contact_point),
                        axis_offset=one(UnitQuaternion), color=RGBA(1.0, 0.0, 0.0, 0.5))
                    visshape = convert_shape(contact_shape)
                    subvisshape = nothing
                    if visshape !== nothing
                        subvisshape = vis[name][:contacts][Symbol(contact.name, "_$(jd)")]
                        setobject!(subvisshape,visshape,contact_shape,transparent=false)
                    end
                end
            end
        end
    end

    # Origin
    id = origin.id
    shape = origin.shape
    visshape = convert_shape(shape)
    if visshape !== nothing
        subvisshape = vis[name][:bodies][Symbol(:origin, "_id")]
        setobject!(subvisshape,visshape,shape,transparent=show_contact)
    end
    return vis
end

function set_robot(vis::Visualizer, mechanism::Mechanism, z::Vector{T};
        show_contact::Bool=true, name::Symbol=:robot) where {T,N}
    (length(z) == minimal_dimension(mechanism)) && (z = minimal_to_maximal(mechanism, z))
    bodies = mechanism.bodies
    origin = mechanism.origin

    # Bodies and Contacts
    for (id,body) in enumerate(bodies)
        x, _, q, _ = unpack_maximal_state(z, id)
        shape = body.shape
        visshape = convert_shape(shape)
        subvisshape = nothing
        showshape = false
        if visshape !== nothing
            subvisshape = vis[name][:bodies][Symbol(body.name, "_id_$id")]
            showshape = true
        end
        set_node!(x, q, id, shape, subvisshape, showshape)

        if show_contact
            for (jd, contact) in enumerate(mechanism.contacts)
                if contact.parent_id == body.id
                    radius = abs(contact.model.offset[3])
                    (radius == 0.0) && (radius = 0.01)
                    contact_shape = Sphere(radius,
                        xoffset=(contact.model.contact_point),
                        axis_offset=one(UnitQuaternion), color=RGBA(1.0, 0.0, 0.0, 0.5))
                    visshape = convert_shape(contact_shape)
                    subvisshape = nothing
                    showshape = false
                    if visshape !== nothing
                        subvisshape = vis[name][:contacts][Symbol(contact.name, "_$(jd)")]
                        showshape = true
                    end
                    set_node!(x, q, id, contact_shape, subvisshape, showshape)
                end
            end
        end
    end

    # Origin
    id = origin.id
    shape = origin.shape
    visshape = convert_shape(shape)
    if visshape !== nothing
        subvisshape = vis[name][:bodies][Symbol(:origin, "_id")]
        shapetransform = transform(szeros(T,3), one(UnitQuaternion{T}), shape)
        settransform!(subvisshape, shapetransform)
    end
    return vis
end

function visualize(mechanism::Mechanism, storage::Storage{T,N}; vis::Visualizer=Visualizer(),
        build::Bool=true, show_contact=false, animation=nothing, color=nothing, name::Symbol=:robot) where {T,N}

    storage = deepcopy(storage)
    bodies = mechanism.bodies
    origin = mechanism.origin

    # Build robot in the visualizer
    build && build_robot(mechanism, vis=vis, show_contact=show_contact, color=color, name=name)

    # Create animations
    framerate = Int64(round(1/mechanism.timestep))
    (animation == nothing) && (animation =
        MeshCat.Animation(Dict{MeshCat.SceneTrees.Path,MeshCat.AnimationClip}(), framerate))

    # Bodies and Contacts
    for (id,body) in enumerate(bodies)
        shape = body.shape
        visshape = convert_shape(shape)
        subvisshape = nothing
        showshape = false
        if visshape !== nothing
            subvisshape = vis[name][:bodies][Symbol(body.name, "_id_$id")]
            showshape = true
        end
        animate_node!(storage, id, shape, animation, subvisshape, showshape)

        if show_contact
            for (jd, contact) in enumerate(mechanism.contacts)
                if contact.parent_id == body.id
                    radius = abs(contact.model.offset[3])
                    (radius == 0.0) && (radius = 0.01)
                    contact_shape = Sphere(radius,
                        xoffset=(contact.model.contact_point),
                        axis_offset=one(UnitQuaternion), color=RGBA(1.0, 0.0, 0.0, 0.5))
                    visshape = convert_shape(contact_shape)
                    subvisshape = nothing
                    showshape = false
                    if visshape !== nothing
                        subvisshape = vis[name][:contacts][Symbol(contact.name, "_$(jd)")]
                        showshape = true
                    end
                    animate_node!(storage, id, contact_shape, animation, subvisshape, showshape)
                end
            end
        end
    end

    # Origin
    id = origin.id
    shape = origin.shape
    visshape = convert_shape(shape)
    if visshape !== nothing
        subvisshape = vis[name][:bodies][Symbol(:origin, "_id")]
        shapetransform = transform(szeros(T,3), one(UnitQuaternion{T}), shape)
        settransform!(subvisshape, shapetransform)
    end

    setanimation!(vis, animation)
    return vis, animation
end

#
# function visualize_maximal(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, z::AbstractVector, vis::Visualizer) where {T,Nn,Ne,Nb,Ni}
# 	storage = Storage(1,Nb)
# 	for t = 1:1
# 		for b = 1:Nb
# 			x2, v15, q2, ϕ15 = unpack_maximal_state(z, b)
# 			storage.x[b][t] = x2
# 			storage.v[b][t] = v15
# 			storage.q[b][t] = q2
# 			storage.ω[b][t] = ϕ15
# 		end
# 	end
# 	visualize(mechanism, storage, vis = vis)
# end


# function prepare_vis!(storage::Storage{T,N}, id, shape, animation, shapevisualizer, showshape) where {T,N}
#     if showshape
#         for i=1:N
#             x = storage.x[id][i]
#             q = storage.q[id][i]
#             # TODO currently setting props directly because MeshCat/Rotations doesn't convert scaled rotation properly.
#             # If this changes, do similarily to origin
#             atframe(animation, i) do
#                 setprop!(shapevisualizer, "scale", MeshCat.js_scaling(shape.scale))
#                 setprop!(shapevisualizer, "position", MeshCat.js_position(x + vrotate(shape.xoffset, q)))
#                 setprop!(shapevisualizer, "quaternion", MeshCat.js_quaternion(q * shape.axis_offset))
#             end
#         end
#     end
#     return
# end
#
# function set_robot(vis::Visualizer, mechanism::Mechanism, z::Vector{T}; name::Symbol=:robot) where {T,N}
#     bodies = mechanism.bodies
#     origin = mechanism.origin
#
#     i = 1
#     for (id,body) in enumerate(bodies)
#         shape = body.shape
#         visshape = convert_shape(shape)
#         subvisshape = vis[name]["bodies/body:"*string(id)]
#
#         x = z[(i-1) * 13 .+ (1:3)]
#         q = UnitQuaternion(z[(i-1) * 13 + 6 .+ (1:4)]...)
#
#         if visshape !== nothing
#             setprop!(subvisshape, "scale", MeshCat.js_scaling(shape.scale))
#             setprop!(subvisshape, "position", MeshCat.js_position(x + vrotate(shape.xoffset, q)))
#             setprop!(subvisshape, "quaternion", MeshCat.js_quaternion(q * shape.axis_offset))
#         end
#         i += 1
#     end
#
#     id = origin.id
#     shape = origin.shape
#     visshape = convert_shape(shape)
#     subvisshape = vis[name]["bodies/origin:"*string(id)]
#     if visshape !== nothing
#         shapetransform = transform(szeros(T,3), one(UnitQuaternion{T}), shape)
#         settransform!(subvisshape, shapetransform)
#     end
#
#     return vis
# end

# function visualize(mechanism::Mechanism, storage::Storage{T,N}; vis::Visualizer=Visualizer(),
#         showframes::Bool=false, show_contact=false, animation=nothing, name::Symbol=:robot) where {T,N}
#
#     storage = deepcopy(storage)
#     bodies = mechanism.bodies
#     origin = mechanism.origin
#     if showframes
#         triads = [Triad(0.33) for i=1:length(bodies)]
#     end
#
#     set_background!(vis)
#     set_light!(vis)
#     set_floor!(vis)
#
#     # Create animations
#     framerate = Int64(round(1/mechanism.timestep))
#     (animation == nothing) && (animation =
#         MeshCat.Animation(Dict{MeshCat.SceneTrees.Path,MeshCat.AnimationClip}(), framerate))
#
#     # Bodies and Contacts
#     for (id,body) in enumerate(bodies)
#         shape = body.shape
#         visshape = convert_shape(shape)
#         subvisshape = nothing
#         subvisframe = nothing
#         showshape = false
#         if visshape !== nothing
#             subvisshape = vis[name][:bodies][Symbol(body.name, "_id_$id")]
#             setobject!(subvisshape,visshape,shape,transparent=show_contact)
#             showshape = true
#         end
#
#         prepare_vis!(storage, id, shape, animation, subvisshape, subvisframe, showshape, showframes)###########
#
#         if show_contact
#             for (jd, contact) in enumerate(mechanism.contacts)
#                 if contact.parent_id == body.id
#                     radius = abs(contact.model.offset[3])
#                     (radius == 0.0) && (radius = 0.01)
#                     contact_shape = Sphere(radius,
#                         xoffset=(contact.model.contact_point),
#                         axis_offset=one(UnitQuaternion), color=RGBA(1.0, 0.0, 0.0, 0.5))
#                     visshape = convert_shape(contact_shape)
#                     subvisshape = nothing
#                     subvisframe = nothing
#                     showshape = false
#                     if visshape !== nothing
#                         subvisshape = vis[name][:contacts][Symbol(contact.name, "_$(jd)")]
#                         setobject!(subvisshape,visshape,contact_shape,transparent=false)
#                         showshape = true
#                     end
#                     prepare_vis!(storage, id, contact_shape, animation, subvisshape, subvisframe, showshape, showframes)################
#                 end
#             end
#         end
#     end
#
#     # Origin
#     id = origin.id
#     shape = origin.shape
#     visshape = convert_shape(shape)
#     if visshape !== nothing
#         subvisshape = vis[name][:bodies][Symbol(:origin, "_id")]
#         setobject!(subvisshape,visshape,shape,transparent=show_contact)
#         shapetransform = transform(szeros(T,3), one(UnitQuaternion{T}), shape)################################
#         settransform!(subvisshape, shapetransform)################################
#     end
#
#     setanimation!(vis, animation)
#     return vis, animation
# end


# function build_robot(vis::Visualizer, mechanism::Mechanism; name::Symbol=:robot, color=nothing) where {T,N}
#
#     # bodies = mechanism.bodies
#     # origin = mechanism.origin
#     #
#     # setprop!(vis["/Background"], "top_color", RGBA(1.0, 1.0, 1.0))
#     # setprop!(vis["/Background"], "bottom_color", RGBA(1.0, 1.0, 1.0))
#     # setvisible!(vis["/Axes"],false)
#
#     for (id,body) in enumerate(bodies)
#         shape = deepcopy(body.shape)
#
#         # if color !== nothing
#         #     if shape isa Shapes
#         #         for i = 1:length(shape.shape)
#         #             shape.shape[i].color = color
#         #         end
#         #     else
#         #         if shape isa EmptyShape
#         #             nothing
#         #         else
#         #             shape.color = color
#         #         end
#         #     end
#         end
#         visshape = convert_shape(shape)
#         subvisshape = nothing
#         subvisframe = nothing
#         showshape = false
#         if visshape !== nothing
#             subvisshape = vis[name]["bodies/body:"*string(id)]
#             setobject!(subvisshape,visshape,shape)
#             showshape = true
#         end
#     end
#
#     id = origin.id
#     shape = origin.shape
#     visshape = convert_shape(shape)
#     if visshape !== nothing
#         subvisshape = vis[name]["bodies/origin:"*string(id)]
#         setobject!(subvisshape,visshape,shape)
#     end
#
#    return vis
# end
