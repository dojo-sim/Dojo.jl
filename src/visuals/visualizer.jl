"""
    visualize(mechanism, storage; vis, build, show_contact, animation, color, name)

    visualize mechanism using trajectory from storage 

    mechanism: Mechanism 
    storage: Storage 
    vis: Visualizer 
    build: flag to construct mechanism visuals (only needs to be built once)
    show_contact: flag to show contact locations on system 
    color: RGBA 
    name: unique identifier for mechanism
"""
function visualize(mechanism::Mechanism, storage::Storage{T,N}; vis::Visualizer=Visualizer(),
    build::Bool=true, 
    show_joint=false,
    joint_radius=0.1,
    show_contact=false,
    show_frame=false, 
    animation=nothing, 
    color=nothing, 
    name::Symbol=:robot,
    return_animation=false,
    visualize_floor=true) where {T,N}

    storage = deepcopy(storage)
    bodies = mechanism.bodies
    origin = mechanism.origin

    # Build robot in the visualizer
    build && build_robot(mechanism; 
        vis, show_joint, show_contact, show_frame, color, name, visualize_floor)

    # Create animations
    framerate = Int64(round(1/mechanism.timestep))
    (animation === nothing) && (animation =
        MeshCat.Animation(Dict{MeshCat.SceneTrees.Path,MeshCat.AnimationClip}(), framerate))

    # Bodies and Contacts
    for (id, body) in enumerate(bodies)
        shape = body.shape
        visshape = convert_shape(shape)
        subvisshape = nothing
        showshape = false
        if visshape !== nothing
            subvisshape = vis[name][:bodies][Symbol(body.name, "__id_$id")]
            showshape = true
        end

        animate_node!(storage, id, shape, animation, subvisshape, showshape)

        if show_joint
            for (jd, joint) in enumerate(mechanism.joints)
                if joint.child_id == body.id
                    radius = 0.1
                    joint_shape = Sphere(radius, 
                        position_offset=joint.translational.vertices[2],
                        color=RGBA(0.0, 0.0, 1.0, 0.5))
                    visshape = convert_shape(joint_shape)
                    subvisshape = nothing
                    showshape = false
                    if visshape !== nothing
                        subvisshape = vis[name][:joints][Symbol(joint.name, "__id_$(jd)")]
                        showshape = true
                    end
                    animate_node!(storage, id, joint_shape, animation, subvisshape, showshape)
                end
            end
        end

        if show_contact
            for (jd, contact) in enumerate(mechanism.contacts)
                if contact.parent_id == body.id
                    radius = abs(contact.model.collision.contact_radius)
                    (radius == 0.0) && (radius = 0.01)
                    contact_shape = Sphere(radius,
                        position_offset=contact.model.collision.contact_origin, #TODO: generalize for collision checking
                        orientation_offset=one(Quaternion), 
                        color=RGBA(1.0, 0.0, 0.0, 0.5))
                    visshape = convert_shape(contact_shape)
                    subvisshape = nothing
                    showshape = false
                    if visshape !== nothing
                        subvisshape = vis[name][:contacts][Symbol(contact.name, "__id_$(jd)")]
                        showshape = true
                    end
                    animate_node!(storage, id, contact_shape, animation, subvisshape, showshape)
                end
            end
        end

        if show_frame
            frame_shape = FrameShape(scale=0.33*ones(3))
            visshape = convert_shape(frame_shape)
            subvisshape = vis[name][:frames][Symbol(body.name, "__id_$id")]
            showshape = true
            animate_node!(storage, id, frame_shape, animation, subvisshape, showshape)
        end
    end

    # Origin
    id = origin.id
    shape = origin.shape
    visshape = convert_shape(shape)
    if visshape !== nothing
        subvisshape = vis[name][:bodies][Symbol(:origin, "_id")]
        shapetransform = transform(szeros(T,3), one(Quaternion{T}), shape)
        settransform!(subvisshape, shapetransform)
    end

    setanimation!(vis, animation)
    return_animation ? (return vis, animation) : (return vis) 
end

"""
    build_robot(mechanism; vis, show_contact, name, color)

    construct visuals for mechanism 

    mechanism: Mechanism 
    vis: Visualizer 
    show_contact: flag to show contact locations on mechanism 
    name: unique identifier 
    color: RGBA
"""
function build_robot(mechanism::Mechanism; 
    vis::Visualizer=Visualizer(),
    show_joint=false,
    joint_radius=0.1,
    show_contact=false, 
    show_frame=false,
    name::Symbol=:robot, 
    color=nothing,
    visualize_floor=true)

    bodies = mechanism.bodies
    origin = mechanism.origin
    set_background!(vis)
    set_light!(vis)
    visualize_floor && set_floor!(vis)

    # Bodies and Contacts
    for (id, body) in enumerate(bodies)
        if color !== nothing
            shape = deepcopy(body.shape)
            set_color!(shape, color)
        else
            shape = body.shape
        end
        visshape = convert_shape(shape)
        subvisshape = nothing
        if visshape !== nothing
            subvisshape = vis[name][:bodies][Symbol(body.name, "__id_$id")]
            setobject!(subvisshape, visshape, shape, 
                transparent=(show_joint || show_contact))
        end

        if show_joint
            for (jd, joint) in enumerate(mechanism.joints)
                if joint.child_id == body.id
                    radius = joint_radius
                    joint_shape = Sphere(radius,
                        position_offset=joint.translational.vertices[2],
                        color=RGBA(0.0, 0.0, 1.0, 0.5))
                    visshape = convert_shape(joint_shape)
                    subvisshape = nothing
                    if visshape !== nothing
                        subvisshape = vis[name][:joints][Symbol(joint.name, "__id_$(jd)")]
                        setobject!(subvisshape, visshape, joint_shape, 
                            transparent=false)
                    end
                end
            end
        end

        if show_contact
            for (jd, contact) in enumerate(mechanism.contacts)
                if contact.parent_id == body.id
                    radius = abs(contact.model.collision.contact_radius)
                    (radius == 0.0) && (radius = 0.01)
                    contact_shape = Sphere(radius,
                        position_offset=(contact.model.collision.contact_origin),
                        orientation_offset=one(Quaternion), color=RGBA(1.0, 0.0, 0.0, 0.5))
                    visshape = convert_shape(contact_shape)
                    subvisshape = nothing
                    if visshape !== nothing
                        subvisshape = vis[name][:contacts][Symbol(contact.name, "__id_$(jd)")]
                        setobject!(subvisshape,visshape,contact_shape,transparent=false)
                    end
                end
            end
        end

        if show_frame
            frame_shape = FrameShape(scale=0.33*ones(3))
            visshape = convert_shape(frame_shape)
            subvisshape = vis[name][:frames][Symbol(body.name, "__id_$id")]
            setobject!(subvisshape, visshape, frame_shape, 
                transparent=false)
        end
    end

    # Origin
    id = origin.id
    shape = origin.shape
    visshape = convert_shape(shape)
    if visshape !== nothing
        subvisshape = vis[name][:bodies][Symbol(:origin, "_id")]
        setobject!(subvisshape,visshape,shape,transparent=(show_joint || show_contact))
    end
    return vis
end

"""
    set_robot(vis, mechanism, z; show_contact, name)

    visualze mechanism configuration from maximal representation 

    vis: Visualizer 
    mechanism: Mechanism 
    z: maximal state 
    show_contact: flag to show contact locations on mechanism 
    name: unique identifier
"""
function set_robot(vis::Visualizer, mechanism::Mechanism, z::Vector{T};
    show_joint::Bool=false,
    joint_radius=0.1,
    show_contact::Bool=true, 
    name::Symbol=:robot) where {T}

    (length(z) == minimal_dimension(mechanism)) && (z = minimal_to_maximal(mechanism, z))
    bodies = mechanism.bodies
    origin = mechanism.origin

    # Bodies and Contacts
    for (id, body) in enumerate(bodies)
        x, _, q, _ = unpack_maximal_state(z, id)
        shape = body.shape
        visshape = convert_shape(shape)
        subvisshape = nothing
        showshape = false
        if visshape !== nothing
            subvisshape = vis[name][:bodies][Symbol(body.name, "__id_$id")]
            showshape = true
        end

        set_node!(x, q, id, shape, subvisshape, showshape)

        if show_joint
            for (jd, joint) in enumerate(mechanism.joints)
                if joint.child_id == body.id
                    radius = joint_radius
                    joint_shape = Sphere(radius,
                        position_offset=joint.translational.vertices[2],
                        color=RGBA(0.0, 0.0, 1.0, 0.5))
                    visshape = convert_shape(joint_shape)
                    subvisshape = nothing
                    showshape = false
                    if visshape !== nothing
                        subvisshape = vis[name][:joints][Symbol(joint.name, "__id_$(jd)")]
                        showshape = true
                    end
                    set_node!(x, q, id, joint_shape, subvisshape, showshape)
                end
            end
        end

        if show_contact
            for (jd, contact) in enumerate(mechanism.contacts)
                if contact.parent_id == body.id
                    radius = abs(contact.model.collision.contact_radius)
                    (radius == 0.0) && (radius = 0.01)
                    contact_shape = Sphere(radius,
                        position_offset=(contact.model.collision.contact_origin),
                        orientation_offset=one(Quaternion), color=RGBA(1.0, 0.0, 0.0, 0.5))
                    visshape = convert_shape(contact_shape)
                    subvisshape = nothing
                    showshape = false
                    if visshape !== nothing
                        subvisshape = vis[name][:contacts][Symbol(contact.name, "__id_$(jd)")]
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
        shapetransform = transform(szeros(T,3), one(Quaternion{T}), shape)
        settransform!(subvisshape, shapetransform)
    end
    return vis
end

function transform(x, q, shape)
    scale_transform = MeshCat.LinearMap(diagm(shape.scale))
    x_transform = MeshCat.Translation(x + vector_rotate(shape.position_offset, q))
    q_transform = MeshCat.LinearMap(q * shape.orientation_offset)
    return MeshCat.compose(x_transform, q_transform, scale_transform)
end

MeshCat.js_scaling(s::AbstractVector) = s
MeshCat.js_position(p::AbstractVector) = p

function set_node!(x, q, id, shape, shapevisualizer, showshape)
    if showshape
        # TODO currently setting props directly because MeshCat/Rotations doesn't convert scaled rotation properly.
        # If this changes, do similarily to origin
        setprop!(shapevisualizer, "scale", MeshCat.js_scaling(shape.scale))
        setprop!(shapevisualizer, "position", MeshCat.js_position(x + vector_rotate(shape.position_offset, q)))
        setprop!(shapevisualizer, "quaternion", MeshCat.js_quaternion(q * shape.orientation_offset))
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

function MeshCat.setobject!(subvisshape, visshapes::Vector, shape::CombinedShapes; transparent=false)
    for (i,visshape) in enumerate(visshapes) 
        v = subvisshape["shape"*string(i)]
        s = shape.shapes[i]
        setobject!(v, visshape, s; transparent)
        scale_transform = MeshCat.LinearMap(diagm(s.scale))
        x_transform = MeshCat.Translation(s.position_offset)
        q_transform = MeshCat.LinearMap(s.orientation_offset)
        t = MeshCat.compose(x_transform, q_transform, scale_transform)
        settransform!(v, t)
    end
end

function MeshCat.setobject!(subvisshape, visshape, shape::Shape; transparent=false)
    setobject!(subvisshape, visshape, MeshPhongMaterial(color=(transparent ? RGBA(0.75, 0.75, 0.75, 0.5) : shape.color)))
end

function MeshCat.setobject!(subvisshape, visshape, shape::FrameShape; transparent=false)
    setobject!(subvisshape, visshape)
end

function MeshCat.setobject!(subvisshape, visshape, shape::Mesh; transparent=false)
    if visshape.mtl_library == ""
        visshape = MeshFileGeometry(visshape.contents, visshape.format)
        setobject!(subvisshape, visshape, MeshPhongMaterial(color=(transparent ? RGBA(0.75, 0.75, 0.75, 0.5) : shape.color)))
    else
        setobject!(subvisshape, visshape)
    end
end
