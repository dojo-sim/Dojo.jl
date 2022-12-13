### Before parsing

unsafeattribute(::Nothing, ::Core.AbstractString) = nothing
unsafeattribute(x::LightXML.XMLElement, name::Core.AbstractString) = attribute(x, name)

function parse_scalar(xel, name::String, T; default::Union{String,Nothing}=nothing)
    scalarstr = unsafeattribute(xel, name)
    if scalarstr === nothing
        if default === nothing
            @error "no parsable scalar found"
        else
            scalarstr = default
        end
    end

    return parse(T, scalarstr)
end

function parse_vector(xel, name::String, T; default::Union{String,Nothing}=nothing)
    vectorstr = unsafeattribute(xel, name)
    if vectorstr === nothing
        if default === nothing
            @error "no parsable vector found"
        else
            vectorstr = default
        end
    end

    return parse.(T,split(vectorstr))
end

function parse_inertiamatrix(xinertia, T)
    if xinertia === nothing
        J = zeros(T, 3, 3)
    else
        ixx = parse_scalar(xinertia, "ixx", T)
        ixy = parse_scalar(xinertia, "ixy", T; default = "0")
        ixz = parse_scalar(xinertia, "ixz", T; default = "0")
        iyy = parse_scalar(xinertia, "iyy", T)
        iyz = parse_scalar(xinertia, "iyz", T; default = "0")
        izz = parse_scalar(xinertia, "izz", T)
        J = [ixx ixy ixz; ixy iyy iyz; ixz iyz izz]
    end

    return J
end

function parse_pose(xpose, T)
    if xpose === nothing
        x, q = zeros(T, 3), one(Quaternion{T})
    else
        x = parse_vector(xpose, "xyz", T; default = "0 0 0")
        rpy = parse_vector(xpose, "rpy", T; default = "0 0 0")
        q = RotZ(rpy[3])*RotY(rpy[2])*RotX(rpy[1])
    end

    return x, q
end

function parse_inertia(xinertial, T)
    if xinertial === nothing
        x = zeros(T, 3)
        q = one(Quaternion{T})
        m = zero(T)
        J = zeros(T, 3, 3)
    else
        x, q = parse_pose(find_element(xinertial, "origin"), T)
        J = parse_inertiamatrix(find_element(xinertial, "inertia"), T)
        m = parse_scalar(find_element(xinertial, "mass"), "value", T; default = "0")
    end

    return x, q, m, J
end

function parse_robotmaterials(xroot, T)
    xmaterials = get_elements_by_tagname(xroot, "material")

    mdict = Dict{String,Vector{T}}()

    for xmaterial in xmaterials
        name = attribute(xmaterial, "name")
        colornode = find_element(xmaterial, "color")
        colorvec = parse_vector(colornode, "rgba", T)
        mdict[name] = colorvec
    end

    return mdict
end

function parse_xmaterial(xmaterial, materialdict, T)
    if xmaterial === nothing
        color = RGBA(0.75, 0.75, 0.75)
    else
        name = unsafeattribute(xmaterial, "name")
        if haskey(materialdict, name)
            colorvec = materialdict[name]
        else
            colorvec = [0.75; 0.75; 0.75; 1.0]
        end
        colornode = find_element(xmaterial, "color")
        colorvec = parse_vector(colornode, "rgba", T; default = string(colorvec[1]," ",colorvec[2]," ",colorvec[3]," ",colorvec[4]))

        color = RGBA(colorvec...)
    end

    return color
end

# function parse_shape(xvisual, materialdict, T, xb, qb)
function parse_shape(xvisual, materialdict, T; path_prefix)

    if xvisual === nothing
        shape = nothing
    else
        xgeometry = find_element(xvisual, "geometry")
        @assert xgeometry !== nothing

        color = parse_xmaterial(find_element(xvisual, "material"), materialdict, T)
        x, q = parse_pose(find_element(xvisual, "origin"), T)

        shapenodes = LightXML.XMLElement[]
        for node in child_nodes(xgeometry)  # node is an instance of XMLNode
            if is_elementnode(node)
                push!(shapenodes, XMLElement(node))
            end
        end

        if length(shapenodes) == 0
            shape = nothing
        else
            if length(shapenodes) > 1
                @warn "Multiple geometries."
            end

            shapenode = shapenodes[1]
            shape = get_shape(shapenode, x, q, color, T; path_prefix)
        end
    end
    return shape
end

function get_shape(shapenode, x, q, color, T; path_prefix)
    if name(shapenode) == "box"
        xyz = parse_vector(shapenode, "size", T; default = "1 1 1")
        shape = Box(xyz..., zero(T); color, position_offset = x, orientation_offset = q)
    elseif name(shapenode) == "cylinder"
        r = parse_scalar(shapenode, "radius", T; default = "0.5")
        l = parse_scalar(shapenode, "length", T; default = "1")
        shape = Cylinder(r, l, zero(T); color, position_offset = x, orientation_offset = q)
    elseif name(shapenode) == "pyramid"
        w = parse_scalar(shapenode, "width", T; default = "1")
        h = parse_scalar(shapenode, "height", T; default = "1")
        shape = Cylinder(w, h, zero(T); color, position_offset = x, orientation_offset = q)
    elseif name(shapenode) == "sphere"
        r = parse_scalar(shapenode, "radius", T; default = "0.5")
        shape = Sphere(r, zero(T); color, position_offset = x, orientation_offset = q)
    elseif name(shapenode) == "mesh"
        path = attribute(shapenode, "filename")
        scale = parse_vector(shapenode, "scale", T; default = "1 1 1")
        shape = Mesh(normpath(joinpath(path_prefix, path)), zero(T), zeros(T, 3, 3); scale, color, position_offset = x, orientation_offset = q)
    elseif name(shapenode) == "capsule"
        r = parse_scalar(shapenode, "radius", T; default = "0.5")
        l = parse_scalar(shapenode, "length", T; default = "1")
        shape = Capsule(r, l, zero(T); color, position_offset = x, orientation_offset = q)
    else
        @info "Unknown geometry."
        shape = nothing
    end
end

function parse_link(xlink, materialdict, T; path_prefix)
    x, q, m, J = parse_inertia(find_element(xlink, "inertial"), T)
    xvisuals = get_elements_by_tagname(xlink, "visual")
    shapes = [parse_shape(xvisual, materialdict, T; path_prefix) for xvisual in xvisuals]
    # shapes = [parse_shape(xvisual, materialdict, T, x, q) for xvisual in xvisuals]
    if length(shapes) == 0
        shape = nothing
    elseif length(shapes) > 1
        s = [s.shape for s in shapes]
        shape = CombinedShapes(s, 0.0, diagm([0.0; 0.0; 0.0]))
    else
        shape = shapes[1]
    end
    name = attribute(xlink, "name")

    if shape === nothing
        link = Body(m, J, name=Symbol(name))
    else
        link = shape
        link.mass = m
        link.inertia = J
        link.name = Symbol(name)
    end

    link.state.x2 = x
    link.state.q2 = q


    return link
end

function parse_links(xlinks, materialdict, T; path_prefix)
    ldict = Dict{Symbol,Body{T}}()

    for xlink in xlinks
        link = parse_link(xlink, materialdict, T; path_prefix)
        ldict[link.name] = link
    end

    return ldict
end

# TODO: fix axis
function joint_selector(joint_type, pbody, cbody, T;
        axis = SA{T}[1;0;0], parent_vertex = szeros(T,3), child_vertex = szeros(T,3), orientation_offset = one(Quaternion{T}), name = Symbol("joint_" * randstring(4)), damper = zero(T))

    # TODO @warn "this is not great"
    axis = vector_rotate(axis, orientation_offset) # inv(orientation_offset) * axis

    # TODO limits for revolute joint?
    if joint_type == "revolute" || joint_type == "continuous"
        joint = JointConstraint(Revolute(pbody, cbody, axis; parent_vertex, child_vertex, orientation_offset, damper); name)
    elseif joint_type == "prismatic"
        joint = JointConstraint(Prismatic(pbody, cbody, axis; parent_vertex, child_vertex, orientation_offset, damper); name)
    elseif joint_type == "planar"
        joint = JointConstraint(Planar(pbody, cbody, axis; parent_vertex, child_vertex, orientation_offset, damper); name)
    elseif joint_type == "planarfree"
        joint = JointConstraint(PlanarFree(pbody, cbody, axis; parent_vertex, child_vertex, damper); name)
    elseif joint_type == "fixed"
        joint = JointConstraint(Fixed(pbody, cbody; parent_vertex, child_vertex, orientation_offset); name)
    elseif joint_type == "floating"
        joint = JointConstraint(Floating(pbody, cbody; damper); name)
    elseif joint_type == "orbital"
        joint = JointConstraint(Orbital(pbody, cbody, axis; parent_vertex, child_vertex, orientation_offset, damper); name)
    elseif joint_type == "ball"
        joint = JointConstraint(Spherical(pbody, cbody; parent_vertex, child_vertex, orientation_offset, damper); name)
    elseif joint_type == "fixedorientation"
        joint = JointConstraint(FixedOrientation(pbody, cbody; orientation_offset, damper); name)
    elseif joint_type == "cylindrical"
        joint = JointConstraint(Cylindrical(pbody, cbody, axis; parent_vertex, child_vertex, orientation_offset, damper); name)
    elseif joint_type == "cylindricalfree"
        joint = JointConstraint(CylindricalFree(pbody, cbody, axis; parent_vertex, child_vertex, damper); name)
    elseif joint_type == "planaraxis"
        joint = JointConstraint(PlanarAxis(pbody, cbody, axis; parent_vertex, child_vertex, orientation_offset, damper); name)
    else
        @error "Unknown joint type"
    end

    return joint
end

function parse_joint(xjoint, plink, clink, T, parse_damper)
    joint_type = attribute(xjoint, "type")
    x, q = parse_pose(find_element(xjoint, "origin"), T)
    axis = parse_vector(find_element(xjoint, "axis"), "xyz", T; default = "1 0 0")
    parent_vertex = x
    name = Symbol(attribute(xjoint, "name"))
    parse_damper ? damper = parse_scalar(find_element(xjoint, "dynamics"), "damping", T; default = "0") : damper = 0

    return joint_selector(joint_type, plink, clink, T; axis, parent_vertex, orientation_offset = q, name, damper)
end

function parse_loop_joint(xjoint, pbody, cbody, T, parse_damper)
    find_element(xjoint, "link1")
    find_element(xjoint, "link2")

    joint_type = attribute(xjoint, "type")
    axis = parse_vector(find_element(xjoint, "axis"), "xyz", T, default = "1 0 0")
    x1, q1 = parse_pose(find_element(xjoint, "link1"), T)
    x2, _ = parse_pose(find_element(xjoint, "link2"), T) # The orientation q2 of the second body is ignored because it is determined by the mechanism's structure
    parent_vertex = x1
    child_vertex = x2
    name = Symbol(attribute(xjoint, "name"))
    parse_damper ? damper = parse_scalar(find_element(xjoint, "dynamics"), "damping", T; default = "0") : damper = 0

    return joint_selector(joint_type, pbody, cbody, T; axis, parent_vertex, child_vertex, orientation_offset = q1, name, damper)
end

function parse_joints(xjoints, ldict, floating, T, parse_damper)
    origins = Origin{T}[]
    links = Body{T}[]
    joints = JointConstraint{T}[]
    floatingname = ""

    for name in keys(ldict)
        childflag = false
        for xjoint in xjoints
            xchild = find_element(xjoint, "child")
            childname = attribute(xchild, "link")
            if Symbol(childname) == name
                childflag = true
                break
            end
        end
        if childflag
            push!(links, ldict[name])
        else
            origin = Origin{T}()
            if floating # keep current link and create new origin
                push!(links, ldict[name])
                floatingname = name
            else # make current link origin
                origin = Origin(ldict[name])
            end
            push!(origins, origin)
        end
    end

    @assert length(origins) == 1 "Multiple origins"
    origin = origins[1]

    for xjoint in xjoints
        xplink = find_element(xjoint, "parent")
        xclink = find_element(xjoint, "child")
        clink = ldict[Symbol(attribute(xclink, "link"))]

        plink = ldict[Symbol(attribute(xplink, "link"))]
        if plink.id == origin.id
            plink = origin
            joint = parse_joint(xjoint, plink, clink, T, parse_damper)
            joints = [joint; joints] # For proper parsing the first joint must be connected to the origin
        else
            joint = parse_joint(xjoint, plink, clink, T, parse_damper)
            push!(joints, joint)
        end
    end

    if floating
        originjoint = JointConstraint(Floating(origin, ldict[Symbol(floatingname)]), name=:floating_base)
        joints = [originjoint; joints] # For proper parsing the first joint must be connected to the origin
    end

    return origin, links, joints
end

# TODO This might be missing the detection of a direct loop, i.e., only two links connected by two joints
# TODO Also only works for a single loop closure in a cycle (so no ladders)
function parse_loop_joints(xloopjoints, origin, joints, ldict, T, parse_damper)
    loopjoints = JointConstraint{T}[]

    for xloopjoint in xloopjoints
        xpbody = find_element(xloopjoint, "link1")
        xcbody = find_element(xloopjoint, "link2")
        pbody = ldict[Symbol(attribute(xpbody, "link"))]
        cbody = ldict[Symbol(attribute(xcbody, "link"))]

        predlist = Tuple{Int64,Int64}[]
        jointlist = [(joints[i].id,joints[i].parent_id, joints[i].child_id) for i=1:length(joints)]
        linkid = pbody.id

        while true # create list of predecessor joints and parent links for pbody
            for (i,jointdata) in enumerate(jointlist)
                if linkid ∈ jointdata[3]
                    push!(predlist,(jointdata[1],jointdata[2]))
                    linkid = jointdata[2]
                    deleteat!(jointlist,i)
                    break
                end
            end
            if linkid == origin.id
                break
            end
        end

        jointlist = [(joints[i].id, joints[i].parent_id, joints[i].child_id) for i=1:length(joints)]
        linkid = cbody.id
        joint1id = 0
        joint2id = 0
        foundflag = false

        while true # check which predecessor link of cbody is also a predecessor link of pbody
            for (i,jointdata) in enumerate(jointlist)
                if linkid ∈ jointdata[3]
                    joint2id = jointdata[1]
                    linkid = jointdata[2]
                    deleteat!(jointlist,i)
                    break
                end
            end
            for el in predlist
                if linkid == el[2]
                    joint1id = el[1]
                    foundflag = true
                    break
                end
            end
            foundflag && break
        end

        loopjoint = parse_loop_joint(xloopjoint, pbody, cbody, T, parse_damper)
        push!(loopjoints, loopjoint)
    end

    return joints, loopjoints
end

function parse_urdf(filename, floating, ::Type{T}, parse_damper) where T
    xdoc = LightXML.parse_file(filename)
    xroot = LightXML.root(xdoc)
    @assert LightXML.name(xroot) == "robot"

    xlinks = get_elements_by_tagname(xroot, "link")
    xjoints = get_elements_by_tagname(xroot, "joint")
    xloopjoints = get_elements_by_tagname(xroot, "loop_joint")

    materialdict = parse_robotmaterials(xroot, T)
    ldict = parse_links(xlinks, materialdict, T; path_prefix=dirname(filename))
    origin, links, joints = parse_joints(xjoints, ldict, floating, T, parse_damper)

    joints, loopjoints = parse_loop_joints(xloopjoints, origin, joints, ldict, T, parse_damper)

    free(xdoc)

    return origin, links, joints, loopjoints
end


### After parsing

function set_parsed_values!(mechanism::Mechanism{T}, loopjoints) where T
    system = mechanism.system
    timestep= mechanism.timestep
    xjointlist = Dict{Int64,SVector{3,T}}() # stores id, x in world frame
    qjointlist = Dict{Int64,Quaternion{T}}() # stores id, q in world frame

    for id in root_to_leaves_ordering(mechanism, exclude_origin=true, exclude_loop_joints=true)
        node = get_node(mechanism, id)
        !(node isa JointConstraint) && continue # only for joints

        # Parent joint --> Parent body --> Child joint --> Child body
        # Child joint (joint of interest here)
        cjoint = node
        x_cjoint = cjoint.translational.vertices[1] # stored in parent_vertex
        q_cjoint = cjoint.rotational.orientation_offset # stored in orientation_offset
        # axis_pjoint = #

        # Child body (body of interest here)
        cnode = get_node(mechanism, cjoint.child_id) # x and q are stored in x2[1] and q2[1]
        shape = cnode.shape
        # x_cnode = cnode.state.x2
        # q_cnode = cnode.state.q2
        xbodylocal = cnode.state.x2
        qbodylocal = cnode.state.q2

        # Parent body
        pnode = get_node(mechanism, cjoint.parent_id, origin=true) # x and q are stored in x2[1] and q2[1]
        # x_pnode = pnode.state.x2
        # q_pnode = pnode.state.q2
        xparentbody = pnode.state.x2
        qparentbody = pnode.state.q2

        # Parent joint
        if pnode.id == 0 # parent body is origin
            # x_pjoint = SA{T}[0; 0; 0]
            # q_pjoint = one(Quaternion{T})
            # # axis_pjoint = SA{T}[1; 0; 0]
            xparentjoint = SA{T}[0; 0; 0]
            qparentjoint = one(Quaternion{T})
        else
            pjoint = get_node(mechanism, get_parent_id(mechanism, pnode.id, loopjoints))
            # x_pjoint = pjoint.translational.vertices[1] # stored in parent_vertex
            # q_pjoint = pjoint.rotational.orientation_offset # stored in orientation_offset
            # # axis_pjoint = #
            xparentjoint = xjointlist[pjoint.id] # in world frame
            qparentjoint = qjointlist[pjoint.id] # in world frame
        end

        # urdf joint's x and q in parent's (parentbody) frame
        # xjointlocal = vector_rotate(x_pjoint + vector_rotate(x_cjoint, q_pjoint) - x_pnode, inv(q_pnode))
        xjointlocal = vector_rotate(xparentjoint + vector_rotate(x_cjoint, qparentjoint) - xparentbody, inv(qparentbody))
        qjointlocal = qparentbody \ qparentjoint * q_cjoint

        # store joint's x and q in world frame
        xjoint = xparentbody + vector_rotate(xjointlocal, qparentbody)
        qjoint = qparentbody * qjointlocal
        xjointlist[cjoint.id] = xjoint
        qjointlist[cjoint.id] = qjoint

        # difference to parent body (parentbody)
        orientation_offset = qjointlocal * qbodylocal

        # actual joint properties
        parent_vertex = xjointlocal # in parent's (parentbody) frame
        child_vertex = vector_rotate(-xbodylocal, inv(qbodylocal)) # in body frame (xbodylocal and qbodylocal are both relative to the same (joint) frame -> rotationg by inv(body.q) gives body frame)
        cjoint.translational.vertices = (parent_vertex, child_vertex)

        V3 = vector_rotate(cjoint.rotational.axis_mask3', qjointlocal) # in parent's (parentbody) frame
        V1 = (svd(skew(V3)).Vt)[1:1,:]
        V2 = (svd(skew(V3)).Vt)[2:2,:]

        cjoint.rotational.axis_mask3 = SVector{3}(V3)'
        cjoint.rotational.axis_mask1 = SVector{3}(V1)'
        cjoint.rotational.axis_mask2 = SVector{3}(V2)'
        cjoint.rotational.orientation_offset = orientation_offset # in parent's (parentbody) frame

        # actual body properties
        set_maximal_configurations!(cnode) # set everything to zero
        set_maximal_configurations!(pnode, cnode; parent_vertex, child_vertex, Δq = orientation_offset)
        xbody = cnode.state.x2
        qbody = cnode.state.q2

        # shape relative
        if !(typeof(shape) <: EmptyShape)
            shape.position_offset = vector_rotate(xjoint + vector_rotate(shape.position_offset, qjoint) - xbody, inv(qbody))
            shape.orientation_offset = orientation_offset \ qjointlocal * shape.orientation_offset
        end
    end

    for (i,constraint) in enumerate(loopjoints)

        parent_id1 = constraint.parent_id
        parent_id2 = constraint.child_id
        if parent_id1 == 0 # predecessor is origin
            parentpbody = mechanism.origin

            xparentpbody = SA{T}[0; 0; 0]
            qparentpbody = one(Quaternion{T})

            xparentjoint1 = SA{T}[0; 0; 0]
            qparentjoint1 = one(Quaternion{T})
        else
            parentpbody = get_body(mechanism, parent_id1)

            grandparent_id1 = get_parent_id(mechanism, parent_id1, loopjoints)
            parentconstraint1 = get_joint(mechanism, grandparent_id1)

            xparentpbody = parentpbody.state.x2 # in world frame
            qparentpbody = parentpbody.state.q2 # in world frame

            xparentjoint1 = xjointlist[parentconstraint1.id] # in world frame
            qparentjoint1 = qjointlist[parentconstraint1.id] # in world frame
        end
        parentcbody = get_body(mechanism, parent_id2)

        grandparent_id2 = get_parent_id(mechanism, parent_id2, loopjoints)
        parentconstraint2 = get_joint(mechanism, grandparent_id2)

        xparentcbody = parentcbody.state.x2 # in world frame
        qparentcbody = parentcbody.state.q2 # in world frame

        xparentjoint2 = xjointlist[parentconstraint2.id] # in world frame
        qparentjoint2 = qjointlist[parentconstraint2.id] # in world frame


        ind1 = 1
        ind2 = ind1+1

        # urdf joint's x and q in parent's (parentbody) frame
        xjointlocal1 = vector_rotate(xparentjoint1 + vector_rotate(constraint.translational.vertices[1], qparentjoint1) - xparentpbody, inv(qparentpbody))
        xjointlocal2 = vector_rotate(xparentjoint2 + vector_rotate(constraint.translational.vertices[2], qparentjoint2) - xparentcbody, inv(qparentcbody))
        qjointlocal1 = qparentpbody \ qparentjoint1 * constraint.rotational.orientation_offset

        # difference to parent body (parentbody)
        orientation_offset1 = qjointlocal1 * qparentcbody #  qparentcbody = body in for loop above

        # actual joint properties
        parent_vertex = xjointlocal1 # in parent's (parentpbody) frame
        child_vertex = xjointlocal2 # in parent's (parentcbody) frame
        constraint.translational.vertices = (parent_vertex, child_vertex)

        V3 = vector_rotate(constraint.rotational.axis_mask3', qjointlocal1) # in parent's (parentpbody) frame
        V1 = (svd(skew(V3)).Vt)[1:1,:]
        V2 = (svd(skew(V3)).Vt)[2:2,:]

        constraint.rotational.axis_mask3 = SVector{3}(V3)'
        constraint.rotational.axis_mask1 = SVector{3}(V1)'
        constraint.rotational.axis_mask2 = SVector{3}(V2)'

        constraint.rotational.orientation_offset = orientation_offset1 # in parent's (parentpbody) frame
    end
end

function get_parent_id(mechanism, id, loopjoints)
    system = mechanism.system
    conns = connections(system, id)
    for connsid in conns
        constraint = get_joint(mechanism, connsid)
        if constraint ∉ loopjoints && id == constraint.child_id
            return connsid
        end
    end

    return nothing
end
