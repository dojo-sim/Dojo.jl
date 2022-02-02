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
        ixy = parse_scalar(xinertia, "ixy", T, default = "0")
        ixz = parse_scalar(xinertia, "ixz", T, default = "0")
        iyy = parse_scalar(xinertia, "iyy", T)
        iyz = parse_scalar(xinertia, "iyz", T, default = "0")
        izz = parse_scalar(xinertia, "izz", T)
        J = [ixx ixy ixz; ixy iyy iyz; ixz iyz izz]
    end

    return J
end

function parse_pose(xpose, T)
    if xpose === nothing
        x, q = zeros(T, 3), one(UnitQuaternion{T})
    else
        x = parse_vector(xpose, "xyz", T, default = "0 0 0")
        rpy = parse_vector(xpose, "rpy", T, default = "0 0 0")
        q = UnitQuaternion(RotZYX(rpy[3], rpy[2], rpy[1]))
    end

    return x, q
end

function parse_inertia(xinertial, T)
    if xinertial === nothing
        x = zeros(T, 3)
        q = one(UnitQuaternion{T})
        m = zero(T)
        J = zeros(T, 3, 3)
    else
        x, q = parse_pose(find_element(xinertial, "origin"), T)
        J = parse_inertiamatrix(find_element(xinertial, "inertia"), T)
        m = parse_scalar(find_element(xinertial, "mass"), "value", T, default = "0")
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
        colorvec = parse_vector(colornode, "rgba", T, default = string(colorvec[1]," ",colorvec[2]," ",colorvec[3]," ",colorvec[4]))

        color = RGBA(colorvec...)
    end

    return color
end

function parse_shape(xvisual, materialdict, T)
# function parse_shape(xvisual, materialdict, T, xb, qb)

    if xvisual === nothing
        shape = nothing
    else
        xgeometry = find_element(xvisual, "geometry")
        @assert xgeometry !== nothing

        color = parse_xmaterial(find_element(xvisual, "material"), materialdict, T)
        x, q = parse_pose(find_element(xvisual, "origin"), T)
        # x = vrotate(x - xb, inv(qb))
        # q = inv(qb) * q

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
            shape = get_shape(shapenode, x, q, color, T)
        end
    end
    return shape
end

function get_shape(shapenode, x, q, color, T)
    if name(shapenode) == "box"
        xyz = parse_vector(shapenode, "size", T, default = "1 1 1")
        shape = Box(xyz..., zero(T), color = color, xoffset = x, qoffset = q)
    elseif name(shapenode) == "cylinder"
        r = parse_scalar(shapenode, "radius", T, default = "0.5")
        l = parse_scalar(shapenode, "length", T, default = "1")
        shape = Cylinder(r, l, zero(T), color = color, xoffset = x, qoffset = q)
    elseif name(shapenode) == "pyramid"
        w = parse_scalar(shapenode, "width", T, default = "1")
        h = parse_scalar(shapenode, "height", T, default = "1")
        shape = Cylinder(w, h, zero(T), color = color, xoffset = x, qoffset = q)
    elseif name(shapenode) == "sphere"
        r = parse_scalar(shapenode, "radius", T, default = "0.5")
        shape = Sphere(r, zero(T), color = color, xoffset = x, qoffset = q)
    elseif name(shapenode) == "mesh"
        path = attribute(shapenode, "filename")
        scale = parse_vector(shapenode, "scale", T, default = "1 1 1")
        shape = Mesh(path, zero(T), zeros(T, 3, 3), scale=scale, color = color, xoffset = x, qoffset = q)
    elseif name(shapenode) == "capsule"
        r = parse_scalar(shapenode, "radius", T, default = "0.5")
        l = parse_scalar(shapenode, "length", T, default = "1")
        shape = Capsule(r, l, zero(T), color = color, xoffset = x, qoffset = q)
    else
        @info "Unknown geometry."
        shape = nothing
    end
end

function parse_link(xlink, materialdict, T)
    x, q, m, J = parse_inertia(find_element(xlink, "inertial"), T)
    xvisuals = get_elements_by_tagname(xlink, "visual")
    shapes = [parse_shape(xvisual, materialdict, T) for xvisual in xvisuals]
    # shapes = [parse_shape(xvisual, materialdict, T, x, q) for xvisual in xvisuals]
    if length(shapes) == 0
        shape = nothing
    elseif length(shapes) > 1
        s = [s.shape for s in shapes]
        shape = Shapes(s, 0.0, diagm([0.0; 0.0; 0.0]))
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

    link.state.x2[1] = x
    link.state.q2[1] = q


    return link
end

function parse_links(xlinks, materialdict, T)
    ldict = Dict{Symbol,Body{T}}()

    for xlink in xlinks
        link = parse_link(xlink, materialdict, T)
        ldict[link.name] = link
    end

    return ldict
end

function joint_selector(jointtype, body1, body2, T;
        axis = SA{T}[1;0;0], p1 = szeros(T,3), p2 = szeros(T,3), qoffset = one(UnitQuaternion{T}), name = Symbol("joint_" * randstring(4)))

    # TODO limits for revolute joint?
    if jointtype == "revolute" || jointtype == "continuous"
        joint = JointConstraint(Revolute(body1, body2, axis; p1=p1, p2=p2, qoffset = qoffset), name=name)
    elseif jointtype == "prismatic"
        joint = JointConstraint(Prismatic(body1, body2, axis; p1=p1, p2=p2, qoffset = qoffset), name=name)
    elseif jointtype == "planar"
        joint = JointConstraint(Planar(body1, body2, axis; p1=p1, p2=p2, qoffset = qoffset), name=name)
    elseif jointtype == "planarfree"
        joint = JointConstraint(PlanarFree(body1, body2, axis; p1=p1, p2=p2), name=name)
    elseif jointtype == "fixed"
        joint = JointConstraint(Fixed(body1, body2; p1=p1, p2=p2, qoffset = qoffset), name=name)
    elseif jointtype == "floating"
        joint = JointConstraint(Floating(body1, body2), name=name)
    elseif jointtype == "orbital"
        joint = JointConstraint(Orbital(body1, body2, axis; p1=p1, p2=p2, qoffset = qoffset), name=name)
    elseif jointtype == "ball"
        joint = JointConstraint(Spherical(body1, body2; p1=p1, p2=p2, qoffset = qoffset), name=name)
    elseif jointtype == "fixedorientation"
        joint = JointConstraint(FixedOrientation(body1, body2; qoffset = qoffset), name=name)
    elseif jointtype == "cylindrical"
        joint = JointConstraint(Cylindrical(body1, body2, axis; p1=p1, p2=p2, qoffset = qoffset), name=name)
    elseif jointtype == "cylindricalfree"
        joint = JointConstraint(CylindricalFree(body1, body2, axis; p1=p1, p2=p2), name=name)
    elseif jointtype == "planaraxis"
        joint = JointConstraint(PlanarAxis(body1, body2, axis; p1=p1, p2=p2, qoffset = qoffset), name=name)
    else
        @error "Unknown joint type"
    end

    return joint
end

function parse_joint(xjoint, plink, clink, T)
    jointtype = attribute(xjoint, "type")
    x, q = parse_pose(find_element(xjoint, "origin"), T)
    axis = parse_vector(find_element(xjoint, "axis"), "xyz", T, default = "1 0 0")
    p1 = x
    name = Symbol(attribute(xjoint, "name"))

    return joint_selector(jointtype, plink, clink, T, axis = axis, p1 = p1, qoffset = q, name = name)
end

function parse_loop_joint(xjoint, body1, body2, T)
    find_element(xjoint, "body1")
    find_element(xjoint, "body2")

    jointtype = attribute(xjoint, "type")
    axis = parse_vector(find_element(xjoint, "axis"), "xyz", T, default = "1 0 0")
    x1, q1 = parse_pose(find_element(xjoint, "body1"), T)
    x2, _ = parse_pose(find_element(xjoint, "body2"), T) # The orientation q2 of the second body is ignored because it is determined by the mechanism's structure
    p1 = x1
    p2 = x2
    name = Symbol(attribute(xjoint, "name"))

    return joint_selector(jointtype, body1, body2, T, axis = axis, p1 = p1, p2 = p2, qoffset = q1, name = name)
end

function parse_joints(xjoints, ldict, floating, T)
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
            joint = parse_joint(xjoint, plink, clink, T)
            joints = [joint; joints] # For proper parsing the first joint must be connected to the origin
        else
            joint = parse_joint(xjoint, plink, clink, T)
            push!(joints, joint)
        end
    end

    if floating
        originjoint = JointConstraint(Floating(origin, ldict[Symbol(floatingname)]), name=:auto_generated_floating_joint)
        joints = [originjoint; joints] # For proper parsing the first joint must be connected to the origin
    end

    return origin, links, joints
end

# TODO This might be missing the detection of a direct loop, i.e. only two links connected by two joints
# TODO Also only works for a single loop closure in a cycle (so no ladders)
function parse_loop_joints(xloopjoints, origin, joints, ldict, T)
    loopjoints = JointConstraint{T}[]


    for xloopjoint in xloopjoints
        xbody1 = find_element(xloopjoint, "body1")
        xbody2 = find_element(xloopjoint, "body2")
        body1 = ldict[attribute(xbody1, "link")]
        body2 = ldict[attribute(xbody2, "link")]

        predlist = Tuple{Int64,Int64}[]
        jointlist = [(joints[i].id,joints[i].parent_id, joints[i].child_id) for i=1:length(joints)]
        linkid = body1.id

        while true # create list of predecessor joints and parent links for body1
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
        linkid = body2.id
        joint1id = 0
        joint2id = 0
        foundflag = false

        while true # check which predecessor link of body2 is also a predecessor link of body1
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

        # Find and remove joints to combine them
        joint1 = 0
        joint2 = 0
        for (i,joint) in enumerate(joints)
            if joint.id == joint1id
                joint1 = joint
                deleteat!(joints,i)
                break
            end
        end
        # if joint1id == joint2id # already combined joint
        #     joint2 = joint
        for (i,joint) in enumerate(joints)
            if joint.id == joint2id
                joint2 = joint
                deleteat!(joints,i)
                break
            end
        end

        joint = cat(joint1,joint2)
        push!(joints,joint)
        loopjoint = parse_loop_joint(xloopjoint, body1, body2, T)
        push!(loopjoints, loopjoint)
    end

    return joints, loopjoints
end

function parse_urdf(filename, floating, ::Type{T}) where T
    xdoc = LightXML.parse_file(filename)
    xroot = LightXML.root(xdoc)
    @assert LightXML.name(xroot) == "robot"

    xlinks = get_elements_by_tagname(xroot, "link")
    xjoints = get_elements_by_tagname(xroot, "joint")
    xloopjoints = get_elements_by_tagname(xroot, "loop_joint")

    materialdict = parse_robotmaterials(xroot, T)
    ldict = parse_links(xlinks, materialdict, T)
    origin, links, joints = parse_joints(xjoints, ldict, floating, T)

    joints, loopjoints = parse_loop_joints(xloopjoints, origin, joints, ldict, T)

    free(xdoc)

    return origin, links, joints, loopjoints
end


### After parsing

function set_parsed_values!(mechanism::Mechanism{T}, loopjoints) where T
    system = mechanism.system
    xjointlist = Dict{Int64,SVector{3,T}}() # stores id, x in world frame
    qjointlist = Dict{Int64,UnitQuaternion{T}}() # stores id, q in world frame

    for id in reverse(system.dfs_list) # from root to leaves
        node = get_node(mechanism, id)
        !(node isa Body) && continue # only for bodies

        body = node
        xbodylocal = body.state.x2[1]
        qbodylocal = body.state.q2[1]
        shape = body.shape

        parent_id = get_parent_id(mechanism, id, loopjoints)
        constraint = get_joint_constraint(mechanism, parent_id)

        grandparent_id = constraint.parent_id
        if grandparent_id == 0 # predecessor is origin
            parentbody = mechanism.origin

            xparentbody = SA{T}[0; 0; 0]
            qparentbody = one(UnitQuaternion{T})

            xparentjoint = SA{T}[0; 0; 0]
            qparentjoint = one(UnitQuaternion{T})
        else
            parentbody = get_body(mechanism, grandparent_id)

            grandgrandparent_id = get_parent_id(mechanism, grandparent_id, loopjoints)
            parentconstraint = get_joint_constraint(mechanism, grandgrandparent_id)

            xparentbody = parentbody.state.x2[1] # in world frame
            qparentbody = parentbody.state.q2[1] # in world frame

            xparentjoint = xjointlist[parentconstraint.id] # in world frame
            qparentjoint = qjointlist[parentconstraint.id] # in world frame
        end

        ind1 = findfirst(x -> x == id, constraint.child_id)
        ind2 = ind1+1

        # urdf joint's x and q in parent's (parentbody) frame
        xjointlocal = vrotate(xparentjoint + vrotate(constraint.constraints[ind1].vertices[1], qparentjoint) - xparentbody, inv(qparentbody))
        qjointlocal = qparentbody \ qparentjoint * constraint.constraints[ind2].qoffset

        # store joint's x and q in world frame
        xjoint = xparentbody + vrotate(xjointlocal, qparentbody)
        qjoint = qparentbody * qjointlocal
        xjointlist[constraint.id] = xjoint
        qjointlist[constraint.id] = qjoint

        # difference to parent body (parentbody)
        qoffset = qjointlocal * qbodylocal

        # actual joint properties
        p1 = xjointlocal # in parent's (parentbody) frame
        p2 = vrotate(-xbodylocal, inv(qbodylocal)) # in body frame (xbodylocal and qbodylocal are both relative to the same (joint) frame -> rotationg by inv(body.q) gives body frame)
        constraint.constraints[ind1].vertices = (p1, p2)

        V3 = vrotate(constraint.constraints[ind2].V3', qjointlocal) # in parent's (parentbody) frame
        V12 = (svd(skew(V3)).Vt)[1:2,:]
        constraint.constraints[ind2].V3 = V3'
        constraint.constraints[ind2].V12 = V12
        constraint.constraints[ind2].qoffset = qoffset # in parent's (parentbody) frame

        # actual body properties
        set_position!(body) # set everything to zero
        set_position!(parentbody, body, p1 = p1, p2 = p2, Δq = qoffset)
        xbody = body.state.x2[1]
        qbody = body.state.q2[1]

        # shape relative
        if !(typeof(shape) <: EmptyShape)
            shape.xoffset = vrotate(xjoint + vrotate(shape.xoffset, qjoint) - xbody, inv(qbody))
            shape.qoffset = qoffset \ qjointlocal * shape.qoffset
        end
    end
    for (i,constraint) in enumerate(loopjoints)

        parent_id1 = constraint.parent_id
        parent_id2 = constraint.child_id
        if parent_id1 == 0 # predecessor is origin
            parentbody1 = mechanism.origin

            xparentbody1 = SA{T}[0; 0; 0]
            qparentbody1 = one(UnitQuaternion{T})

            xparentjoint1 = SA{T}[0; 0; 0]
            qparentjoint1 = one(UnitQuaternion{T})
        else
            parentbody1 = get_body(mechanism, parent_id1)

            grandparent_id1 = get_parent_id(mechanism, parent_id1, loopjoints)
            parentconstraint1 = get_joint_constraint(mechanism, grandparent_id1)

            xparentbody1 = parentbody1.state.x2[1] # in world frame
            qparentbody1 = parentbody1.state.q2[1] # in world frame

            xparentjoint1 = xjointlist[parentconstraint1.id] # in world frame
            qparentjoint1 = qjointlist[parentconstraint1.id] # in world frame
        end
        parentbody2 = get_body(mechanism, parent_id2)

        grandparent_id2 = get_parent_id(mechanism, parent_id2, loopjoints)
        parentconstraint2 = get_joint_constraint(mechanism, grandparent_id2)

        xparentbody2 = parentbody2.state.x2[1] # in world frame
        qparentbody2 = parentbody2.state.q2[1] # in world frame

        xparentjoint2 = xjointlist[parentconstraint2.id] # in world frame
        qparentjoint2 = qjointlist[parentconstraint2.id] # in world frame


        ind1 = 1
        ind2 = ind1+1

        # urdf joint's x and q in parent's (parentbody) frame
        xjointlocal1 = vrotate(xparentjoint1 + vrotate(constraint.constraints[ind1].vertices[1], qparentjoint1) - xparentbody1, inv(qparentbody1))
        xjointlocal2 = vrotate(xparentjoint2 + vrotate(constraint.constraints[ind1].vertices[2], qparentjoint2) - xparentbody2, inv(qparentbody2))
        qjointlocal1 = qparentbody1 \ qparentjoint1 * constraint.constraints[ind2].qoffset

        # difference to parent body (parentbody)
        qoffset1 = qjointlocal1 * qparentbody2 #  qparentbody2 = body in for loop above

        # actual joint properties
        p1 = xjointlocal1 # in parent's (parentbody1) frame
        p2 = xjointlocal2 # in parent's (parentbody2) frame
        constraint.constraints[ind1].vertices = (p1, p2)

        V3 = vrotate(constraint.constraints[ind2].V3', qjointlocal1) # in parent's (parentbody1) frame
        V12 = (svd(skew(V3)).Vt)[1:2,:]
        constraint.constraints[ind2].V3 = V3'
        constraint.constraints[ind2].V12 = V12
        constraint.constraints[ind2].qoffset = qoffset1 # in parent's (parentbody1) frame
    end
end

function get_parent_id(mechanism, id, loopjoints)
    system = mechanism.system
    conns = connections(system, id)
    for connsid in conns
        constraint = get_joint_constraint(mechanism, connsid)
        if constraint ∉ loopjoints && id == constraint.child_id
            return connsid
        end
    end

    return nothing
end
