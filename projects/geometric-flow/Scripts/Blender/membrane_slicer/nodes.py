"""Node groups and modifier-input helpers for Membrane Slicer.

Two Geometry Nodes groups do all the work, so there are no frame handlers and
nothing touches the mesh data that the OBJ-sequence addon rewrites each frame:

  MS Cutter  lives on the slicer object.  It builds the closed wedge solid
             (apex line = local Z, opening centred on local +X) and stores the
             slicer settings (enabled, frame range, ...) as point attributes on
             that solid, so one set of controls drives every target.
  MS Cut     is the modifier on each membrane.  It reads the slicer object's
             evaluated geometry (Object Info, relative transform), reads the
             settings back from it, and subtracts the wedge with an exact mesh
             boolean.  The cap faces the boolean creates on the cutter walls are
             tagged and deleted, so the membrane stays an open surface.
"""

import math

import bpy

CUTTER_GROUP = "MS Cutter"
CUT_GROUP = "MS Cut"
CUTTER_MOD = "Slicer Cutter"
CUT_MOD = "Slicer Cut"
VERSION = 1

CAP_ATTR = "ms_cap"
# modifier-input name on the cutter -> point attribute that carries it to targets
SETTINGS_ATTR = {
    "Enabled": "ms_enabled",
    "Limit to Frames": "ms_limit",
    "Frame Start": "ms_start",
    "Frame End": "ms_end",
    "Keep Cut Faces": "ms_keepcaps",
    "Self Intersections": "ms_self",
    "Hole Tolerant": "ms_hole",
}


# ---------------------------------------------------------------------------
# Modifier inputs (Blender 5.x exposes them as modifier.properties.inputs.<id>.value,
# earlier versions as modifier["<id>"])
# ---------------------------------------------------------------------------

def input_identifier(mod, name):
    ng = mod.node_group
    if ng is None:
        return None
    for item in ng.interface.items_tree:
        if item.item_type == 'SOCKET' and item.in_out == 'INPUT' and item.name == name:
            return item.identifier
    return None


def get_input(mod, name):
    ident = input_identifier(mod, name)
    if ident is None:
        return None
    props = getattr(mod, "properties", None)
    if props is not None and hasattr(props, "inputs"):
        return getattr(props.inputs, ident).value
    return mod[ident]


def set_input(mod, name, value):
    ident = input_identifier(mod, name)
    if ident is None:
        return False
    props = getattr(mod, "properties", None)
    if props is not None and hasattr(props, "inputs"):
        getattr(props.inputs, ident).value = value
    else:
        mod[ident] = value
    mod.id_data.update_tag()
    return True


def input_data_path(mod, name):
    """Animation data path (from the object) of a modifier input."""
    ident = input_identifier(mod, name)
    if ident is None:
        return None
    props = getattr(mod, "properties", None)
    if props is not None and hasattr(props, "inputs"):
        return f'modifiers["{mod.name}"].properties.inputs.{ident}.value'
    return f'modifiers["{mod.name}"]["{ident}"]'


def draw_input(layout, mod, name, text=None, **kw):
    ident = input_identifier(mod, name)
    if ident is None:
        return
    props = getattr(mod, "properties", None)
    if props is not None and hasattr(props, "inputs"):
        layout.prop(getattr(props.inputs, ident), "value", text=text or name, **kw)
    else:
        layout.prop(mod, f'["{ident}"]', text=text or name, **kw)


# ---------------------------------------------------------------------------
# Finding slicers and their targets
# ---------------------------------------------------------------------------

def _role_modifier(obj, role):
    if obj is None or obj.type != 'MESH':
        return None
    for mod in obj.modifiers:
        if mod.type == 'NODES' and mod.node_group is not None \
                and mod.node_group.get("ms_role") == role:
            return mod
    return None


def cutter_modifier(obj):
    """The MS Cutter modifier if `obj` is a slicer, else None."""
    return _role_modifier(obj, "cutter")


def cut_modifier(obj):
    """The MS Cut modifier if `obj` is a cut membrane, else None."""
    return _role_modifier(obj, "cut")


def slicer_of(obj):
    """The slicer object that cuts `obj`, or `obj` itself if it is a slicer."""
    if cutter_modifier(obj) is not None:
        return obj
    mod = cut_modifier(obj)
    return get_input(mod, "Slicer") if mod is not None else None


def targets_of(slicer):
    return [o for o in bpy.data.objects
            if cut_modifier(o) is not None and get_input(cut_modifier(o), "Slicer") == slicer]


# ---------------------------------------------------------------------------
# Small builders
# ---------------------------------------------------------------------------

def _sock(sockets, name, idx=0):
    """The idx-th enabled socket called `name` (typed nodes keep hidden twins)."""
    return [s for s in sockets if s.name == name and s.enabled][idx]


def _iface(ng, name, stype, io='INPUT', default=None, mn=None, mx=None,
           subtype=None, desc=None):
    item = ng.interface.new_socket(name=name, in_out=io, socket_type=stype)
    if default is not None:
        item.default_value = default
    if mn is not None:
        item.min_value = mn
    if mx is not None:
        item.max_value = mx
    if subtype:
        item.subtype = subtype
    if desc:
        item.description = desc
    return item


class _Builder:
    """Tiny helper: place nodes on a grid and wire them by name."""

    def __init__(self, ng):
        self.ng = ng
        self.nodes = ng.nodes
        self.links = ng.links

    def node(self, bl_idname, x, y, label=None, **props):
        n = self.nodes.new(bl_idname)
        n.location = (x * 220, y * 160)
        if label:
            n.label = label
        for k, v in props.items():
            setattr(n, k, v)
        return n

    def link(self, out_sock, in_sock):
        self.links.new(out_sock, in_sock)

    def math(self, op, a, b=None, x=0, y=0):
        n = self.node("ShaderNodeMath", x, y, operation=op)
        self.link_or_set(a, n.inputs[0])
        if b is not None:
            self.link_or_set(b, n.inputs[1])
        return n.outputs[0]

    def link_or_set(self, src, dst):
        if isinstance(src, bpy.types.NodeSocket):
            self.links.new(src, dst)
        else:
            dst.default_value = src


def _clear(ng):
    ng.nodes.clear()
    ng.interface.clear()


# ---------------------------------------------------------------------------
# MS Cutter  (on the slicer object)
# ---------------------------------------------------------------------------

def ensure_cutter_group():
    ng = bpy.data.node_groups.get(CUTTER_GROUP)
    if ng is not None and ng.get("ms_version") == VERSION:
        return ng
    if ng is None:
        ng = bpy.data.node_groups.new(CUTTER_GROUP, "GeometryNodeTree")
    else:
        _clear(ng)
    ng["ms_version"] = VERSION
    ng["ms_role"] = "cutter"

    _iface(ng, "Geometry", "NodeSocketGeometry", io='OUTPUT')
    _iface(ng, "Angle", "NodeSocketFloat", default=math.radians(90),
           mn=math.radians(1), mx=math.radians(359), subtype='ANGLE',
           desc="Opening angle of the removed slice. 180 degrees is a plane cut")
    _iface(ng, "Radius", "NodeSocketFloat", default=10.0, mn=1e-4,
           desc="How far the slice reaches from its apex line")
    _iface(ng, "Height", "NodeSocketFloat", default=10.0, mn=1e-4,
           desc="Half-height of the slice along the apex line")
    _iface(ng, "Enabled", "NodeSocketBool", default=True,
           desc="Switch the cut off without removing it")
    _iface(ng, "Limit to Frames", "NodeSocketBool", default=False,
           desc="Only cut between Frame Start and Frame End")
    _iface(ng, "Frame Start", "NodeSocketInt", default=1)
    _iface(ng, "Frame End", "NodeSocketInt", default=250)
    _iface(ng, "Keep Cut Faces", "NodeSocketBool", default=False,
           desc="Keep the faces the boolean creates on the cutter walls")
    _iface(ng, "Self Intersections", "NodeSocketBool", default=False,
           desc="Exact boolean handles membranes that intersect themselves (slower)")
    _iface(ng, "Hole Tolerant", "NodeSocketBool", default=True,
           desc="Exact boolean copes with open membranes")

    b = _Builder(ng)
    gi = b.node("NodeGroupInput", 0, 0)
    go = b.node("NodeGroupOutput", 11, 0)

    # pie outline: arc from -angle/2 to +angle/2, closed through the centre
    arc = b.node("GeometryNodeCurveArc", 3, 2, mode='RADIUS')
    _sock(arc.inputs, "Resolution").default_value = 32
    _sock(arc.inputs, "Connect Center").default_value = True
    b.link(b.math('MULTIPLY', gi.outputs["Angle"], -0.5, 2, 3), _sock(arc.inputs, "Start Angle"))
    b.link(gi.outputs["Angle"], _sock(arc.inputs, "Sweep Angle"))
    b.link(gi.outputs["Radius"], _sock(arc.inputs, "Radius"))

    fill = b.node("GeometryNodeFillCurve", 4, 2)
    b.link(arc.outputs["Curve"], fill.inputs["Curve"])

    # extrude gives the sides and top, the flipped fill is the bottom
    ext = b.node("GeometryNodeExtrudeMesh", 5, 3, mode='FACES')
    ext.inputs["Offset"].default_value = (0.0, 0.0, 1.0)
    ext.inputs["Individual"].default_value = False
    b.link(fill.outputs["Mesh"], ext.inputs["Mesh"])
    b.link(b.math('MULTIPLY', gi.outputs["Height"], 2.0, 4, 4), ext.inputs["Offset Scale"])

    flip = b.node("GeometryNodeFlipFaces", 5, 1)
    b.link(fill.outputs["Mesh"], flip.inputs["Mesh"])

    join = b.node("GeometryNodeJoinGeometry", 6, 2)
    b.link(ext.outputs["Mesh"], join.inputs[0])
    b.link(flip.outputs["Mesh"], join.inputs[0])

    merge = b.node("GeometryNodeMergeByDistance", 7, 2)
    merge.inputs["Distance"].default_value = 1e-6
    b.link(join.outputs[0], merge.inputs["Geometry"])

    lift = b.node("ShaderNodeCombineXYZ", 7, 4)
    b.link(b.math('MULTIPLY', gi.outputs["Height"], -1.0, 6, 4), lift.inputs["Z"])
    move = b.node("GeometryNodeTransform", 8, 2)
    b.link(merge.outputs["Geometry"], move.inputs["Geometry"])
    b.link(lift.outputs[0], move.inputs["Translation"])

    # carry the settings to the targets as float point attributes
    prev = move.outputs["Geometry"]
    for i, (name, attr) in enumerate(SETTINGS_ATTR.items()):
        store = b.node("GeometryNodeStoreNamedAttribute", 9, 2 - i * 0.8,
                       data_type='FLOAT', domain='POINT')
        _sock(store.inputs, "Name").default_value = attr
        b.link(prev, store.inputs["Geometry"])
        b.link(gi.outputs[name], _sock(store.inputs, "Value"))
        prev = store.outputs["Geometry"]
    b.link(prev, go.inputs["Geometry"])
    return ng


# ---------------------------------------------------------------------------
# MS Cut  (modifier on each membrane)
# ---------------------------------------------------------------------------

def ensure_cut_group():
    ng = bpy.data.node_groups.get(CUT_GROUP)
    if ng is not None and ng.get("ms_version") == VERSION:
        return ng
    if ng is None:
        ng = bpy.data.node_groups.new(CUT_GROUP, "GeometryNodeTree")
    else:
        _clear(ng)
    ng["ms_version"] = VERSION
    ng["ms_role"] = "cut"

    _iface(ng, "Geometry", "NodeSocketGeometry", io='OUTPUT')
    _iface(ng, "Geometry", "NodeSocketGeometry", io='INPUT')
    _iface(ng, "Slicer", "NodeSocketObject",
           desc="Slicer object that defines the cut")

    b = _Builder(ng)
    gi = b.node("NodeGroupInput", 0, 0)
    go = b.node("NodeGroupOutput", 13, 0)

    info = b.node("GeometryNodeObjectInfo", 1, 2, transform_space='RELATIVE')
    b.link(gi.outputs["Slicer"], info.inputs["Object"])
    cutter = info.outputs["Geometry"]

    def setting(name, y):
        na = b.node("GeometryNodeInputNamedAttribute", 2, y, data_type='FLOAT')
        _sock(na.inputs, "Name").default_value = SETTINGS_ATTR[name]
        st = b.node("GeometryNodeAttributeStatistic", 3, y, data_type='FLOAT', domain='POINT')
        b.link(cutter, st.inputs["Geometry"])
        b.link(na.outputs[0], _sock(st.inputs, "Attribute"))
        return st.outputs["Max"]

    def flag(name, y):
        cmp = b.node("FunctionNodeCompare", 4, y, data_type='FLOAT', operation='GREATER_THAN')
        b.link(setting(name, y), _sock(cmp.inputs, "A"))
        _sock(cmp.inputs, "B").default_value = 0.5
        return cmp.outputs["Result"]

    enabled = flag("Enabled", 4)
    limit = flag("Limit to Frames", 3)
    keep = flag("Keep Cut Faces", 2)
    self_x = flag("Self Intersections", 1)
    hole = flag("Hole Tolerant", 0)

    frame = b.node("GeometryNodeInputSceneTime", 4, -1).outputs["Frame"]
    ge = b.node("FunctionNodeCompare", 5, -1.5, data_type='FLOAT', operation='GREATER_EQUAL')
    le = b.node("FunctionNodeCompare", 5, -2.5, data_type='FLOAT', operation='LESS_EQUAL')
    b.link(frame, _sock(ge.inputs, "A"))
    b.link(setting("Frame Start", -1.5), _sock(ge.inputs, "B"))
    b.link(frame, _sock(le.inputs, "A"))
    b.link(setting("Frame End", -2.5), _sock(le.inputs, "B"))

    def bool_math(op, a, c=None, x=6, y=0):
        n = b.node("FunctionNodeBooleanMath", x, y, operation=op)
        b.link(a, n.inputs[0])
        if c is not None:
            b.link(c, n.inputs[1])
        return n.outputs[0]

    in_range = bool_math('AND', ge.outputs["Result"], le.outputs["Result"], 6, -2)
    free_run = bool_math('NOT', limit, None, 6, 3)
    window = bool_math('OR', free_run, in_range, 7, 1)
    active = bool_math('AND', enabled, window, 8, 2)

    # tag the cutter faces; the boolean carries the tag onto the cap faces
    tag = b.node("GeometryNodeStoreNamedAttribute", 5, 4, data_type='BOOLEAN', domain='FACE')
    _sock(tag.inputs, "Name").default_value = CAP_ATTR
    _sock(tag.inputs, "Value").default_value = True
    b.link(cutter, tag.inputs["Geometry"])

    boolean = b.node("GeometryNodeMeshBoolean", 8, 0, operation='DIFFERENCE', solver='EXACT')
    b.link(gi.outputs["Geometry"], boolean.inputs["Mesh 1"])
    b.link(tag.outputs["Geometry"], boolean.inputs["Mesh 2"])
    b.link(self_x, boolean.inputs["Self Intersection"])
    b.link(hole, boolean.inputs["Hole Tolerant"])

    cap = b.node("GeometryNodeInputNamedAttribute", 9, -1, data_type='BOOLEAN')
    _sock(cap.inputs, "Name").default_value = CAP_ATTR
    drop = bool_math('AND', cap.outputs[0], bool_math('NOT', keep, None, 9, -2), 10, -1)
    delete = b.node("GeometryNodeDeleteGeometry", 10, 0, domain='FACE', mode='ALL')
    b.link(boolean.outputs["Mesh"], delete.inputs["Geometry"])
    b.link(drop, delete.inputs["Selection"])

    # the boolean may merge clipped pieces into n-gons; keep the output all triangles
    tris = b.node("GeometryNodeTriangulate", 11, 2)
    b.link(delete.outputs["Geometry"], tris.inputs["Mesh"])

    clean = b.node("GeometryNodeRemoveAttribute", 11, 1)
    clean.inputs["Name"].default_value = CAP_ATTR
    b.link(tris.outputs["Mesh"], clean.inputs["Geometry"])

    switch = b.node("GeometryNodeSwitch", 12, 0, input_type='GEOMETRY')
    b.link(active, switch.inputs["Switch"])
    b.link(gi.outputs["Geometry"], switch.inputs["False"])
    b.link(clean.outputs["Geometry"], switch.inputs["True"])
    b.link(switch.outputs["Output"], go.inputs["Geometry"])
    return ng
