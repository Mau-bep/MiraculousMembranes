"""Node group, materials and modifier-input helpers for Membrane Normals.

Each arrow set is one object whose mesh is just a point per membrane vertex,
carrying the data `ops.py` measured or read from a file:

  ne_normal     unit vertex normal (the offset goes along it)
  ne_direction  unit vector the arrow points along
  ne_magnitude  value that colours (and optionally sizes) the arrow: |H| for
                mean curvature normals, |F| for forces

The `Normal Arrows` Geometry Nodes modifier turns those points into arrows, so
offset, length, thickness, density and the colour range stay live.  Instances
are realized, which keeps the per-arrow `ne_value` (magnitude mapped to 0..1)
readable by the material.
"""

import math

import bpy

ARROWS_GROUP = "Normal Arrows"
ARROWS_MOD = "Normal Arrows"
NORMAL_MAT = "Normal Arrows"
COLORED_MAT = "Magnitude Arrows"
VERSION = 2

NORMAL_ATTR = "ne_normal"
DIRECTION_ATTR = "ne_direction"
MAGNITUDE_ATTR = "ne_magnitude"
VALUE_ATTR = "ne_value"

HEAD_WIDTH = 2.4          # head radius / shaft radius
SIDES = 12                # sides of the shaft and the head

# 9 sRGB stops, sampled from matplotlib
COLORMAPS = {
    "VIRIDIS": ("Viridis", [(0.267, 0.0049, 0.3294), (0.2788, 0.1755, 0.4834), (0.2297, 0.3224, 0.5457),
                            (0.1727, 0.4488, 0.5579), (0.1276, 0.5669, 0.5506), (0.1579, 0.6838, 0.5017),
                            (0.3692, 0.7889, 0.3829), (0.6785, 0.8637, 0.1895), (0.9932, 0.9062, 0.1439)]),
    "PLASMA": ("Plasma", [(0.0504, 0.0298, 0.528), (0.2999, 0.0096, 0.6316), (0.4949, 0.012, 0.6579),
                          (0.6651, 0.1386, 0.5856), (0.7982, 0.2802, 0.4695), (0.9018, 0.4251, 0.3597),
                          (0.9734, 0.5858, 0.2515), (0.993, 0.7717, 0.1548), (0.94, 0.9752, 0.1313)]),
    "INFERNO": ("Inferno", [(0.0015, 0.0005, 0.0139), (0.1293, 0.0473, 0.2908), (0.3415, 0.0623, 0.4294),
                            (0.5409, 0.1347, 0.4151), (0.7357, 0.2159, 0.3302), (0.8943, 0.3534, 0.1936),
                            (0.9784, 0.5579, 0.0349), (0.9746, 0.7977, 0.2063), (0.9884, 0.9984, 0.6449)]),
    "TURBO": ("Turbo", [(0.19, 0.0718, 0.2322), (0.2763, 0.4212, 0.8912), (0.1584, 0.7355, 0.9231),
                        (0.1966, 0.949, 0.5947), (0.6436, 0.99, 0.2336), (0.933, 0.8124, 0.2267),
                        (0.9836, 0.4929, 0.1285), (0.8161, 0.1846, 0.0181), (0.4796, 0.0158, 0.0106)]),
    "COOLWARM": ("Cool to Warm", [(0.2298, 0.2987, 0.7537), (0.3837, 0.5102, 0.9178), (0.5543, 0.6901, 0.9955),
                                  (0.724, 0.8149, 0.9757), (0.8674, 0.8644, 0.8626), (0.9595, 0.767, 0.6741),
                                  (0.9567, 0.598, 0.4773), (0.8654, 0.3711, 0.2958), (0.7057, 0.0156, 0.1502)]),
}
DEFAULT_COLORMAP = "VIRIDIS"


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
# Finding arrow objects
# ---------------------------------------------------------------------------

def arrow_modifier(obj):
    """The Normal Arrows modifier if `obj` is an arrow set, else None."""
    if obj is None or obj.type != 'MESH':
        return None
    for mod in obj.modifiers:
        if mod.type == 'NODES' and mod.node_group is not None \
                and mod.node_group.get("na_role") == "arrows":
            return mod
    return None


def source_of(arrows):
    return arrows.get("ne_source")


def arrow_sets_of(source):
    return [o for o in bpy.data.objects
            if arrow_modifier(o) is not None and source_of(o) == source]


# ---------------------------------------------------------------------------
# Materials
# ---------------------------------------------------------------------------

def _srgb_to_linear(c):
    return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4


def _new_material(name):
    mat = bpy.data.materials.new(name)
    if not mat.node_tree:
        mat.use_nodes = True
    mat.node_tree.nodes.clear()
    return mat


def _principled(mat, x=0):
    nt = mat.node_tree
    bsdf = nt.nodes.new("ShaderNodeBsdfPrincipled")
    bsdf.location = (x, 0)
    out = nt.nodes.new("ShaderNodeOutputMaterial")
    out.location = (x + 300, 0)
    nt.links.new(bsdf.outputs["BSDF"], out.inputs["Surface"])
    return bsdf


def ensure_normal_material():
    mat = bpy.data.materials.get(NORMAL_MAT)
    if mat is not None:
        return mat
    mat = _new_material(NORMAL_MAT)
    color = (0.95, 0.32, 0.08, 1.0)
    _principled(mat).inputs["Base Color"].default_value = color
    mat.diffuse_color = color
    return mat


def set_colormap(mat, key):
    """Write the colormap stops into the colour ramp of the magnitude material."""
    name, stops = COLORMAPS[key]
    ramp = mat.node_tree.nodes["Colormap"].color_ramp
    while len(ramp.elements) > 1:
        ramp.elements.remove(ramp.elements[-1])
    ramp.interpolation = 'LINEAR'
    n = len(stops)
    for i, rgb in enumerate(stops):
        el = ramp.elements[0] if i == 0 else ramp.elements.new(i / (n - 1))
        el.position = i / (n - 1)
        el.color = (*(_srgb_to_linear(c) for c in rgb), 1.0)
    mat["na_colormap"] = key
    mat.diffuse_color = (*stops[n // 2], 1.0)


def ensure_colored_material():
    mat = bpy.data.materials.get(COLORED_MAT)
    if mat is not None:
        return mat
    mat = _new_material(COLORED_MAT)
    nt = mat.node_tree
    bsdf = _principled(mat, 400)
    attr = nt.nodes.new("ShaderNodeAttribute")
    attr.attribute_type = 'GEOMETRY'
    attr.attribute_name = VALUE_ATTR
    attr.location = (-100, 0)
    ramp = nt.nodes.new("ShaderNodeValToRGB")
    ramp.name = ramp.label = "Colormap"
    ramp.location = (150, 0)
    nt.links.new(attr.outputs["Fac"], ramp.inputs["Fac"])
    nt.links.new(ramp.outputs["Color"], bsdf.inputs["Base Color"])
    set_colormap(mat, DEFAULT_COLORMAP)
    return mat


# ---------------------------------------------------------------------------
# Small builders
# ---------------------------------------------------------------------------

def _sock(sockets, name, idx=0):
    """The idx-th enabled socket called `name` (typed nodes keep hidden twins)."""
    return [s for s in sockets if s.name == name and s.enabled][idx]


def _iface(ng, name, stype, io='INPUT', default=None, mn=None, mx=None,
           subtype=None, desc=None, parent=None):
    item = ng.interface.new_socket(name=name, in_out=io, socket_type=stype, parent=parent)
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

    def put(self, src, dst):
        """Link a socket, or set a constant."""
        if isinstance(src, bpy.types.NodeSocket):
            self.links.new(src, dst)
        else:
            dst.default_value = src

    def math(self, op, a, b=None, x=0, y=0):
        n = self.node("ShaderNodeMath", x, y, operation=op)
        self.put(a, n.inputs[0])
        if b is not None:
            self.put(b, n.inputs[1])
        return n.outputs[0]

    def switch(self, kind, cond, false, true, x=0, y=0):
        n = self.node("GeometryNodeSwitch", x, y, input_type=kind)
        self.put(cond, n.inputs["Switch"])
        self.put(false, _sock(n.inputs, "False"))
        self.put(true, _sock(n.inputs, "True"))
        return n.outputs["Output"]

    def attribute(self, name, kind, x=0, y=0):
        n = self.node("GeometryNodeInputNamedAttribute", x, y, data_type=kind)
        _sock(n.inputs, "Name").default_value = name
        return n.outputs[0]


def _clear(ng):
    ng.nodes.clear()
    ng.interface.clear()


# ---------------------------------------------------------------------------
# Normal Arrows  (modifier on every arrow set)
# ---------------------------------------------------------------------------

def ensure_arrows_group():
    ng = bpy.data.node_groups.get(ARROWS_GROUP)
    if ng is not None and ng.get("na_version") == VERSION:
        return ng
    if ng is None:
        ng = bpy.data.node_groups.new(ARROWS_GROUP, "GeometryNodeTree")
    else:
        _clear(ng)
    ng["na_version"] = VERSION
    ng["na_role"] = "arrows"

    _iface(ng, "Geometry", "NodeSocketGeometry", io='OUTPUT')
    _iface(ng, "Geometry", "NodeSocketGeometry")
    _iface(ng, "Color by Magnitude", "NodeSocketBool", default=False,
           desc="Colour the arrows by their magnitude (|H| or |F|) with the colormap. Off: one colour")
    p = ng.interface.new_panel("Arrow")
    _iface(ng, "Offset", "NodeSocketFloat", default=0.0, subtype='DISTANCE', parent=p,
           desc="Distance of the arrow base from the membrane, along the vertex normal. 0 touches the mesh")
    _iface(ng, "Length", "NodeSocketFloat", default=0.1, mn=0.0, subtype='DISTANCE', parent=p,
           desc="Length of every arrow")
    _iface(ng, "Thickness", "NodeSocketFloat", default=0.05, mn=0.0, mx=0.5, parent=p,
           desc="Shaft radius as a fraction of the length")
    _iface(ng, "Head Length", "NodeSocketFloat", default=0.3, mn=0.0, mx=1.0, parent=p,
           desc="Length of the head as a fraction of the length")
    _iface(ng, "Flip", "NodeSocketBool", default=False, parent=p,
           desc="Reverse the arrows (the offset still goes along the normal)")
    _iface(ng, "Size by Magnitude", "NodeSocketBool", default=False, parent=p,
           desc="Scale each arrow by its magnitude on the colour range: Length is the arrow at Range Max")
    _iface(ng, "Density", "NodeSocketFloat", default=1.0, mn=0.0, mx=1.0, subtype='FACTOR', parent=p,
           desc="Fraction of the vertices that get an arrow")
    _iface(ng, "Seed", "NodeSocketInt", default=0, parent=p,
           desc="Which vertices stay when Density is below 1")
    p = ng.interface.new_panel("Color")
    _iface(ng, "Range Min", "NodeSocketFloat", default=0.0, parent=p,
           desc="Magnitude that gets the first colour of the colormap")
    _iface(ng, "Range Max", "NodeSocketFloat", default=1.0, parent=p,
           desc="Magnitude that gets the last colour of the colormap")
    p = ng.interface.new_panel("Materials", default_closed=True)
    _iface(ng, "Normal Material", "NodeSocketMaterial", parent=p)
    _iface(ng, "Colored Material", "NodeSocketMaterial", parent=p)

    b = _Builder(ng)
    gi = b.node("NodeGroupInput", 0, 0)
    go = b.node("NodeGroupOutput", 16, 0)

    # --- per-point fields -------------------------------------------------
    normal = b.attribute(NORMAL_ATTR, 'FLOAT_VECTOR', 1, 6)
    along = b.attribute(DIRECTION_ATTR, 'FLOAT_VECTOR', 1, 5)
    magnitude = b.attribute(MAGNITUDE_ATTR, 'FLOAT', 1, 4)

    sign_flip = b.switch('FLOAT', gi.outputs["Flip"], 1.0, -1.0, 4, 3)
    direction = b.node("ShaderNodeVectorMath", 6, 5, operation='SCALE')
    b.link(along, direction.inputs[0])
    b.link(sign_flip, _sock(direction.inputs, "Scale"))

    offset = b.node("ShaderNodeVectorMath", 3, 6, operation='SCALE')
    b.link(normal, offset.inputs[0])
    b.link(gi.outputs["Offset"], _sock(offset.inputs, "Scale"))

    value = b.node("ShaderNodeMapRange", 3, 2, data_type='FLOAT', interpolation_type='LINEAR', clamp=True)
    b.link(magnitude, _sock(value.inputs, "Value"))
    b.link(gi.outputs["Range Min"], _sock(value.inputs, "From Min"))
    b.link(gi.outputs["Range Max"], _sock(value.inputs, "From Max"))
    result = _sock(value.outputs, "Result")

    # an arrow of zero size is not worth drawing, but keep it findable
    size = b.switch('FLOAT', gi.outputs["Size by Magnitude"], 1.0, b.math('MAXIMUM', result, 0.02, 4, 1), 5, 1)
    scale = b.node("ShaderNodeCombineXYZ", 6, 1)
    for axis in ("X", "Y", "Z"):
        b.link(size, scale.inputs[axis])

    rand = b.node("FunctionNodeRandomValue", 3, -1, data_type='FLOAT')
    _sock(rand.inputs, "Min").default_value = 0.0
    _sock(rand.inputs, "Max").default_value = 1.0
    b.link(gi.outputs["Seed"], rand.inputs["Seed"])
    drop = b.node("FunctionNodeCompare", 4, -1, data_type='FLOAT', operation='GREATER_THAN')
    b.link(_sock(rand.outputs, "Value"), _sock(drop.inputs, "A"))
    b.link(gi.outputs["Density"], _sock(drop.inputs, "B"))

    # --- points: offset, tag with the colour value, thin out ----------------
    move = b.node("GeometryNodeSetPosition", 5, 0)
    b.link(gi.outputs["Geometry"], move.inputs["Geometry"])
    b.link(offset.outputs["Vector"], move.inputs["Offset"])

    store = b.node("GeometryNodeStoreNamedAttribute", 6, 0, data_type='FLOAT', domain='POINT')
    _sock(store.inputs, "Name").default_value = VALUE_ATTR
    b.link(move.outputs["Geometry"], store.inputs["Geometry"])
    b.link(result, _sock(store.inputs, "Value"))

    thin = b.node("GeometryNodeDeleteGeometry", 7, 0, domain='POINT', mode='ALL')
    b.link(store.outputs["Geometry"], thin.inputs["Geometry"])
    b.link(drop.outputs["Result"], thin.inputs["Selection"])

    # --- the arrow: shaft + head along +Z, base at the origin ------------------
    head = b.math('MULTIPLY', gi.outputs["Length"], gi.outputs["Head Length"], 2, -4)
    shaft = b.math('SUBTRACT', gi.outputs["Length"], head, 3, -4)
    radius = b.math('MULTIPLY', gi.outputs["Length"], gi.outputs["Thickness"], 2, -5)

    def seat(mesh, base_z, x, y):
        """Move `mesh` along Z so its lowest point sits at `base_z` (primitives differ in where they are anchored)."""
        low = b.node("GeometryNodeBoundBox", x, y)
        b.link(mesh, low.inputs["Geometry"])
        z = b.node("ShaderNodeSeparateXYZ", x + 1, y)
        b.link(low.outputs["Min"], z.inputs["Vector"])
        lift = b.node("ShaderNodeCombineXYZ", x + 2, y)
        b.link(b.math('SUBTRACT', base_z, z.outputs["Z"], x + 1, y - 1), lift.inputs["Z"])
        move = b.node("GeometryNodeTransform", x + 3, y)
        b.link(mesh, move.inputs["Geometry"])
        b.link(lift.outputs[0], move.inputs["Translation"])
        return move.outputs["Geometry"]

    cyl = b.node("GeometryNodeMeshCylinder", 4, -3)
    cyl.inputs["Vertices"].default_value = SIDES
    b.link(radius, cyl.inputs["Radius"])
    b.link(shaft, cyl.inputs["Depth"])
    shaft_mesh = seat(cyl.outputs["Mesh"], 0.0, 5, -3)

    cone = b.node("GeometryNodeMeshCone", 4, -5)
    cone.inputs["Vertices"].default_value = SIDES
    _sock(cone.inputs, "Radius Top").default_value = 0.0
    b.link(b.math('MULTIPLY', radius, HEAD_WIDTH, 3, -5), _sock(cone.inputs, "Radius Bottom"))
    b.link(head, cone.inputs["Depth"])
    head_mesh = seat(cone.outputs["Mesh"], shaft, 5, -5)

    arrow = b.node("GeometryNodeJoinGeometry", 7, -4)
    b.link(shaft_mesh, arrow.inputs[0])
    b.link(head_mesh, arrow.inputs[0])

    # --- one arrow per point, +Z turned onto the direction ------------------
    align = b.node("FunctionNodeAlignRotationToVector", 8, 3, axis='Z')
    b.link(direction.outputs["Vector"], align.inputs["Vector"])

    inst = b.node("GeometryNodeInstanceOnPoints", 9, 0)
    b.link(thin.outputs["Geometry"], inst.inputs["Points"])
    b.link(arrow.outputs[0], inst.inputs["Instance"])
    b.link(align.outputs["Rotation"], inst.inputs["Rotation"])
    b.link(scale.outputs[0], inst.inputs["Scale"])

    realize = b.node("GeometryNodeRealizeInstances", 10, 0)
    b.link(inst.outputs["Instances"], realize.inputs["Geometry"])

    material = b.switch('MATERIAL', gi.outputs["Color by Magnitude"],
                        gi.outputs["Normal Material"], gi.outputs["Colored Material"], 10, -2)
    setmat = b.node("GeometryNodeSetMaterial", 11, 0)
    b.link(realize.outputs["Geometry"], setmat.inputs["Geometry"])
    b.link(material, setmat.inputs["Material"])
    b.link(setmat.outputs["Geometry"], go.inputs["Geometry"])
    return ng
