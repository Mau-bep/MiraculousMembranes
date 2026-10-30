"""Node group, modifier-input helpers and the settings applier for RenderReady.

For a mesh object (the source) RenderReady makes two more objects that SHARE the source's mesh data and
keep their own modifiers:

  <name> wire    a Wireframe modifier after whatever modifiers the source has
  <name> points  a `RenderReady Points` Geometry Nodes modifier (Mesh to Points + Set Material)

Because the mesh data is shared and nothing is applied, they follow everything that changes the source's mesh:
edit mode, shape keys, modifiers copied from the source, and the frame changes of an OBJ sequence (the sequence
add-on rewrites the one mesh datablock in place).
"""

import bpy

POINTS_GROUP = "RenderReady Points"
POINTS_MOD = "RenderReady Points"
WIRE_MOD = "RenderReady Wire"
VERSION = 1


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


# ---------------------------------------------------------------------------
# The points group
# ---------------------------------------------------------------------------

def ensure_points_group():
    for ng in bpy.data.node_groups:
        if ng.get("rr_role") == "points" and ng.get("rr_version") == VERSION:
            return ng
    ng = bpy.data.node_groups.new(POINTS_GROUP, "GeometryNodeTree")
    ng["rr_role"] = "points"
    ng["rr_version"] = VERSION

    ng.interface.new_socket(name="Geometry", in_out='OUTPUT', socket_type="NodeSocketGeometry")
    ng.interface.new_socket(name="Geometry", in_out='INPUT', socket_type="NodeSocketGeometry")
    radius = ng.interface.new_socket(name="Radius", in_out='INPUT', socket_type="NodeSocketFloat")
    radius.subtype = 'DISTANCE'
    radius.default_value = 0.02
    radius.min_value = 0.0
    radius.description = "Radius of the points"
    ng.interface.new_socket(name="Material", in_out='INPUT', socket_type="NodeSocketMaterial")

    gi = ng.nodes.new("NodeGroupInput")
    gi.location = (0, 0)
    go = ng.nodes.new("NodeGroupOutput")
    go.location = (700, 0)
    points = ng.nodes.new("GeometryNodeMeshToPoints")
    points.mode = 'VERTICES'
    points.location = (230, 0)
    material = ng.nodes.new("GeometryNodeSetMaterial")
    material.location = (460, 0)
    ng.links.new(gi.outputs["Geometry"], points.inputs["Mesh"])
    ng.links.new(gi.outputs["Radius"], points.inputs["Radius"])
    ng.links.new(points.outputs["Points"], material.inputs["Geometry"])
    ng.links.new(gi.outputs["Material"], material.inputs["Material"])
    ng.links.new(material.outputs["Geometry"], go.inputs["Geometry"])
    return ng


# ---------------------------------------------------------------------------
# Finding the parts
# ---------------------------------------------------------------------------

def parts_of(src):
    """(wire, points) objects that belong to `src`; either is None when missing. Checked both ways, because
    duplicating a source copies its custom properties but not the ownership."""
    parts = []
    for key in ("rr_wire", "rr_points"):
        part = src.get(key)
        parts.append(part if part is not None and part.get("rr_source") == src else None)
    return tuple(parts)


def source_of(obj):
    """The source object `obj` is, or the source whose wire or points copy it is; None for any other object."""
    if obj is None:
        return None
    if any(parts_of(obj)):
        return obj
    src = obj.get("rr_source")
    if src is not None and obj in parts_of(src):
        return src
    return None


def wire_modifier(wire):
    return wire.modifiers.get(WIRE_MOD) if wire is not None else None


def points_modifier(points):
    return points.modifiers.get(POINTS_MOD) if points is not None else None


# ---------------------------------------------------------------------------
# Applying the settings
# ---------------------------------------------------------------------------

def set_materials(obj, material):
    """Give every material slot of `obj` the material, on the object only: the mesh data is shared with the
    other parts, so their slots (and the source's own) are not touched. No material goes back to the data's."""
    if material is None:
        for slot in obj.material_slots:
            slot.link = 'DATA'
        return
    if not obj.material_slots:
        obj.data.materials.append(None)          # the one empty slot every part needs, on the shared mesh
    for slot in obj.material_slots:
        slot.link = 'OBJECT'
        slot.material = material


def apply_settings(src):
    """Push the settings of `src` (width, radius, materials) onto its wire and points copies and onto itself."""
    s = src.render_ready
    wire, points = parts_of(src)
    mod = wire_modifier(wire)
    if mod is not None:
        mod.thickness = s.wire_width
        mod.use_even_offset = s.even_thickness
        set_materials(wire, s.wire_material)
    mod = points_modifier(points)
    if mod is not None:
        set_input(mod, "Radius", s.point_radius)
        set_input(mod, "Material", s.wire_material)
    set_materials(src, s.mesh_material)
