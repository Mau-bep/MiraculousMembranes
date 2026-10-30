"""Operators for RenderReady."""

import bpy
import numpy as np
from bpy.types import Operator
from mathutils import Matrix

from . import nodes

COLLECTION_PREFIX = "RenderReady - "
# first width and radius, as a fraction of the mean edge length of the mesh
WIRE_FRACTION = 0.06
POINT_FRACTION = 0.10


def mean_edge_length(mesh):
    if len(mesh.edges) == 0:
        return 0.1
    co = np.empty(len(mesh.vertices) * 3)
    mesh.vertices.foreach_get("co", co)
    co = co.reshape(-1, 3)
    ends = np.empty(len(mesh.edges) * 2, dtype=np.int64)
    mesh.edges.foreach_get("vertices", ends)
    ends = ends.reshape(-1, 2)
    length = float(np.linalg.norm(co[ends[:, 0]] - co[ends[:, 1]], axis=1).mean())
    return length if length > 0.0 else 0.1


def _make_part(context, src, collection, suffix):
    """A copy of `src` that shares its mesh data and modifiers, sits where `src` sits and follows it."""
    part = src.copy()                            # Object.copy() keeps the mesh data: linked, not duplicated
    part.name = f"{src.name} {suffix}"
    part.animation_data_clear()
    for constraint in list(part.constraints):
        part.constraints.remove(constraint)
    part.parent = src
    part.matrix_parent_inverse = Matrix.Identity(4)
    part.location = (0.0, 0.0, 0.0)
    part.rotation_euler = (0.0, 0.0, 0.0)
    part.rotation_quaternion = (1.0, 0.0, 0.0, 0.0)
    part.rotation_axis_angle = (0.0, 0.0, 1.0, 0.0)
    part.scale = (1.0, 1.0, 1.0)
    part.delta_location = (0.0, 0.0, 0.0)
    part.delta_rotation_euler = (0.0, 0.0, 0.0)
    part.delta_rotation_quaternion = (1.0, 0.0, 0.0, 0.0)
    part.delta_scale = (1.0, 1.0, 1.0)
    # An OBJ sequence updates the mesh of every object flagged as a sequence: the copies share the mesh the
    # source's own sequence updates, so they must not update it again
    sequence = getattr(part, "smo_sequence", None)
    if sequence is not None and sequence.is_sequence:
        sequence.is_sequence = False
    part["rr_source"] = src
    collection.objects.link(part)
    return part


def make_render_ready(context, src):
    name = COLLECTION_PREFIX + src.name
    collection = bpy.data.collections.get(name)
    if collection is None:
        collection = bpy.data.collections.new(name)
    if collection not in context.scene.collection.children_recursive[:]:
        context.scene.collection.children.link(collection)

    wire = _make_part(context, src, collection, "wire")
    points = _make_part(context, src, collection, "points")

    mod = wire.modifiers.new(nodes.WIRE_MOD, 'WIREFRAME')      # after the modifiers copied from the source
    mod.use_replace = True
    mod.use_boundary = True
    mod.use_even_offset = src.render_ready.even_thickness      # off: it makes spikes on the slivers of a Membrane Slicer cut
    mod.use_relative_offset = False
    mod = points.modifiers.new(nodes.POINTS_MOD, 'NODES')
    mod.node_group = nodes.ensure_points_group()

    src["rr_wire"] = wire
    src["rr_points"] = points
    edge = mean_edge_length(src.data)
    settings = src.render_ready
    settings.wire_width = WIRE_FRACTION * edge                 # each assignment pushes onto the parts
    settings.point_radius = POINT_FRACTION * edge
    return wire, points


def _select_only(context, obj):
    for o in context.selected_objects:
        o.select_set(False)
    obj.select_set(True)
    context.view_layer.objects.active = obj


class RR_OT_create(Operator):
    """Make a wireframe copy and a points copy of the selected meshes, sharing their mesh (and so following an OBJ sequence)"""
    bl_idname = "render_ready.create"
    bl_label = "Make Render Ready"
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(cls, context):
        return any(o.type == 'MESH' for o in context.selected_objects)

    def execute(self, context):
        made = None
        for obj in list(context.selected_objects):
            if obj.type != 'MESH' or nodes.source_of(obj) not in (None, obj):
                continue                                       # not a mesh, or already a wire/points copy
            if nodes.source_of(obj) is obj:
                self.report({'INFO'}, f"{obj.name} is already render ready")
                made = obj
                continue
            make_render_ready(context, obj)
            made = obj
        if made is None:
            self.report({'ERROR'}, "Select a mesh that is not a RenderReady copy")
            return {'CANCELLED'}
        _select_only(context, made)
        return {'FINISHED'}


class RR_OT_remove(Operator):
    """Delete the wireframe and points copies (the mesh keeps the Mesh Material)"""
    bl_idname = "render_ready.remove"
    bl_label = "Remove Copies"
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(cls, context):
        return nodes.source_of(context.active_object) is not None

    def execute(self, context):
        src = nodes.source_of(context.active_object)
        collections = set()
        for part in nodes.parts_of(src):
            if part is not None:
                collections.update(c for c in part.users_collection if c.name.startswith(COLLECTION_PREFIX))
                bpy.data.objects.remove(part)                  # the mesh data stays: it is the source's
        for key in ("rr_wire", "rr_points"):
            if key in src:
                del src[key]
        for collection in collections:
            if not collection.objects and not collection.children:
                bpy.data.collections.remove(collection)
        _select_only(context, src)
        return {'FINISHED'}


CLASSES = (RR_OT_create, RR_OT_remove)
