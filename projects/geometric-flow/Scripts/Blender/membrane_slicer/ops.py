"""Operators for Membrane Slicer."""

import math

import bpy
from bpy.props import IntProperty, StringProperty, EnumProperty
from bpy.types import Operator
from mathutils import Matrix, Vector

from . import nodes


def _meshes(objs):
    return [o for o in objs if o.type == 'MESH' and nodes.cutter_modifier(o) is None]


def _world_bounds(objs):
    """(centre, largest extent) of the objects' world-space bounding boxes."""
    pts = [o.matrix_world @ Vector(c) for o in objs for c in o.bound_box]
    lo = Vector(map(min, zip(*pts)))
    hi = Vector(map(max, zip(*pts)))
    return (lo + hi) / 2, max(hi - lo)


def _view_rotation(context):
    """Slicer orientation that opens toward the viewer, hinge pointing up."""
    rv3d = getattr(context.space_data, "region_3d", None)
    if rv3d is None:
        return Matrix.Identity(3)
    v = rv3d.view_matrix.inverted().to_quaternion().to_matrix()   # pure rotation
    right, up, back = v.col[0], v.col[1], v.col[2]
    return Matrix((back, right, up)).transposed()   # columns: x=back, y=right, z=up


def attach(obj, slicer):
    """Put the MS Cut modifier (first in the stack) on `obj`, pointing at `slicer`."""
    mod = nodes.cut_modifier(obj)
    if mod is None:
        mod = obj.modifiers.new(nodes.CUT_MOD, 'NODES')
        mod.node_group = nodes.ensure_cut_group()
        obj.modifiers.move(len(obj.modifiers) - 1, 0)
    nodes.set_input(mod, "Slicer", slicer)
    return mod


def detach(obj):
    mod = nodes.cut_modifier(obj)
    while mod is not None:
        obj.modifiers.remove(mod)
        mod = nodes.cut_modifier(obj)


def make_slicer(context, targets):
    ng = nodes.ensure_cutter_group()
    obj = bpy.data.objects.new("Slicer", bpy.data.meshes.new("Slicer"))
    context.collection.objects.link(obj)
    mod = obj.modifiers.new(nodes.CUTTER_MOD, 'NODES')
    mod.node_group = ng

    centre, extent = (Vector(context.scene.cursor.location), 4.0)
    if targets:
        centre, extent = _world_bounds(targets)
    reach = max(1.5 * extent, 1e-3)
    nodes.set_input(mod, "Radius", reach)
    nodes.set_input(mod, "Height", reach)
    nodes.set_input(mod, "Frame Start", context.scene.frame_start)
    nodes.set_input(mod, "Frame End", context.scene.frame_end)

    obj.matrix_world = Matrix.Translation(centre) @ _view_rotation(context).to_4x4()
    obj.display_type = 'WIRE'
    obj.hide_render = True
    obj.visible_shadow = False
    obj.visible_camera = False
    return obj


class MS_OT_create(Operator):
    """Create a slicer and cut the selected membrane objects with it"""
    bl_idname = "membrane_slicer.create"
    bl_label = "Create Slicer"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        targets = _meshes(context.selected_objects)
        slicer = make_slicer(context, targets)
        for obj in targets:
            attach(obj, slicer)
        for obj in context.selected_objects:
            obj.select_set(False)
        slicer.select_set(True)
        context.view_layer.objects.active = slicer
        self.report({'INFO'}, f"Slicer cuts {len(targets)} object(s)")
        return {'FINISHED'}


class MS_OT_add_targets(Operator):
    """Cut the other selected objects with the active slicer"""
    bl_idname = "membrane_slicer.add_targets"
    bl_label = "Add Selected"
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(cls, context):
        return nodes.cutter_modifier(context.active_object) is not None

    def execute(self, context):
        slicer = context.active_object
        targets = [o for o in _meshes(context.selected_objects) if o != slicer]
        for obj in targets:
            attach(obj, slicer)
        self.report({'INFO'}, f"Added {len(targets)} object(s)")
        return {'FINISHED'}


class MS_OT_remove_target(Operator):
    """Stop cutting this object (the slicer stays)"""
    bl_idname = "membrane_slicer.remove_target"
    bl_label = "Remove From Slicer"
    bl_options = {'REGISTER', 'UNDO'}

    name: StringProperty()

    def execute(self, context):
        obj = bpy.data.objects.get(self.name)
        if obj is not None:
            detach(obj)
        return {'FINISHED'}


class MS_OT_delete(Operator):
    """Delete the slicer and uncut everything it was cutting"""
    bl_idname = "membrane_slicer.delete"
    bl_label = "Delete Slicer"
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(cls, context):
        return nodes.cutter_modifier(context.active_object) is not None

    def execute(self, context):
        slicer = context.active_object
        for obj in nodes.targets_of(slicer):
            detach(obj)
        mesh = slicer.data
        bpy.data.objects.remove(slicer)
        if mesh.users == 0:
            bpy.data.meshes.remove(mesh)
        return {'FINISHED'}


class _SlicerOp(Operator):
    @classmethod
    def poll(cls, context):
        return nodes.slicer_of(context.active_object) is not None

    def modifier(self, context):
        return nodes.cutter_modifier(nodes.slicer_of(context.active_object))


class MS_OT_plane(_SlicerOp):
    """Make the cut a single plane (180 degree slice)"""
    bl_idname = "membrane_slicer.plane"
    bl_label = "Plane (180°)"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        nodes.set_input(self.modifier(context), "Angle", math.pi)
        return {'FINISHED'}


class MS_OT_keyframe(_SlicerOp):
    """Insert keyframes at the current frame"""
    bl_idname = "membrane_slicer.keyframe"
    bl_label = "Keyframe Slicer"
    bl_options = {'REGISTER', 'UNDO'}

    what: EnumProperty(items=[
        ('ANGLE', "Angle", "Keyframe the wedge angle"),
        ('POSE', "Angle + Position", "Keyframe the wedge angle, location and rotation"),
    ])

    def execute(self, context):
        slicer = nodes.slicer_of(context.active_object)
        paths = [nodes.input_data_path(nodes.cutter_modifier(slicer), "Angle")]
        if self.what == 'POSE':
            rot = {'QUATERNION': "rotation_quaternion", 'AXIS_ANGLE': "rotation_axis_angle"}
            paths += ["location", rot.get(slicer.rotation_mode, "rotation_euler")]
        for path in paths:
            slicer.keyframe_insert(data_path=path)
        return {'FINISHED'}


class MS_OT_fit(_SlicerOp):
    """Resize the cutter so it reaches all of its targets"""
    bl_idname = "membrane_slicer.fit"
    bl_label = "Fit to Targets"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        slicer = nodes.slicer_of(context.active_object)
        targets = nodes.targets_of(slicer)
        if not targets:
            return {'CANCELLED'}
        mod = nodes.cutter_modifier(slicer)
        pts = [slicer.matrix_world.inverted() @ o.matrix_world @ Vector(c)
               for o in targets for c in o.bound_box]
        radius = max(math.hypot(p.x, p.y) for p in pts)
        height = max(abs(p.z) for p in pts)
        nodes.set_input(mod, "Radius", 1.1 * radius + 1e-3)
        nodes.set_input(mod, "Height", 1.1 * height + 1e-3)
        return {'FINISHED'}


class MS_OT_range(_SlicerOp):
    """Set the frame range the cut is active in"""
    bl_idname = "membrane_slicer.range"
    bl_label = "Set Frame Range"
    bl_options = {'REGISTER', 'UNDO'}

    which: EnumProperty(items=[
        ('SCENE', "Scene Range", "Use the scene start and end frame"),
        ('START', "Start = Now", "Start the cut at the current frame"),
        ('END', "End = Now", "End the cut at the current frame"),
    ])

    def execute(self, context):
        mod = self.modifier(context)
        scene = context.scene
        if self.which == 'SCENE':
            nodes.set_input(mod, "Frame Start", scene.frame_start)
            nodes.set_input(mod, "Frame End", scene.frame_end)
        elif self.which == 'START':
            nodes.set_input(mod, "Frame Start", scene.frame_current)
        else:
            nodes.set_input(mod, "Frame End", scene.frame_current)
        nodes.set_input(mod, "Limit to Frames", True)
        return {'FINISHED'}


CLASSES = (
    MS_OT_create,
    MS_OT_add_targets,
    MS_OT_remove_target,
    MS_OT_delete,
    MS_OT_plane,
    MS_OT_keyframe,
    MS_OT_fit,
    MS_OT_range,
)
