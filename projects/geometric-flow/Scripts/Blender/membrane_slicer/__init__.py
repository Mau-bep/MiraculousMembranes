"""Membrane Slicer: cut away a plane or a cake slice of OBJ-sequence membranes."""

bl_info = {
    "name": "Membrane Slicer",
    "author": "Mau",
    "version": (0, 1, 1),
    "blender": (4, 2, 0),
    "location": "3D Viewport > Sidebar > Slicer",
    "description": "Non-destructive plane / wedge cut for membranes, works on every frame of an OBJ sequence",
    "category": "Object",
}

import bpy

from . import nodes, ops, ui


def register():
    for cls in ops.CLASSES + ui.CLASSES:
        bpy.utils.register_class(cls)


def unregister():
    for cls in reversed(ops.CLASSES + ui.CLASSES):
        bpy.utils.unregister_class(cls)
