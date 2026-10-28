"""Membrane Normals: render-ready normal, mean curvature and force arrows for a membrane mesh."""

bl_info = {
    "name": "Membrane Normals",
    "author": "Mau",
    "version": (0, 2, 0),
    "blender": (4, 2, 0),
    "location": "3D Viewport > Sidebar > Normals",
    "description": "Arrows for the normals, the mean curvature normals or the dumped simulation forces of a mesh, colour coded by magnitude",
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
