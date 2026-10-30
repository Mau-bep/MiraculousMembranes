"""RenderReady: a wireframe and a set of vertex points on top of a mesh, with materials you choose."""

bl_info = {
    "name": "RenderReady",
    "author": "Mau",
    "version": (0, 1, 1),
    "blender": (4, 2, 0),
    "location": "3D Viewport > Sidebar > RenderReady",
    "description": "Two live copies of a mesh, one as a wireframe and one as points, plus a material for the original",
    "category": "Object",
}

import bpy

from . import ops, props, ui

CLASSES = (props.RR_Settings,) + ops.CLASSES + ui.CLASSES


def register():
    for cls in CLASSES:
        bpy.utils.register_class(cls)
    bpy.types.Object.render_ready = bpy.props.PointerProperty(type=props.RR_Settings)


def unregister():
    del bpy.types.Object.render_ready
    for cls in reversed(CLASSES):
        bpy.utils.unregister_class(cls)
