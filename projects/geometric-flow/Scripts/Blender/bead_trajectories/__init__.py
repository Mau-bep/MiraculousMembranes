"""Bead Trajectories: spheres (and force arrows) for the beads of a simulation, keyframed from Bead_<i>_data.txt."""

bl_info = {
    "name": "Bead Trajectories",
    "author": "Mau",
    "version": (0, 1, 2),
    "blender": (4, 2, 0),
    "location": "3D Viewport > Sidebar > Beads",
    "description": "Animated spheres and force arrows for the beads of a results folder, on the frames of its OBJ sequence",
    "category": "Import-Export",
}

import bpy

from . import ops, props, ui

CLASSES = (props.BT_Settings,) + ops.CLASSES + ui.CLASSES


def register():
    for cls in CLASSES:
        bpy.utils.register_class(cls)
    bpy.types.Scene.bead_traj = bpy.props.PointerProperty(type=props.BT_Settings)


def unregister():
    del bpy.types.Scene.bead_traj
    for cls in reversed(CLASSES):
        bpy.utils.unregister_class(cls)
