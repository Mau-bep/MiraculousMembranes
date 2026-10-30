"""Sidebar panel of Bead Trajectories."""

import bpy
from bpy.types import Panel

from . import ops


class BT_PT_main(Panel):
    bl_label = "Bead Trajectories"
    bl_idname = "BT_PT_main"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = "Beads"

    def draw(self, context):
        layout = self.layout
        s = context.scene.bead_traj
        layout.prop(s, "folder")
        col = layout.column(align=True)
        col.prop(s, "start_frame")
        col.prop(s, "subdivisions")
        col.prop(s, "radius")
        layout.prop(s, "arrows")
        if s.arrows:
            col = layout.column(align=True)
            col.prop(s, "arrow_length")
            col.prop(s, "arrow_cap")
        exists = bool(s.folder) and bpy.data.collections.get(ops.collection_name(s.folder)) is not None
        layout.operator("bead_traj.load", text="Update Beads" if exists else "Load Beads", icon='MESH_UVSPHERE')
        if exists:
            box = layout.box()
            box.label(text="Same material for all beads")
            box.template_ID(s, "material", new="material.new")
            box.operator("bead_traj.same_material", icon='MATERIAL')
            box = layout.box()
            box.label(text="Move everything together")
            box.label(text="Select the membrane (and slicer), then:")
            box.operator("bead_traj.pivot", icon='EMPTY_AXIS')
            layout.operator("bead_traj.remove", icon='TRASH')


CLASSES = (BT_PT_main,)
