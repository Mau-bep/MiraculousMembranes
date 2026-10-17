"""Sidebar panel for Membrane Slicer."""

import bpy
from bpy.types import Panel

from . import nodes


class MS_PT_main(Panel):
    bl_label = "Membrane Slicer"
    bl_idname = "MS_PT_main"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = "Slicer"

    def draw(self, context):
        layout = self.layout
        slicer = nodes.slicer_of(context.active_object)

        layout.operator("membrane_slicer.create", icon='MOD_BOOLEAN')
        if slicer is None:
            layout.label(text="Select membrane objects, then Create", icon='INFO')
            return

        mod = nodes.cutter_modifier(slicer)
        layout.label(text=slicer.name, icon='MESH_CONE')
        if context.active_object != slicer:
            layout.label(text="Select the slicer object to edit it", icon='RESTRICT_SELECT_OFF')
            return

        col = layout.column(align=True)
        nodes.draw_input(col, mod, "Angle")
        col.operator("membrane_slicer.plane", icon='AXIS_TOP')
        row = col.row(align=True)
        row.operator("membrane_slicer.keyframe", text="Key Angle", icon='KEYFRAME_HLT').what = 'ANGLE'
        row.operator("membrane_slicer.keyframe", text="Key Pose", icon='KEYFRAME').what = 'POSE'
        col.separator()
        nodes.draw_input(col, mod, "Radius")
        nodes.draw_input(col, mod, "Height")
        col.operator("membrane_slicer.fit", icon='FULLSCREEN_ENTER')

        box = layout.box()
        nodes.draw_input(box, mod, "Enabled")
        nodes.draw_input(box, mod, "Limit to Frames")
        if nodes.get_input(mod, "Limit to Frames"):
            row = box.row(align=True)
            nodes.draw_input(row, mod, "Frame Start", text="Start")
            nodes.draw_input(row, mod, "Frame End", text="End")
            row = box.row(align=True)
            op = row.operator("membrane_slicer.range", text="Scene")
            op.which = 'SCENE'
            op = row.operator("membrane_slicer.range", text="Start = Now")
            op.which = 'START'
            op = row.operator("membrane_slicer.range", text="End = Now")
            op.which = 'END'
        else:
            box.operator("membrane_slicer.range", text="Limit to Scene Range").which = 'SCENE'

        layout.label(text="Cut Membranes")
        box = layout.box()
        for obj in nodes.targets_of(slicer):
            row = box.row()
            row.label(text=obj.name, icon='MESH_DATA')
            row.operator("membrane_slicer.remove_target", text="", icon='X').name = obj.name
        if not nodes.targets_of(slicer):
            box.label(text="None yet")
        layout.operator("membrane_slicer.add_targets", icon='ADD')

        col = layout.column(align=True)
        col.label(text="Options")
        nodes.draw_input(col, mod, "Keep Cut Faces")
        nodes.draw_input(col, mod, "Hole Tolerant")
        nodes.draw_input(col, mod, "Self Intersections")
        layout.operator("membrane_slicer.delete", icon='TRASH')


CLASSES = (MS_PT_main,)
