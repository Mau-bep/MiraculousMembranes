"""Sidebar panel for Membrane Normals."""

import os

import bpy
import numpy as np
from bpy.types import Panel

from . import nodes


class NA_PT_main(Panel):
    bl_label = "Membrane Normals"
    bl_idname = "NA_PT_main"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = "Normals"

    def draw(self, context):
        layout = self.layout
        obj = context.active_object
        mod = nodes.arrow_modifier(obj)

        col = layout.column(align=True)
        row = col.row(align=True)
        row.operator("normal_arrows.add", text="Normals", icon='ORIENTATION_NORMAL').mean_curvature = False
        row.operator("normal_arrows.add", text="Mean Curvature", icon='FORCE_FORCE').mean_curvature = True
        col.operator("normal_arrows.load_forces", text="Load Forces", icon='IMPORT')
        if mod is None:
            layout.label(text="Select membrane meshes, then add arrows", icon='INFO')
            return

        kind = obj.get("ne_kind", 'NORMALS')
        source = nodes.source_of(obj)
        layout.label(text=f"On {source.name if source else 'a deleted mesh'}", icon='MESH_DATA')
        if kind == 'FORCE':
            layout.label(text=f"{obj.get('ne_term')}, step {obj.get('ne_step')}: "
                              f"{os.path.basename(str(obj.get('ne_file')))}", icon='FILE_TEXT')

        col = layout.column(align=True)
        nodes.draw_input(col, mod, "Offset")
        nodes.draw_input(col, mod, "Length")
        nodes.draw_input(col, mod, "Thickness")
        nodes.draw_input(col, mod, "Head Length")
        nodes.draw_input(col, mod, "Flip")
        colored = bool(nodes.get_input(mod, "Color by Magnitude"))
        if colored:
            nodes.draw_input(col, mod, "Size by Magnitude")

        col = layout.column(align=True)
        nodes.draw_input(col, mod, "Density", slider=True)
        if nodes.get_input(mod, "Density") < 1.0:
            nodes.draw_input(col, mod, "Seed")

        if colored:
            symbol = "|F|" if kind == 'FORCE' else "|H|"
            box = layout.box()
            box.label(text=f"Color by {symbol}")
            row = box.row(align=True)
            nodes.draw_input(row, mod, "Range Min", text="Min")
            nodes.draw_input(row, mod, "Range Max", text="Max")
            nodes.draw_input(box, mod, "Hide Below")
            row = box.row(align=True)
            row.operator("normal_arrows.fit_range", icon='FULLSCREEN_ENTER')
            mat = bpy.data.materials.get(nodes.COLORED_MAT)
            key = mat.get("na_colormap") if mat is not None else None
            row.operator_menu_enum("normal_arrows.colormap", "colormap",
                                   text=nodes.COLORMAPS[key][0] if key in nodes.COLORMAPS else "Colormap")
            stats = obj.get("ne_stats")
            if stats is not None:
                box.label(text=f"Mesh {symbol}: {stats[0]:.3g} to {stats[2]:.3g}")
                cutoff = nodes.get_input(mod, "Hide Below")
                if cutoff:
                    values = np.empty(len(obj.data.vertices))
                    obj.data.attributes[nodes.MAGNITUDE_ATTR].data.foreach_get("value", values)
                    hidden = int(np.count_nonzero(values < cutoff))
                    box.label(text=f"Hidden below the cutoff: {hidden} of {len(values)}",
                              icon='INFO' if hidden == len(values) else 'NONE')
            view = context.scene.view_settings
            if view.view_transform != 'Standard':
                box.label(text="Not Standard: colors shift", icon='ERROR')
                box.prop(view, "view_transform", text="")

        row = layout.row(align=True)
        row.operator("normal_arrows.update", text="Reload File" if kind == 'FORCE' else "Update from Mesh",
                     icon='FILE_REFRESH')
        row.operator("normal_arrows.delete", text="", icon='TRASH')


CLASSES = (NA_PT_main,)
