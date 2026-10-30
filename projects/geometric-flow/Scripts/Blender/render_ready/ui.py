"""Sidebar panel for RenderReady."""

from bpy.types import Panel

from . import nodes


class RR_PT_main(Panel):
    bl_label = "RenderReady"
    bl_idname = "RR_PT_main"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = "RenderReady"

    def draw(self, context):
        layout = self.layout
        layout.operator("render_ready.create", icon='MOD_WIREFRAME')
        src = nodes.source_of(context.active_object)
        if src is None:
            layout.label(text="Select a mesh first", icon='INFO')
            return

        settings = src.render_ready
        layout.label(text=f"{src.name}: wireframe + points", icon='MESH_DATA')

        col = layout.column(align=True)
        col.label(text="Mesh Material")
        col.template_ID(settings, "mesh_material", new="material.new")
        col = layout.column(align=True)
        col.label(text="Wire and Points Material")
        col.template_ID(settings, "wire_material", new="material.new")

        col = layout.column(align=True)
        col.prop(settings, "wire_width")
        col.prop(settings, "point_radius")
        col.prop(settings, "even_thickness")

        layout.operator("render_ready.remove", icon='TRASH')


CLASSES = (RR_PT_main,)
