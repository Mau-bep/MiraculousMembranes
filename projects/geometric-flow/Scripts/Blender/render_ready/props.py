"""The settings RenderReady keeps on the source object."""

import bpy
from bpy.props import BoolProperty, FloatProperty, PointerProperty
from bpy.types import PropertyGroup

from . import nodes


def _apply(self, context):
    nodes.apply_settings(self.id_data)


class RR_Settings(PropertyGroup):
    wire_material: PointerProperty(
        name="Wire Material", type=bpy.types.Material, update=_apply,
        description="Material of the wireframe and of the points")
    mesh_material: PointerProperty(
        name="Mesh Material", type=bpy.types.Material, update=_apply,
        description="Material of the original mesh")
    wire_width: FloatProperty(
        name="Wire Width", default=0.01, min=0.0, soft_max=1.0, precision=4, subtype='DISTANCE', update=_apply,
        description="Thickness of the wireframe")
    even_thickness: BoolProperty(
        name="Even Thickness", default=False, update=_apply,
        description="Wireframe Even Thickness: keeps the wire width constant at sharp corners, but on the thin sliver triangles a cut leaves "
                    "(Membrane Slicer) it throws wire vertices far off the surface: long stray lines. Leave it off for cut membranes")
    point_radius: FloatProperty(
        name="Point Radius", default=0.02, min=0.0, soft_max=1.0, precision=4, subtype='DISTANCE', update=_apply,
        description="Radius of the points")
