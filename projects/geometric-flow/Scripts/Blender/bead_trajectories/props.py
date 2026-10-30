"""The settings of the add-on, kept on the scene."""

import bpy
from bpy.props import BoolProperty, FloatProperty, IntProperty, PointerProperty, StringProperty
from bpy.types import PropertyGroup


class BT_Settings(PropertyGroup):
    folder: StringProperty(
        name="Results Folder", subtype='DIR_PATH', default="",
        description="Folder of one run: the Bead_<i>_data.txt files, Input_file.json (radii) and the membrane_<step>.obj meshes")
    start_frame: IntProperty(
        name="Start Frame", default=1,
        description="Frame of the first saved mesh: set it to the first frame of the OBJ sequence of the membrane")
    subdivisions: IntProperty(
        name="Sphere Detail", default=3, min=1, max=6,
        description="Subdivisions of the icosphere of every bead (3 = 642 vertices, 4 = 2562)")
    radius: FloatProperty(
        name="Radius", default=0.0, min=0.0, soft_max=10.0, precision=3, subtype='DISTANCE',
        description="Bead radius; 0 reads the radius of every bead from Input_file.json (1 when there is none)")
    material: PointerProperty(
        name="Material", type=bpy.types.Material,
        description="Material for Apply to All Beads (the spheres start with a colour of their own each)")
    arrows: BoolProperty(
        name="Force Arrows", default=False,
        description="An arrow on every bead along the force on it (the Fx Fy Fz columns), starting at its surface; zero forces show nothing")
    arrow_length: FloatProperty(
        name="Arrow Length", default=2.0, min=0.0, soft_max=10.0, precision=2,
        description="Length, in bead radii, of the arrow of a typical force (the 95th percentile of all the forces of the run)")
    arrow_cap: FloatProperty(
        name="Arrow Cap", default=3.0, min=1.0, soft_max=10.0, precision=2,
        description="The longest arrow is this many times the Arrow Length, so a few huge forces do not fill the picture")
