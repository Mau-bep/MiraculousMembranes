"""Operators of Bead Trajectories: build the animated beads of a results folder, remove them."""

import colorsys
import os

import bmesh
import bpy
import numpy as np
from bpy.types import Operator
from mathutils import Matrix, Quaternion, Vector

from . import data

COLLECTION_PREFIX = "Beads - "
PIVOT_NAME = "Scene Pivot"
ARROW_MATERIAL = "Bead Force Arrows"
SHAFT = 0.72          # the unit arrow: a shaft from 0 to SHAFT, a cone head from SHAFT to 1
HIDE_BELOW = 1e-3     # an arrow shorter than this fraction of the typical force is not drawn


def collection_name(folder):
    folder = os.path.normpath(os.path.abspath(bpy.path.abspath(folder)))
    return COLLECTION_PREFIX + os.path.basename(os.path.dirname(folder)) + "/" + os.path.basename(folder)


def sphere_mesh(subdivisions):
    """One unit icosphere shared by every bead (the radius is the object scale)."""
    name = "Bead Sphere {}".format(subdivisions)
    mesh = bpy.data.meshes.get(name)
    if mesh is None:
        mesh = bpy.data.meshes.new(name)
        bm = bmesh.new()
        bmesh.ops.create_icosphere(bm, subdivisions=subdivisions, radius=1.0)
        bm.to_mesh(mesh)
        bm.free()
        for poly in mesh.polygons:
            poly.use_smooth = True
    return mesh


def arrow_mesh():
    """A unit arrow along +Z, from the origin to z = 1 (the object scale is its length)."""
    mesh = bpy.data.meshes.get("Bead Force Arrow")
    if mesh is None:
        mesh = bpy.data.meshes.new("Bead Force Arrow")
        bm = bmesh.new()
        bmesh.ops.create_cone(bm, cap_ends=True, segments=16, radius1=0.035, radius2=0.035, depth=SHAFT,
                              matrix=Matrix.Translation((0.0, 0.0, SHAFT / 2.0)))
        bmesh.ops.create_cone(bm, cap_ends=True, segments=16, radius1=0.09, radius2=0.0, depth=1.0 - SHAFT,
                              matrix=Matrix.Translation((0.0, 0.0, SHAFT + (1.0 - SHAFT) / 2.0)))
        bm.to_mesh(mesh)
        bm.free()
    return mesh


def material(name, color):
    mat = bpy.data.materials.get(name)
    if mat is None:
        mat = bpy.data.materials.new(name)
        mat.diffuse_color = color
        try:
            mat.use_nodes = True
            bsdf = mat.node_tree.nodes.get("Principled BSDF")
            if bsdf is not None:
                bsdf.inputs["Base Color"].default_value = color
        except (AttributeError, RuntimeError):
            pass
    return mat


def _object_material(obj, mat):
    """A material on the object (object-linked slot), so objects that share a mesh can differ."""
    if not obj.data.materials:
        obj.data.materials.append(None)
    obj.material_slots[0].link = 'OBJECT'
    obj.material_slots[0].material = mat


def _child(parent, obj, location=(0.0, 0.0, 0.0)):
    obj.parent = parent
    obj.matrix_parent_inverse = Matrix.Identity(4)
    obj.location = location


def clear_collection(collection):
    for obj in list(collection.objects):
        bpy.data.objects.remove(obj)


def arrow_rotations(direction):
    """Quaternions taking +Z to each force direction (F, 3); the sign is kept continuous in time so interpolating the components is safe.

    A zero force keeps the last direction (the arrow has zero length there anyway).
    """
    z = Vector((0.0, 0.0, 1.0))
    quats, last = [], Quaternion((1.0, 0.0, 0.0, 0.0))
    for d in direction:
        v = Vector(d)
        if v.length > 1e-12:
            q = z.rotation_difference(v.normalized())
            if q.dot(last) < 0.0:
                q.negate()
            last = q
        quats.append(last.copy())
    return quats


def build(context, settings):
    """Make the collection of beads of settings.folder; returns (n_beads, n_frames, first_step, last_step)."""
    folder = bpy.path.abspath(settings.folder)
    steps, pos, force, file_radii = data.load(folder)
    n_beads, n_frames = pos.shape[0], pos.shape[1]
    radii = [settings.radius if settings.radius > 0.0 else (r if r is not None else 1.0) for r in file_radii]
    frames = settings.start_frame + np.arange(n_frames)

    name = collection_name(folder)
    collection = bpy.data.collections.get(name)
    if collection is None:
        collection = bpy.data.collections.new(name)
    if collection not in context.scene.collection.children_recursive[:]:
        context.scene.collection.children.link(collection)
    clear_collection(collection)

    mesh = sphere_mesh(settings.subdivisions)
    pivot = bpy.data.objects.get(PIVOT_NAME)
    if pivot is not None:
        context.view_layer.update()                            # its world matrix may be stale after a move made by a script
    mag = np.linalg.norm(force, axis=2)
    f_ref = float(np.percentile(mag[mag > 1e-12], 95)) if np.any(mag > 1e-12) else 1.0
    r_typ = float(np.mean(radii))
    length = np.minimum(mag * (settings.arrow_length * r_typ / f_ref), settings.arrow_cap * settings.arrow_length * r_typ)
    length[mag < HIDE_BELOW * f_ref] = 0.0

    edit = context.preferences.edit
    old_interpolation = edit.keyframe_new_interpolation_type
    edit.keyframe_new_interpolation_type = 'LINEAR'          # the beads move in straight lines between the saved steps
    try:
        for i in range(n_beads):
            hue = (0.62 + 0.8 * i / max(n_beads, 1)) % 1.0
            color = colorsys.hsv_to_rgb(hue, 0.55, 0.9) + (1.0,)

            bead = bpy.data.objects.new("Bead {}".format(i), None)
            bead.empty_display_type = 'PLAIN_AXES'
            bead.empty_display_size = radii[i]
            collection.objects.link(bead)
            if pivot is not None:                              # keep the grouping when the beads are rebuilt
                bead.parent = pivot
                bead.matrix_parent_inverse = pivot.matrix_world.inverted()
            sphere = bpy.data.objects.new("Bead {} sphere".format(i), mesh)
            sphere.scale = (radii[i],) * 3
            collection.objects.link(sphere)
            _child(bead, sphere)
            _object_material(sphere, material("Bead {}".format(i), color))

            for f, p in zip(frames, pos[i]):
                bead.location = p
                bead.keyframe_insert("location", frame=int(f))

            if settings.arrows:
                pivot = bpy.data.objects.new("Bead {} force".format(i), None)
                pivot.empty_display_type = 'SINGLE_ARROW'
                pivot.empty_display_size = 0.5 * radii[i]
                pivot.rotation_mode = 'QUATERNION'
                collection.objects.link(pivot)
                _child(bead, pivot)
                arrow = bpy.data.objects.new("Bead {} arrow".format(i), arrow_mesh())
                collection.objects.link(arrow)
                _child(pivot, arrow, (0.0, 0.0, radii[i]))      # starts at the surface of the bead
                _object_material(arrow, material(ARROW_MATERIAL, (0.85, 0.1, 0.08, 1.0)))
                for f, q, L in zip(frames, arrow_rotations(force[i]), length[i]):
                    pivot.rotation_quaternion = q
                    pivot.keyframe_insert("rotation_quaternion", frame=int(f))
                    arrow.scale = (L, L, L)
                    arrow.keyframe_insert("scale", frame=int(f))
    finally:
        edit.keyframe_new_interpolation_type = old_interpolation

    return n_beads, n_frames, int(steps[0]), int(steps[-1])


class BT_OT_load(Operator):
    """Make a sphere for every bead of the results folder, keyframed on the frames of the membrane OBJ sequence"""
    bl_idname = "bead_traj.load"
    bl_label = "Load Beads"
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(cls, context):
        return bool(context.scene.bead_traj.folder)

    def execute(self, context):
        settings = context.scene.bead_traj
        try:
            n_beads, n_frames, first, last = build(context, settings)
        except (ValueError, OSError) as err:
            self.report({'ERROR'}, str(err))
            return {'CANCELLED'}
        self.report({'INFO'}, "{} beads on {} frames (steps {} to {}), frames {} to {}".format(
            n_beads, n_frames, first, last, settings.start_frame, settings.start_frame + n_frames - 1))
        return {'FINISHED'}


class BT_OT_remove(Operator):
    """Delete the beads of the results folder (the collection and its objects)"""
    bl_idname = "bead_traj.remove"
    bl_label = "Remove Beads"
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(cls, context):
        settings = context.scene.bead_traj
        return bool(settings.folder) and bpy.data.collections.get(collection_name(settings.folder)) is not None

    def execute(self, context):
        collection = bpy.data.collections.get(collection_name(context.scene.bead_traj.folder))
        clear_collection(collection)
        bpy.data.collections.remove(collection)
        return {'FINISHED'}


class BT_OT_same_material(Operator):
    """Give every bead sphere the Material of the panel (the arrows keep theirs); Update Beads brings back one colour per bead"""
    bl_idname = "bead_traj.same_material"
    bl_label = "Apply to All Beads"
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(cls, context):
        settings = context.scene.bead_traj
        return (settings.material is not None and bool(settings.folder)
                and bpy.data.collections.get(collection_name(settings.folder)) is not None)

    def execute(self, context):
        settings = context.scene.bead_traj
        collection = bpy.data.collections.get(collection_name(settings.folder))
        spheres = [o for o in collection.objects if o.type == 'MESH' and o.data.name.startswith("Bead Sphere")]
        for sphere in spheres:
            _object_material(sphere, settings.material)
        self.report({'INFO'}, "{} on {} beads".format(settings.material.name, len(spheres)))
        return {'FINISHED'}


def _world_points(obj):
    """The corners of the bounding box of a mesh (or the origin of anything else) in world space."""
    if obj.type == 'MESH' and len(obj.data.vertices):
        return [obj.matrix_world @ Vector(c) for c in obj.bound_box]
    return [obj.matrix_world.translation.copy()]


class BT_OT_pivot(Operator):
    """Parent the beads and the selected objects (membrane, its RenderReady copies, slicer, camera ...) to one pivot empty, so that
moving, rotating or scaling the pivot moves everything together; the bead keyframes stay relative to it"""
    bl_idname = "bead_traj.pivot"
    bl_label = "Group With Pivot"
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(cls, context):
        settings = context.scene.bead_traj
        return bool(settings.folder) and bpy.data.collections.get(collection_name(settings.folder)) is not None

    def execute(self, context):
        settings = context.scene.bead_traj
        collection = bpy.data.collections.get(collection_name(settings.folder))
        pivot = bpy.data.objects.get(PIVOT_NAME)
        beads = [o for o in collection.objects if o.parent is None or o.parent is pivot]
        beads = [o for o in beads if o.type == 'EMPTY' and o.name.startswith("Bead ") and o.name.split()[-1].isdigit()]
        chosen = [o for o in context.selected_objects if o is not pivot and o not in beads]
        candidates = beads + chosen
        top = [o for o in candidates if o.parent is None or o.parent is pivot]      # children move with their parent
        skipped = len(candidates) - len(top)

        if pivot is None:
            points = [p for o in top for p in _world_points(o)]
            lo = Vector((min(p[k] for p in points) for k in range(3)))
            hi = Vector((max(p[k] for p in points) for k in range(3)))
            pivot = bpy.data.objects.new(PIVOT_NAME, None)
            pivot.empty_display_type = 'PLAIN_AXES'
            pivot.empty_display_size = max((hi - lo).length * 0.25, 1.0)
            pivot.location = (lo + hi) / 2.0
            context.scene.collection.objects.link(pivot)
        context.view_layer.update()                             # a new pivot has no world matrix until it is evaluated
        inverse = pivot.matrix_world.inverted()
        moved = 0
        for obj in top:
            if obj.parent is pivot:
                continue
            obj.parent = pivot
            obj.matrix_parent_inverse = inverse                 # stays exactly where it is now
            moved += 1
        for o in context.selected_objects:
            o.select_set(False)
        pivot.select_set(True)
        context.view_layer.objects.active = pivot
        self.report({'INFO'}, "{} objects under {}{}".format(moved, PIVOT_NAME, "; {} with a parent of their own were left".format(skipped) if skipped else ""))
        return {'FINISHED'}


CLASSES = (BT_OT_load, BT_OT_remove, BT_OT_same_material, BT_OT_pivot)
