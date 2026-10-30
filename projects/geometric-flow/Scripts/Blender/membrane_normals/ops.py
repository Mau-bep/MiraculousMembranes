"""Operators for Membrane Normals."""

import os
import re

import bpy
import numpy as np
from bpy.props import BoolProperty, CollectionProperty, EnumProperty, StringProperty
from bpy.types import Operator, OperatorFileListElement
from bpy_extras.io_utils import ImportHelper

from . import compute, nodes

NORMALS_PREFIX = "Normals - "
FORCES_PREFIX = "Forces - "
COLLECTION_PREFIXES = (NORMALS_PREFIX, FORCES_PREFIX)

# Forces/<term>_<step>.txt, as written by the simulation's "dump_forces"
FORCE_FILE = re.compile(r"^(?P<term>.+)_(?P<step>\d+)\.txt$")
# Colored arrows below this fraction of the top of the color range are hidden at first ("Hide Below"): a force of
# zero, like that of a bead that does not touch the membrane, is not worth an arrow
CUTOFF_FRACTION = 1e-3
# The .obj files have 6 significant digits: a mesh further off than this (fraction of its size) is another state
MESH_TOLERANCE = 2e-5


def _sources(objs):
    return [o for o in objs if o.type == 'MESH' and nodes.arrow_modifier(o) is None]


def _arrow_sets(objs):
    return [o for o in objs if nodes.arrow_modifier(o) is not None]


def _target_mesh(context):
    """The membrane the active object stands for: itself, or the source of an arrow set."""
    obj = context.active_object
    if obj is None:
        return None
    if nodes.arrow_modifier(obj) is not None:
        return nodes.source_of(obj)
    return obj if obj.type == 'MESH' else None


# ---------------------------------------------------------------------------
# Measuring
# ---------------------------------------------------------------------------

def read_mesh(context, source, use_modifiers):
    """Positions, vertex normals and triangles of `source` (object space)."""
    if use_modifiers:
        evaluated = source.evaluated_get(context.evaluated_depsgraph_get())
        mesh = evaluated.to_mesh()
    else:
        if source.mode == 'EDIT':
            source.update_from_editmode()
        mesh = source.data
    try:
        n = len(mesh.vertices)
        V = np.empty(n * 3)
        mesh.vertices.foreach_get("co", V)
        N = np.empty(n * 3)
        mesh.vertex_normals.foreach_get("vector", N)
        mesh.calc_loop_triangles()
        T = np.empty(len(mesh.loop_triangles) * 3, dtype=np.int64)
        mesh.loop_triangles.foreach_get("vertices", T)
    finally:
        if use_modifiers:
            evaluated.to_mesh_clear()
    return V.reshape(-1, 3), N.reshape(-1, 3), T.reshape(-1, 3)


def measure(context, source, kind, use_modifiers):
    """Arrow data for plain normals ('NORMALS') or mean curvature normals ('CURVATURE')."""
    V, N, T = read_mesh(context, source, use_modifiers)
    if kind == 'CURVATURE':
        h = compute.signed_mean_curvature(V, T, N)
        # the arrow follows the area gradient: along the normal where H > 0, against it where H < 0
        direction = N * np.where(h >= 0.0, 1.0, -1.0)[:, None]
        magnitude = np.abs(h)
    else:
        direction, magnitude = N, np.zeros(len(V))
    return {"co": V, "normal": N, "direction": direction, "magnitude": magnitude,
            "edge": compute.mean_edge_length(V, T)}


def read_force_file(path):
    """The (n, 3) forces of a Forces/<term>_<step>.txt file."""
    try:
        F = np.loadtxt(path, comments="#", ndmin=2)
    except (OSError, ValueError) as e:
        raise ValueError(f"Cannot read {os.path.basename(path)}: {e}")
    if F.ndim != 2 or F.shape[1] != 3:
        raise ValueError(f"{os.path.basename(path)} does not have 3 numbers per line")
    return F


def measure_force(context, source, path):
    """Arrow data for a force file on the vertices of `source`'s own mesh (the file's vertex order)."""
    F = read_force_file(path)
    V, N, T = read_mesh(context, source, False)
    if len(F) != len(V):
        raise ValueError(f"{os.path.basename(path)} has {len(F)} forces but {source.name} has {len(V)} "
                         f"vertices: show the mesh saved at that step (membrane_<step>.obj)")
    magnitude = np.linalg.norm(F, axis=1)
    direction = N.copy()                         # a vertex without force keeps a direction
    moving = magnitude > 0.0
    direction[moving] = F[moving] / magnitude[moving, None]
    return {"co": V, "normal": N, "direction": direction, "magnitude": magnitude,
            "edge": compute.mean_edge_length(V, T)}


def mesh_difference(source, force_path, step):
    """How far the mesh of `source` is from the membrane_<step>.obj saved at the same step, as a fraction of
    its size; None when that file is not next to the Forces folder or has another vertex count."""
    obj_path = os.path.join(os.path.dirname(os.path.dirname(force_path)), f"membrane_{step}.obj")
    if not os.path.isfile(obj_path):
        return None
    with open(obj_path) as fh:
        saved = np.array([line.split()[1:4] for line in fh if line.startswith("v ")], dtype=float)
    mesh = source.data
    if len(saved) != len(mesh.vertices) or len(saved) == 0:
        return None
    shown = np.empty(len(saved) * 3)
    mesh.vertices.foreach_get("co", shown)
    extent = float((saved.max(axis=0) - saved.min(axis=0)).max())
    return float(np.abs(shown.reshape(-1, 3) - saved).max()) / max(extent, 1e-30)


def fill_mesh(mesh, data):
    """(Re)write the point cloud and its attributes into `mesh`."""
    mesh.clear_geometry()
    mesh.vertices.add(len(data["co"]))
    mesh.vertices.foreach_set("co", data["co"].ravel())
    for name, key in ((nodes.NORMAL_ATTR, "normal"), (nodes.DIRECTION_ATTR, "direction")):
        attr = mesh.attributes.new(name, 'FLOAT_VECTOR', 'POINT')
        attr.data.foreach_set("vector", data[key].ravel())
    attr = mesh.attributes.new(nodes.MAGNITUDE_ATTR, 'FLOAT', 'POINT')
    attr.data.foreach_set("value", data["magnitude"])
    mesh.update()


def store_stats(obj, data):
    """min, mean, max of the magnitude, and the top of the default colour range."""
    m = data["magnitude"]
    _, fit = compute.default_range(m)
    lo, mean, hi = (float(m.min()), float(m.mean()), float(m.max())) if m.size else (0.0, 0.0, 0.0)
    obj["ne_stats"] = [lo, mean, hi, fit]


# ---------------------------------------------------------------------------
# Arrow sets
# ---------------------------------------------------------------------------

def _collection_for(context, source, prefix):
    name = prefix + source.name
    coll = bpy.data.collections.get(name)
    if coll is None:
        coll = bpy.data.collections.new(name)
        context.scene.collection.children.link(coll)
    elif coll not in context.scene.collection.children_recursive[:]:
        context.scene.collection.children.link(coll)
    return coll


def make_arrows(context, source, kind, name, data, use_modifiers=False, **props):
    """A new arrow set for `source`; `props` become custom properties of the object."""
    mesh = bpy.data.meshes.new(name)
    fill_mesh(mesh, data)
    obj = bpy.data.objects.new(name, mesh)
    prefix = FORCES_PREFIX if kind == 'FORCE' else NORMALS_PREFIX
    _collection_for(context, source, prefix).objects.link(obj)
    obj.parent = source                      # identity local matrix: same place as the membrane
    obj["ne_source"] = source
    obj["ne_kind"] = kind
    obj["ne_use_modifiers"] = use_modifiers
    for key, value in props.items():
        obj[key] = value
    store_stats(obj, data)

    mod = obj.modifiers.new(nodes.ARROWS_MOD, 'NODES')
    mod.node_group = nodes.ensure_arrows_group()
    nodes.set_input(mod, "Color by Magnitude", kind != 'NORMALS')
    nodes.set_input(mod, "Length", 1.5 * data["edge"])
    lo, hi = compute.default_range(data["magnitude"])
    nodes.set_input(mod, "Range Min", lo)
    nodes.set_input(mod, "Range Max", hi)
    nodes.set_input(mod, "Hide Below", CUTOFF_FRACTION * hi)
    nodes.set_input(mod, "Normal Material", nodes.ensure_normal_material())
    nodes.set_input(mod, "Colored Material", nodes.ensure_colored_material())
    return obj


def refresh(context, arrows):
    """Measure the membrane (or read the force file) again into an existing arrow set."""
    source = nodes.source_of(arrows)
    kind = arrows.get("ne_kind", 'NORMALS')
    if kind == 'FORCE':
        data = measure_force(context, source, arrows["ne_file"])
    else:
        data = measure(context, source, kind, bool(arrows.get("ne_use_modifiers", True)))
    fill_mesh(arrows.data, data)
    store_stats(arrows, data)


def _select_only(context, obj):
    for o in context.selected_objects:
        o.select_set(False)
    obj.select_set(True)
    context.view_layer.objects.active = obj


class NA_OT_add(Operator):
    """Add normal arrows to the selected membrane meshes, in their own collection"""
    bl_idname = "normal_arrows.add"
    bl_label = "Add Arrows"
    bl_options = {'REGISTER', 'UNDO'}

    mean_curvature: BoolProperty(
        name="Mean Curvature Normals", default=False,
        description="Arrows follow the area gradient and are coloured by the mean curvature magnitude")
    use_modifiers: BoolProperty(
        name="Use Modifiers", default=True,
        description="Take the mesh as evaluated (a slicer cut, subdivision, ...) instead of the base mesh")

    def execute(self, context):
        sources = _sources(context.selected_objects)
        if not sources:
            # an arrow set is active: add the other kind for its membrane
            source = _target_mesh(context)
            sources = [source] if source is not None else []
        if not sources:
            self.report({'ERROR'}, "Select a membrane mesh first")
            return {'CANCELLED'}

        kind = 'CURVATURE' if self.mean_curvature else 'NORMALS'
        made = None
        for source in sources:
            existing = [o for o in nodes.arrow_sets_of(source) if o.get("ne_kind") == kind]
            if existing:
                made = existing[0]
                self.report({'INFO'}, f"{made.name} already exists: updated it")
                refresh(context, made)
            else:
                label = "Mean Curvature" if self.mean_curvature else "Normals"
                data = measure(context, source, kind, self.use_modifiers)
                made = make_arrows(context, source, kind, f"{source.name} {label}", data,
                                   use_modifiers=self.use_modifiers)
        _select_only(context, made)
        return {'FINISHED'}


class NA_OT_load_forces(Operator, ImportHelper):
    """Load force files written by the simulation ("dump_forces") as arrows on the active membrane"""
    bl_idname = "normal_arrows.load_forces"
    bl_label = "Load Forces"
    bl_options = {'REGISTER', 'UNDO'}

    filename_ext = ".txt"
    filter_glob: StringProperty(default="*.txt", options={'HIDDEN'})
    files: CollectionProperty(type=OperatorFileListElement, options={'HIDDEN', 'SKIP_SAVE'})
    directory: StringProperty(subtype='DIR_PATH', options={'HIDDEN'})
    all_terms: BoolProperty(
        name="All Terms of This Step", default=False,
        description="Also load every other <term>_<step>.txt of the same step found in that folder")

    @classmethod
    def poll(cls, context):
        return _target_mesh(context) is not None

    def paths(self):
        folder = bpy.path.abspath(self.directory) if self.directory else os.path.dirname(self.filepath)
        chosen = [os.path.join(folder, f.name) for f in self.files if f.name] or [self.filepath]
        if self.all_terms:
            steps = {m.group("step") for m in map(FORCE_FILE.match, map(os.path.basename, chosen)) if m}
            for name in sorted(os.listdir(folder)):
                m = FORCE_FILE.match(name)
                if m and m.group("step") in steps and os.path.join(folder, name) not in chosen:
                    chosen.append(os.path.join(folder, name))
        return chosen

    def execute(self, context):
        source = _target_mesh(context)
        if source is None:
            self.report({'ERROR'}, "Select the membrane mesh first")
            return {'CANCELLED'}
        made = []
        for path in self.paths():
            match = FORCE_FILE.match(os.path.basename(path))
            if not match:
                self.report({'WARNING'}, f"{os.path.basename(path)}: not named <term>_<step>.txt, skipped")
                continue
            term, step = match.group("term"), int(match.group("step"))
            try:
                data = measure_force(context, source, path)
            except ValueError as e:
                self.report({'ERROR'}, str(e))
                continue
            off = mesh_difference(source, path, step)
            if off is not None and off > MESH_TOLERANCE:
                self.report({'WARNING'}, f"{source.name} differs from membrane_{step}.obj by {off:.1e} of its size: "
                                         f"is the mesh on the frame of step {step}?")
            existing = [o for o in nodes.arrow_sets_of(source)
                        if o.get("ne_kind") == 'FORCE' and o.get("ne_term") == term and o.get("ne_step") == step]
            if existing:
                obj = existing[0]
                obj["ne_file"] = path
                fill_mesh(obj.data, data)
                store_stats(obj, data)
            else:
                obj = make_arrows(context, source, 'FORCE', f"{source.name} {term} {step}", data,
                                  ne_file=path, ne_term=term, ne_step=step)
            made.append(obj)
        if not made:
            return {'CANCELLED'}
        _select_only(context, made[0])
        self.report({'INFO'}, f"Loaded {len(made)} force set(s) on {source.name}")
        return {'FINISHED'}


class _ArrowOp(Operator):
    @classmethod
    def poll(cls, context):
        return nodes.arrow_modifier(context.active_object) is not None

    def modifier(self, context):
        return nodes.arrow_modifier(context.active_object)


class NA_OT_update(_ArrowOp):
    """Measure the membrane again (after a frame change or an edit), or read the force file again, and rebuild the arrows"""
    bl_idname = "normal_arrows.update"
    bl_label = "Update from Mesh"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        done = 0
        for arrows in _arrow_sets(set(context.selected_objects) | {context.active_object}):
            if nodes.source_of(arrows) is None:
                self.report({'WARNING'}, f"{arrows.name}: its membrane was deleted")
                continue
            try:
                refresh(context, arrows)
            except ValueError as e:
                self.report({'ERROR'}, str(e))
                continue
            done += 1
        self.report({'INFO'}, f"Updated {done} arrow set(s)")
        return {'FINISHED'} if done else {'CANCELLED'}


class NA_OT_fit_range(_ArrowOp):
    """Set the colour range from the magnitudes of this arrow set (0 to their 98th percentile)"""
    bl_idname = "normal_arrows.fit_range"
    bl_label = "Fit Range"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        stats = context.active_object.get("ne_stats")
        if stats is None:
            return {'CANCELLED'}
        nodes.set_input(self.modifier(context), "Range Min", 0.0)
        nodes.set_input(self.modifier(context), "Range Max", stats[3])
        return {'FINISHED'}


class NA_OT_colormap(Operator):
    """Choose the colormap of the colored arrows (shared by all arrow sets)"""
    bl_idname = "normal_arrows.colormap"
    bl_label = "Colormap"
    bl_options = {'REGISTER', 'UNDO'}

    colormap: EnumProperty(name="Colormap",
                           items=[(k, v[0], "") for k, v in nodes.COLORMAPS.items()])

    def execute(self, context):
        nodes.set_colormap(nodes.ensure_colored_material(), self.colormap)
        return {'FINISHED'}


class NA_OT_delete(_ArrowOp):
    """Delete this arrow set (the membrane stays)"""
    bl_idname = "normal_arrows.delete"
    bl_label = "Delete Arrows"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        obj = context.active_object
        source = nodes.source_of(obj)
        colls = [c for c in obj.users_collection if c.name.startswith(COLLECTION_PREFIXES)]
        mesh = obj.data
        bpy.data.objects.remove(obj)
        if mesh.users == 0:
            bpy.data.meshes.remove(mesh)
        for coll in colls:
            if not coll.objects and not coll.children:
                bpy.data.collections.remove(coll)
        if source is not None:
            _select_only(context, source)
        return {'FINISHED'}


CLASSES = (
    NA_OT_add,
    NA_OT_load_forces,
    NA_OT_update,
    NA_OT_fit_range,
    NA_OT_colormap,
    NA_OT_delete,
)
