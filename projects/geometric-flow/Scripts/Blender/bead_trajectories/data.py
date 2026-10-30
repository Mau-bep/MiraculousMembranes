"""Reading a results folder of the simulation: the bead files, the config and the saved mesh steps."""

import bisect
import json
import os
import re

import numpy as np

MESH_RE = re.compile(r"^membrane_(\d+)\.obj$")


def mesh_steps(folder):
    """Steps of the saved meshes membrane_<step>.obj, sorted (the order of an OBJ sequence)."""
    steps = []
    for name in os.listdir(folder):
        m = MESH_RE.match(name)
        if m:
            steps.append(int(m.group(1)))
    return sorted(steps)


def bead_files(folder):
    """The paths of Bead_0_data.txt, Bead_1_data.txt, ... that exist, in order."""
    paths = []
    while os.path.exists(os.path.join(folder, "Bead_{}_data.txt".format(len(paths)))):
        paths.append(os.path.join(folder, "Bead_{}_data.txt".format(len(paths))))
    return paths


def read_bead_file(path):
    """The rows `step x y z Fx Fy Fz` of one bead file as (steps, values), steps increasing.

    A step that appears twice (a restarted run) keeps its last row; lines that start with # and short lines are skipped.
    """
    rows = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            w = line.split()
            if len(w) < 4:
                continue
            try:
                v = [float(t) for t in w[:7]]
            except ValueError:
                continue
            v += [0.0] * (7 - len(v))
            rows[int(v[0])] = v[1:]
    steps = sorted(rows)
    return np.array(steps, dtype=np.int64), np.array([rows[s] for s in steps], dtype=float).reshape(-1, 6)


def radii(folder, n_beads):
    """The bead radii of Input_file.json (None where the file or the entry is missing)."""
    path = os.path.join(folder, "Input_file.json")
    out = [None] * n_beads
    if os.path.exists(path):
        try:
            with open(path) as f:
                beads = json.load(f).get("Beads", [])
            for i, b in enumerate(beads[:n_beads]):
                if "radius" in b:
                    out[i] = float(b["radius"])
        except (OSError, ValueError):
            pass
    return out


def sample(steps, values, wanted):
    """The row of every wanted step: the last row at or before it (the first row if there is none), like resuming a run does."""
    idx = np.clip(np.searchsorted(steps, wanted, side="right") - 1, 0, len(steps) - 1)
    return values[idx]


def load(folder):
    """(steps, positions, forces, radii) of a results folder: steps (F,), positions and forces (N, F, 3), radii (N,) with None for unknown.

    The steps are those of the saved meshes (so the beads sit where the membrane is on every frame of its OBJ sequence); a folder
    without meshes takes the steps of the first bead file.
    """
    folder = os.path.abspath(bpy_path(folder))
    paths = bead_files(folder)
    if not paths:
        raise ValueError("No Bead_0_data.txt in " + folder)
    data = [read_bead_file(p) for p in paths]
    if any(len(s) == 0 for s, _ in data):
        raise ValueError("A bead file in {} has no data rows".format(folder))
    steps = np.array(mesh_steps(folder), dtype=np.int64)
    if len(steps) == 0:
        steps = data[0][0]
    pos = np.stack([sample(s, v, steps)[:, 0:3] for s, v in data])
    force = np.stack([sample(s, v, steps)[:, 3:6] for s, v in data])
    return steps, pos, force, radii(folder, len(paths))


def bpy_path(path):
    """Resolve a Blender path (//relative, ~) to a normal one."""
    try:
        import bpy
        return bpy.path.abspath(path)
    except ImportError:
        return os.path.expanduser(path)
