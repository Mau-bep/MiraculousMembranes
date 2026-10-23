"""Config for pulling N membrane tubes out of a sphere (a many-tube Tritube).

Each tube has a bead that starts inside the vesicle and a frozen anchor outside,
on the same radial line, joined by a Harmonic spring (like Tritube.json, but with
a rest length). The spring pulls the bead outwards and the bead pushes the
membrane out through its LJ (WCA) repulsion. The anchor sits --anchor_gap beyond
the point where the bead should stop (rescale + --length from the centre) and the
spring rest length is --anchor_gap, so the tube ends where you ask it to, up to
the spring stretch (tube force / --k_spring, about 0.1 with the defaults).

Layouts:
  ring    N tubes evenly spaced on a ring at polar angle --angle from the +z pole
          (90 = the equator, like Tritube.json).
  random  N tubes in random directions that are at least --min_sep apart along the
          surface, so no two beads pull on the same place (--seed to repeat).

    python Create_config_tubes.py -N 6 --layout ring
    python Create_config_tubes.py -N 8 --layout random --seed 3

Run from Scripts/ (like the Create_subjob_*.py scripts). The file goes to
../Config_files/ and its path is printed.
"""
import argparse
import itertools
import json
import math
import os
import random

from jinja2 import Environment, FileSystemLoader


def ring_directions(n, angle_deg, azimuth0_deg=90.0):
    theta = math.radians(angle_deg)
    return [[math.sin(theta) * math.cos(math.radians(azimuth0_deg + 360.0 * k / n)),
             math.sin(theta) * math.sin(math.radians(azimuth0_deg + 360.0 * k / n)),
             math.cos(theta)] for k in range(n)]


def random_directions(n, min_angle, seed, tries=2000, restarts=500):
    """N unit vectors, every pair at least min_angle (radians) apart (rejection sampling)."""
    rng = random.Random(seed)
    cos_max = math.cos(min_angle)
    for _ in range(restarts):
        dirs = []
        for _ in range(n):
            for _ in range(tries):
                v = [rng.gauss(0, 1) for _ in range(3)]
                norm = math.sqrt(sum(c * c for c in v))
                v = [c / norm for c in v]
                if all(sum(a * b for a, b in zip(v, d)) <= cos_max for d in dirs):
                    dirs.append(v)
                    break
            else:
                break  # stuck: start over
        if len(dirs) == n:
            return dirs
    raise SystemExit("Could not place {} tubes at least {:.2f} rad ({:.1f} deg) apart: use fewer tubes, "
                     "a smaller --min_sep or a larger --rescale".format(n, min_angle, math.degrees(min_angle)))


def create_config(n, layout="ring", angle=90.0, azimuth0=90.0, seed=0, rescale=2.0, radius=0.25,
                  length=1.3, start_depth=0.5, anchor_gap=1.0, anchor_radius=0.4, k_spring=200.0,
                  min_sep=None, steps=3000, KA=16.0, KB=1.0, KV=10000.0, bending_energy="Bending",
                  inter_str=10.0, size_max=0.2, save_interval=10, first_dir="../Results/Tube_cluster/",
                  subfolder="1/", out_dir="../Config_files/"):
    if min_sep is None:
        min_sep = 8 * radius
    if start_depth <= 0 or start_depth >= rescale:
        raise SystemExit("--start_depth must be between 0 and --rescale (the bead starts inside the vesicle)")
    if 2 * math.pi * radius < 6 * size_max:
        print("Warning: bead circumference {:.2f} is less than 6 remesher edges (size_max {}): the tube will be "
              "faceted, lower --size_max".format(2 * math.pi * radius, size_max))

    if layout == "ring":
        dirs = ring_directions(n, angle, azimuth0)
    else:
        dirs = random_directions(n, min_sep / rescale, seed)

    # Surface distance between the closest tube bases
    if n > 1:
        dmin = min(rescale * math.acos(max(-1.0, min(1.0, sum(a * b for a, b in zip(u, v)))))
                   for u, v in itertools.combinations(dirs, 2))
        if dmin < min_sep:
            print("Warning: the closest tubes start {:.2f} apart along the surface, less than --min_sep {:.2f}"
                  .format(dmin, min_sep))

    beads = []
    for u in dirs:
        beads.append({"pos": [round((rescale - start_depth) * c, 6) for c in u],
                      "anchor": [round((rescale + length + anchor_gap) * c, 6) for c in u]})

    env = Environment(loader=FileSystemLoader(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../Templates/")))
    text = env.get_template("Tube_cluster.txt").render(
        beads=beads, radius=radius, anchor_radius=anchor_radius, anchor_gap=anchor_gap, k_spring=k_spring,
        inter_str=inter_str, KA=KA, KB=KB, KV=KV or None, bending_energy=bending_energy, rescale=rescale,
        size_max=size_max, timesteps=steps, save_interval=save_interval, first_dir=first_dir, subfolder=subfolder)
    data = json.loads(text)  # fails here if the template renders invalid JSON

    os.makedirs(out_dir, exist_ok=True)
    name = "Tubes_N_{}_{}_L_{}_KA_{}_KB_{}".format(n, layout, length, KA, KB)
    name += "_angle_{}".format(angle) if layout == "ring" else "_seed_{}".format(seed)
    path = os.path.join(out_dir, name + ".json")
    with open(path, "w") as f:
        json.dump(data, f, indent=4)
    return path


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-N", type=int, required=True, help="number of tubes")
    p.add_argument("--layout", choices=["ring", "random"], default="ring")
    p.add_argument("--angle", type=float, default=90.0, help="ring: polar angle from the +z pole (degrees)")
    p.add_argument("--azimuth0", type=float, default=90.0, help="ring: azimuth of the first tube (degrees)")
    p.add_argument("--seed", type=int, default=0, help="random: seed of the positions")
    p.add_argument("--min_sep", type=float, default=None,
                   help="min distance along the surface between tubes (default 8 bead radii)")
    p.add_argument("--rescale", type=float, default=2.0, help="sphere radius")
    p.add_argument("--radius", type=float, default=0.25, help="bead radius (sets the tube radius)")
    p.add_argument("--length", type=float, default=1.3, help="how far past the sphere surface each bead is pulled")
    p.add_argument("--start_depth", type=float, default=0.5, help="how far inside the surface the beads start")
    p.add_argument("--anchor_gap", type=float, default=1.0, help="anchor distance beyond the final bead position "
                   "(also the spring rest length)")
    p.add_argument("--anchor_radius", type=float, default=0.4, help="radius of the frozen anchor beads")
    p.add_argument("--k_spring", type=float, default=200.0, help="spring constant of the bead-anchor bond")
    p.add_argument("--steps", type=int, default=3000, help="number of timesteps")
    p.add_argument("--KA", type=float, default=16.0, help="surface tension")
    p.add_argument("--KB", type=float, default=1.0, help="bending constant")
    p.add_argument("--KV", type=float, default=10000.0, help="Volume_constraint constant (0 = no constraint)")
    p.add_argument("--bending_energy", choices=["Bending", "Bending_tan"], default="Bending")
    p.add_argument("--inter_str", type=float, default=10.0, help="bead-membrane LJ strength (shift = strength)")
    p.add_argument("--size_max", type=float, default=0.2, help="remesher maximum edge length")
    p.add_argument("--save_interval", type=int, default=10)
    p.add_argument("--first_dir", default="../Results/Tube_cluster/")
    p.add_argument("--subfolder", default="1/")
    p.add_argument("--out_dir", default="../Config_files/")
    a = vars(p.parse_args())
    a["n"] = a.pop("N")
    print(create_config(**a))
