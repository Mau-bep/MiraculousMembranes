"""Config for N beads on a ring around the north pole of a sphere.

The beads are evenly spaced in azimuth, all at the same polar angle from the
north pole (+z), and touch the sphere from outside. Every bead has the same
membrane adhesion ("Adhesion" = coverage) and every pair of beads is bonded with
a Shifted_LJ repulsion (bonds_constants [eps, sigma]).

    python Create_config_cluster.py -N 5 --angle 35 --KA 0.5 --KB 10 \
        --inter_str 25 --eps 10

Run from Scripts/ (like the Create_subjob_*.py scripts). The file goes to
../Config_files/ and its path is printed.
"""
import argparse
import itertools
import json
import math
import os

from jinja2 import Environment, FileSystemLoader


def ring_positions(n, angle_deg, distance, azimuth0_deg=90.0):
    """N points at colatitude angle_deg, evenly spaced in azimuth."""
    theta = math.radians(angle_deg)
    positions = []
    for k in range(n):
        phi = math.radians(azimuth0_deg + 360.0 * k / n)
        positions.append([distance * math.sin(theta) * math.cos(phi),
                          distance * math.sin(theta) * math.sin(phi),
                          distance * math.cos(theta)])
    return positions


def create_config(n, angle, KA, KB, inter_str, eps, sigma=1.0, radius=1.0, rescale=7.0,
                  KV=20000.0, timesteps=20000, save_interval=500, switch_time=2000,
                  azimuth0=90.0, first_dir="../Results/Cluster_beads/", subfolder="1/",
                  out_dir="../Config_files/"):
    # The mesh sphere has radius ~1 before rescale, so the bead centre touches it at rescale + radius
    positions = ring_positions(n, angle, rescale + radius, azimuth0)
    beads = [{"pos": [round(c, 6) for c in p], "partners": [j for j in range(n) if j != k]}
             for k, p in enumerate(positions)]

    if n > 1:
        dmin = min(math.dist(a, b) for a, b in itertools.combinations(positions, 2))
        if dmin < 2 * radius:
            print("Warning: the closest beads are {:.3f} apart, less than 2*radius = {:.3f} (they overlap)"
                  .format(dmin, 2 * radius))
        elif dmin < 2 ** (1 / 6) * sigma:
            print("Note: the closest beads are {:.3f} apart, inside the bond cutoff {:.3f}"
                  .format(dmin, 2 ** (1 / 6) * sigma))
    if not 0 <= angle <= 90:
        print("Warning: angle {} is not in the northern hemisphere (0-90 degrees)".format(angle))

    env = Environment(loader=FileSystemLoader(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../Templates/")))
    text = env.get_template("Wrapping_cluster.txt").render(
        beads=beads, radius=radius, inter_str=inter_str, eps=eps, sigma=sigma, KA=KA, KB=KB, KV=KV,
        rescale=rescale, timesteps=timesteps, save_interval=save_interval, switch_time=switch_time,
        first_dir=first_dir, subfolder=subfolder)
    data = json.loads(text)  # fails here if the template renders invalid JSON

    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, "Cluster_N_{}_angle_{}_KA_{}_KB_{}_strg_{}_eps_{}.json"
                        .format(n, angle, KA, KB, inter_str, eps))
    with open(path, "w") as f:
        json.dump(data, f, indent=4)
    return path


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-N", type=int, required=True, help="number of beads")
    p.add_argument("--angle", type=float, required=True, help="polar angle of the beads from the north pole (degrees)")
    p.add_argument("--KA", type=float, required=True, help="surface tension")
    p.add_argument("--KB", type=float, required=True, help="Bending_tan constant")
    p.add_argument("--inter_str", type=float, required=True, help="bead-membrane adhesion strength (all beads)")
    p.add_argument("--eps", type=float, required=True, help="bead-bead Shifted_LJ epsilon (all pairs)")
    p.add_argument("--sigma", type=float, default=1.0, help="bead-bead Shifted_LJ sigma (default 1)")
    p.add_argument("--radius", type=float, default=1.0, help="bead radius")
    p.add_argument("--rescale", type=float, default=7.0, help="sphere scale")
    p.add_argument("--KV", type=float, default=20000.0, help="Volume_constraint constant")
    p.add_argument("--timesteps", type=int, default=20000)
    p.add_argument("--save_interval", type=int, default=500)
    p.add_argument("--switch_time", type=int, default=2000, help="step at which Gradient_descent switches to BFGS")
    p.add_argument("--azimuth0", type=float, default=90.0, help="azimuth of the first bead (degrees)")
    p.add_argument("--first_dir", default="../Results/Cluster_beads/")
    p.add_argument("--subfolder", default="1/")
    p.add_argument("--out_dir", default="../Config_files/")
    a = vars(p.parse_args())
    a["n"] = a.pop("N")
    print(create_config(**a))
