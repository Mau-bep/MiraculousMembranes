"""Config for N beads wrapping on the big planar membrane (the planar sister of Create_config_cluster.py and Create_config_tubes.py).

The membrane is the disk of Big_planar_mem.obj (radius 10 before the rescale, in the plane x = 0, clamped edge) and every bead starts with its centre one radius above
it, at the same height x = radius, like the one-bead runs of Wrapping_planar.txt. All beads have the same adhesion strength (the "Adhesion" interaction, soft shell) and
every pair of beads is bonded with a Shifted_LJ bond (bonds_constants [eps, sigma], sigma = 2 radius by default: the beads do not overlap). By default the bond is the
repulsive core of Lennard-Jones, cut at its minimum 2^(1/6) sigma (2.245 radii): beads further apart do not feel each other. --attract keeps the attractive tail of the
potential out to --bond_cutoff sigma (2.5 by default, bonds_constants [eps, sigma, rc/sigma]): the beads attract each other with a well of depth about eps at 2^(1/6) sigma.
--bond_potential frenkel uses the Frenkel potential instead (bond "Frenkel", bonds_constants [eps, sigma, rc], rc = --bond_cutoff sigma): V = eps alpha (sigma^2/r^2 - 1)
(rc^2/r^2 - 1)^2 is repulsive inside sigma (where it crosses zero), has a well of depth exactly eps at r^2 = 3 sigma^2 rc^2 / (sigma^2 + 2 rc^2) (2.357 radii for
sigma = 2 radii and rc = 2.5 sigma) and goes smoothly to zero, value and slope, at rc. It is always attractive, so --attract is not needed; hex, columns and pyramid
default to its minimum as the neighbour distance.

Layouts (--spacing is the distance between the centres of neighbours, default 3 radii, 2^(1/6) sigma = the minimum of the bead-bead potential for hex):
  line    N beads in a row, evenly spaced and centred on the membrane (default). --line_angle turns the row in the plane (0 = along y).
  ring    N beads evenly spaced on a circle; the radius makes the neighbours --spacing apart unless --ring_radius is given.
  hex     the N sites of a hexagonal packing closest to the centre (7 = a full hexagon, 19 = two shells), neighbours --spacing apart: by default the minimum
          of the bead-bead potential, 2^(1/6) sigma = 2.245 radii, so the beads are in contact range and the adhesion shells overlap.
  columns four columns of beads standing on the membrane, on the corners of a square of side --spacing (default 2^(1/6) sigma, like hex), stacked along the
          normal --layer_spacing apart (default --spacing) and filled layer by layer: N = 8 is a cube of 2 x 2 x 2 beads, 12 a 2 x 2 x 3 tower, 4 a single layer.
          The bottom layer starts one radius above the membrane like every other bead; the upper layers hang on the bead-bead bonds (use --attract or --bond_potential frenkel, without it
          the beads only repel and the upper ones feel nothing). --azimuth0 turns the square in the plane.
  pyramid close-packed pyramid of beads on the membrane, tip up: --base triangle (default, a tetrahedron: 1, 4, 10, 20, 35, 56 beads, layers of 1, 3, 6, 10 ...) or
          --base square (1, 5, 14, 30, 55, 91 beads, layers of 1, 4, 9, 16 ...). N must be one of those numbers so that the pyramid is full. Neighbours are --spacing
          apart (default 2^(1/6) sigma, like hex), every bead sits in the hollow of the layer below, which makes the layers --spacing * sqrt(2/3) (triangle) or
          --spacing / sqrt(2) (square) apart unless --layer_spacing is given. The base layer starts one radius above the membrane. --azimuth0 turns it in the plane.
          --hollow keeps only the surface beads (base layer and the sloped faces, none inside): then N is the number of beads placed, triangle 1, 4, 10, 20, 34, 52, 74,
          square 1, 5, 14, 29, 50, 77 (the first pyramid with an interior is the 5-layer triangular one and the 4-layer square one).
          Use --attract or --bond_potential frenkel: the upper layers hang on the bead-bead bonds.
  random  N beads at random in the annulus --r_inner .. --r_outer (default: a disc about as big as the packing needs), no two closer than --min_sep
          (default --spacing); --seed repeats the positions.

    python Create_config_planar_beads.py -N 5 --KA 0 --KB 1 --inter_str 1.5
    python Create_config_planar_beads.py -N 7 --layout hex --KA 0 --KB 1 --inter_str 1.5 --bead_normal --attract
    python Create_config_planar_beads.py -N 7 --layout hex --KA 0 --KB 1 --inter_str 1.5 --attract --bond_cutoff 2.0 --eps 2
    python Create_config_planar_beads.py -N 10 --layout pyramid --bond_potential frenkel --eps 3 --KA 0 --KB 1 --inter_str 1.5
    python Create_config_planar_beads.py -N 8 --layout columns --KA 0 --KB 1 --inter_str 1.5 --bead_normal --attract
    python Create_config_planar_beads.py -N 10 --layout pyramid --KA 0 --KB 1 --inter_str 1.5 --bead_normal --attract
    python Create_config_planar_beads.py -N 14 --layout pyramid --base square --KA 0 --KB 1 --inter_str 1.5 --attract
    python Create_config_planar_beads.py -N 52 --layout pyramid --hollow --KA 0 --KB 1 --inter_str 1.5 --attract
    python Create_config_planar_beads.py -N 8 --layout random --r_inner 3 --r_outer 6 --seed 2 --KA 0 --KB 1 --inter_str 1.5

The disk is rescaled (--rescale, default: automatic) so that the outermost bead keeps --margin radii (default 9, the one-bead runs of rescale 1) to the clamped edge.
Run from Scripts/ (like the Create_subjob_*.py scripts). The config goes to ../Config_files/ and its path is printed; --subjob also writes the sbatch file
../Subjobs/<name>_subjob.
"""
import argparse
import itertools
import json
import math
import os
import random

from jinja2 import Environment, FileSystemLoader


# Every layout returns the in-plane coordinates (y, z) of the bead centres; the membrane is the plane x = 0 and its centre is the origin.

def line_positions(n, spacing, line_angle_deg=0.0):
    """N points on a line through the origin, evenly spaced, centred."""
    a = math.radians(line_angle_deg)
    return [((k - (n - 1) / 2.0) * spacing * math.cos(a), (k - (n - 1) / 2.0) * spacing * math.sin(a)) for k in range(n)]


def ring_positions(n, spacing, ring_radius=None, azimuth0_deg=0.0):
    """N points evenly spaced on a circle (radius: neighbours spacing apart, a single bead sits in the middle)."""
    if n == 1:
        return [(0.0, 0.0)]
    R = ring_radius if ring_radius is not None else spacing / (2.0 * math.sin(math.pi / n))
    return [(R * math.cos(math.radians(azimuth0_deg) + 2.0 * math.pi * k / n),
             R * math.sin(math.radians(azimuth0_deg) + 2.0 * math.pi * k / n)) for k in range(n)]


def hex_positions(n, spacing):
    """The N sites of the hexagonal lattice closest to a lattice site (rows along y), then centred on their centroid."""
    sites = []
    m = int(math.ceil(math.sqrt(n))) + 2
    for i in range(-m, m + 1):
        for j in range(-m, m + 1):
            y = spacing * (i + 0.5 * j)
            z = spacing * (math.sqrt(3.0) / 2.0) * j
            sites.append((y, z))
    # nearest first; the angle only breaks ties between the sites of one shell so that partial shells fill in order
    sites.sort(key=lambda p: (round(math.hypot(*p), 6), math.atan2(p[1], p[0]) % (2.0 * math.pi)))
    pts = sites[:n]
    cy = sum(p[0] for p in pts) / n
    cz = sum(p[1] for p in pts) / n
    return [(y - cy, z - cz) for y, z in pts]


def column_positions(n, spacing, layer_spacing, rotation_deg=0.0):
    """The in-plane (y, z) and the height above the first layer of N beads in 4 columns on the corners of a square, layer by layer."""
    a = math.radians(rotation_deg)
    corners = [(-spacing / 2.0, -spacing / 2.0), (spacing / 2.0, -spacing / 2.0), (spacing / 2.0, spacing / 2.0), (-spacing / 2.0, spacing / 2.0)]
    yz = [(y * math.cos(a) - z * math.sin(a), y * math.sin(a) + z * math.cos(a)) for y, z in corners]
    return [yz[k % 4] for k in range(n)], [(k // 4) * layer_spacing for k in range(n)]


def pyramid_sites(layers, base, hollow=False):
    """The lattice sites (layer j, a, b) of a pyramid of the given number of layers, base layer j = 0 first. With hollow only the surface sites are kept:
    the base layer and the sloped faces (triangle: a, b or the third barycentric coordinate zero; square: the border of the layer)."""
    sites = []
    for j in range(layers):
        m = layers - j
        for b in range(m):
            for a in range(m - b if base == "triangle" else m):
                if hollow and j > 0:
                    if base == "triangle":
                        interior = a > 0 and b > 0 and m - 1 - a - b > 0
                    else:
                        interior = 0 < a < m - 1 and 0 < b < m - 1
                    if interior:
                        continue
                sites.append((j, a, b))
    return sites


def pyramid_numbers(base, hollow=False, count=8):
    """The first pyramid numbers: tetrahedral n(n+1)(n+2)/6 for a triangular base, square pyramidal n(n+1)(2n+1)/6 for a square one; hollow: the beads of
    the surface only (the full number minus the interior pyramid of n - 4 layers for a triangular base, n - 3 for a square one)."""
    return [len(pyramid_sites(n, base, hollow)) for n in range(1, count + 1)]


def pyramid_positions(n, base, spacing, layer_spacing, rotation_deg=0.0, hollow=False):
    """The in-plane (y, z) and the height above the base layer of the N beads of a full pyramid (base layer first, tip last).

    triangle: layer j (0 = base) is a triangular patch of the hexagonal lattice with m = L - j beads on a side; square: a square grid of m x m beads.
    Every layer is centred on the common axis, which puts each bead in the hollow of the layer below.
    """
    numbers = pyramid_numbers(base, hollow, 24)
    if n not in numbers:
        below = [k for k in numbers if k < n]
        above = [k for k in numbers if k > n]
        raise SystemExit("N = {} is not a {}{} pyramid number: the pyramid would not be full. The numbers are {}; the closest are {} and {}"
                         .format(n, "hollow " if hollow else "", base, ", ".join(str(k) for k in numbers[:8]),
                                 below[-1] if below else "-", above[0] if above else "-"))
    layers = numbers.index(n) + 1
    pts, heights = [], []
    for j, a, b in pyramid_sites(layers, base, hollow):
        m = layers - j
        if base == "triangle":
            y = spacing * (a + 0.5 * b) - spacing * (m - 1) / 2.0
            z = spacing * (math.sqrt(3.0) / 2.0) * b - spacing * (m - 1) * math.sqrt(3.0) / 6.0  # centred on the centroid of the triangle of beads
        else:
            y = spacing * a - spacing * (m - 1) / 2.0
            z = spacing * b - spacing * (m - 1) / 2.0
        pts.append((y, z))
        heights.append(j * layer_spacing)
    ang = math.radians(rotation_deg)
    pts = [(y * math.cos(ang) - z * math.sin(ang), y * math.sin(ang) + z * math.cos(ang)) for y, z in pts]
    return pts, heights


def random_positions(n, min_sep, r_inner, r_outer, seed, tries=3000, restarts=300):
    """N points uniform in the annulus r_inner..r_outer, every pair at least min_sep apart (rejection sampling)."""
    rng = random.Random(seed)
    for _ in range(restarts):
        pts = []
        for _ in range(n):
            for _ in range(tries):
                rho = math.sqrt(rng.uniform(r_inner ** 2, r_outer ** 2))
                phi = rng.uniform(0.0, 2.0 * math.pi)
                p = (rho * math.cos(phi), rho * math.sin(phi))
                if all(math.dist(p, q) >= min_sep for q in pts):
                    pts.append(p)
                    break
            else:
                break  # stuck: start over
        if len(pts) == n:
            return pts
    raise SystemExit("Could not place {} beads at least {:g} apart in the annulus {:g}..{:g}: use fewer beads, a smaller --min_sep "
                     "or a bigger annulus".format(n, min_sep, r_inner, r_outer))


def bond_minimum(potential, sigma, cutoff):
    """Distance of the minimum of the bead-bead potential: 2^(1/6) sigma for LJ, sqrt(3 sigma^2 rc^2 / (sigma^2 + 2 rc^2)) for Frenkel (rc = cutoff * sigma)."""
    if potential == "frenkel":
        rc = cutoff * sigma
        return math.sqrt(3.0 * sigma ** 2 * rc ** 2 / (sigma ** 2 + 2.0 * rc ** 2))
    return 2.0 ** (1.0 / 6.0) * sigma


def create_config(n, KA, KB, inter_str, layout="line", spacing=None, radius=1.0, eps=None, sigma=None,
                  attract=False, bond_cutoff=2.5, bond_potential="lj", line_angle=0.0,
                  ring_radius=None, azimuth0=0.0, layer_spacing=None, base="triangle", hollow=False, r_inner=0.0, r_outer=None, min_sep=None, seed=0, margin=9.0,
                  rescale=None, tension="surface", bead_normal=False, shell_width=None, shell_power=None,
                  shell_profile=None, shell_alpha=None, size_min=0.005, size_max=0.65, refine_angle=0.4,
                  saved_states=15, switch_bfgs=2000, switch_normal=40000, freeze_beads=80000, timesteps=200000,
                  save_interval=500, first_dir="../Results/Planar_beads/", subfolder="1/", out_dir="../Config_files/"):
    if sigma is None:
        sigma = 2.0 * radius
    frenkel = bond_potential == "frenkel"
    attract = attract or frenkel  # the Frenkel potential always has its well
    r_min = bond_minimum(bond_potential, sigma, bond_cutoff)  # where the bead-bead potential is lowest
    if spacing is None:
        # hexagonal packing and columns: neighbours at the minimum of the bead-bead potential; the others 3 radii apart
        spacing = r_min if layout in ("hex", "columns", "pyramid") else 3.0 * radius
    if layer_spacing is None:
        # a pyramid is close packed: the layers are as far apart as the hollows of the layer below allow
        layer_spacing = spacing * (math.sqrt(2.0 / 3.0) if base == "triangle" else math.sqrt(0.5)) if layout == "pyramid" else spacing
    if eps is None:
        eps = 1.0 if attract else 10.0  # a repulsive core is stiff, a well of depth 10 would dominate the adhesion
    if frenkel and bond_cutoff <= 1.0:
        raise SystemExit("--bond_cutoff {:g}: the Frenkel cutoff rc = bond_cutoff * sigma must be larger than sigma".format(bond_cutoff))
    if attract and not frenkel and bond_cutoff < 2.0 ** (1.0 / 6.0):
        raise SystemExit("--bond_cutoff {:g} is below the minimum of the potential, 2^(1/6) = 1.1225 (in units of sigma)".format(bond_cutoff))
    if spacing < 2.0 * radius and n > 1 and layout != "random":
        print("Warning: --spacing {:g} is less than 2*radius = {:g}: neighbouring beads overlap".format(spacing, 2.0 * radius))

    heights = [0.0] * n  # height of each bead above the first layer (only the columns stack)
    if layout == "line":
        yz = line_positions(n, spacing, line_angle)
    elif layout == "ring":
        yz = ring_positions(n, spacing, ring_radius, azimuth0)
    elif layout == "hex":
        yz = hex_positions(n, spacing)
    elif layout == "columns":
        yz, heights = column_positions(n, spacing, layer_spacing, azimuth0)
        if layer_spacing < 2.0 * radius:
            print("Warning: --layer_spacing {:g} is less than 2*radius = {:g}: the beads of a column overlap".format(layer_spacing, 2.0 * radius))
        if n % 4:
            print("Note: {} beads do not fill the {} layers of 4: the top layer has {}".format(n, -(-n // 4), n % 4))
    elif layout == "pyramid":
        yz, heights = pyramid_positions(n, base, spacing, layer_spacing, azimuth0, hollow)
        layers = pyramid_numbers(base, hollow, 24).index(n) + 1
        print("{}{} pyramid of {} layers, base layer of {} beads{}".format("hollow " if hollow else "", base, layers, sum(h == 0.0 for h in heights),
                                                                        ", {} interior beads left out".format(pyramid_numbers(base, False, layers)[-1] - n) if hollow else ""))
    elif layout == "random":
        sep = min_sep if min_sep is not None else spacing
        if r_outer is None:
            # a disc whose area is about 1.5 hexagonal-packing cells per bead (enough room for rejection sampling)
            r_outer = max(r_inner + sep, sep * math.sqrt(1.5 * 0.866 * n / math.pi) + sep / 2.0)
        if r_outer <= r_inner:
            raise SystemExit("--r_outer ({:g}) must be larger than --r_inner ({:g})".format(r_outer, r_inner))
        yz = random_positions(n, sep, r_inner, r_outer, seed)
    else:
        raise SystemExit("unknown layout " + layout)

    xyz = [(h, y, z) for h, (y, z) in zip(heights, yz)]
    if n > 1:
        dmin = min(math.dist(a, b) for a, b in itertools.combinations(xyz, 2))
        if dmin < 2.0 * radius:
            print("Warning: the closest beads are {:.3f} apart, less than 2*radius = {:.3f} (they overlap)".format(dmin, 2.0 * radius))
        elif dmin < r_min - 1e-9:
            print("Note: the closest beads are {:.3f} apart, inside the minimum of the bead-bead potential {:.3f}: they repel"
                  .format(dmin, r_min))
        if attract:
            reach = bond_cutoff * sigma
            print("Bead-bead {} attraction: well of depth {}{:g} at {:.3f}, cutoff {:.3f}; {} of {} pairs are within the cutoff"
                  .format("Frenkel" if frenkel else "LJ", "" if frenkel else "~", eps, r_min, reach, sum(math.dist(a, b) < reach for a, b in itertools.combinations(xyz, 2)),
                          n * (n - 1) // 2))
        elif dmin > r_min:
            print("Note: the closest beads are {:.3f} apart, beyond the repulsion cutoff {:.3f}: no bead-bead force (--attract for the attractive tail)"
                  .format(dmin, r_min))

    # The disk has radius 10 * rescale: the outermost bead keeps `margin` radii to the clamped edge
    extent = max(math.hypot(*p) for p in yz) + radius
    if rescale is None:
        rescale = max(1.0, math.ceil((extent + margin * radius) / 10.0 * 2.0) / 2.0)
    elif 10.0 * rescale < extent + 3.0 * radius:
        print("Warning: --rescale {:g} gives a disk of radius {:g}, the beads reach {:g}: they are on the edge".format(rescale, 10.0 * rescale, extent))
    edge_gap = (10.0 * rescale - extent) / radius

    beads = [{"pos": [round(radius + h, 6), round(y, 6), round(z, 6)], "partners": [j for j in range(n) if j != k]}
             for k, (h, y, z) in enumerate(xyz)]

    env = Environment(loader=FileSystemLoader(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../Templates/")),
                      trim_blocks=True, lstrip_blocks=True)
    # The optional shell settings are left out of the render (the code defaults apply) when not given
    shell = {k: v for k, v in (("shell_width", shell_width), ("shell_power", shell_power), ("shell_profile", shell_profile),
                               ("shell_alpha", shell_alpha)) if v is not None}
    text = env.get_template("Wrapping_planar_beads.txt").render(
        beads=beads, radius=radius, inter_str=inter_str, bond_type="Frenkel" if frenkel else "Shifted_LJ",
        bond=[eps, sigma, bond_cutoff * sigma] if frenkel else [eps, sigma] + ([bond_cutoff] if attract else []), KA=KA, KB=KB, excess=(tension == "excess"),
        bead_normal=bead_normal, rescale=rescale, size_min=size_min, size_max=size_max, refine_angle=refine_angle,
        saved_states=saved_states, switch_bfgs=switch_bfgs, switch_normal=switch_normal, freeze_beads=freeze_beads,
        timesteps=timesteps, save_interval=save_interval, first_dir=first_dir, subfolder=subfolder, **shell)
    data = json.loads(text)  # fails here if the template renders invalid JSON

    os.makedirs(out_dir, exist_ok=True)
    name = "Planar_beads_N_{}_{}_d_{:g}_r_{}_KA_{}_KB_{}_strg_{}".format(n, layout, spacing, radius, KA, KB, inter_str)
    if layout == "columns" and layer_spacing != spacing:
        name += "_dz_{:g}".format(layer_spacing)
    if layout == "pyramid":
        name = name.replace("_pyramid_", "_pyramid_{}{}_".format(base, "_hollow" if hollow else ""), 1)
    if frenkel:
        name += "_frenkel_eps_{}_rc_{}".format(eps, bond_cutoff)
    elif attract:
        name += "_att_eps_{}_rc_{}".format(eps, bond_cutoff)
    if layout == "random":
        name += "_seed_{}".format(seed)
    path = os.path.join(out_dir, name + ".json")
    with open(path, "w") as f:
        json.dump(data, f, indent=4)
    print("{}: {} beads in a {} layout, disk radius {:g} (rescale {:g}), {:.1f} radii between the outermost bead and the edge"
          .format(name, n, layout, 10.0 * rescale, rescale, edge_gap))
    return path


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-N", type=int, required=True, help="number of beads")
    p.add_argument("--KA", type=float, required=True, help="surface tension")
    p.add_argument("--KB", type=float, required=True, help="Bending_tan constant (kappa = KB/2)")
    p.add_argument("--inter_str", type=float, required=True, help="bead-membrane adhesion strength (all beads)")
    p.add_argument("--layout", choices=["line", "ring", "hex", "columns", "pyramid", "random"], default="line")
    p.add_argument("--spacing", type=float, default=None,
                   help="distance between neighbouring bead centres (default 3 radii; hex, columns and pyramid: the minimum of the bead-bead potential, 2^(1/6) sigma for LJ)")
    p.add_argument("--radius", type=float, default=1.0, help="bead radius")
    p.add_argument("--eps", type=float, default=None,
                   help="bead-bead Shifted_LJ epsilon (all pairs; default 10 for the repulsive core, 1 with --attract or the Frenkel potential, the depth of the well)")
    p.add_argument("--attract", action="store_true",
                   help="keep the attractive tail of the bead-bead Lennard-Jones potential out to --bond_cutoff sigma (bonds_constants [eps, sigma, rc/sigma]); "
                        "without it the bond is repulsive, cut at 2^(1/6) sigma")
    p.add_argument("--bond_potential", choices=["lj", "frenkel"], default="lj",
                   help="bead-bead potential: lj = Shifted_LJ bond (repulsive core, or with --attract the Lennard-Jones tail); frenkel = Frenkel bond [eps, sigma, rc], "
                        "zero at sigma, a well of depth eps, smooth zero value and slope at rc = --bond_cutoff sigma (always attractive)")
    p.add_argument("--bond_cutoff", type=float, default=2.5,
                   help="--attract / frenkel: cutoff of the bead-bead potential in units of sigma (2.5 = the usual LJ cutoff, LJ at least 1.1225, Frenkel above 1)")
    p.add_argument("--sigma", type=float, default=None, help="bead-bead Shifted_LJ sigma (default 2 radii: the minimum of the potential is at 2^(1/6) sigma)")
    p.add_argument("--line_angle", type=float, default=0.0, help="line: angle of the row from the y axis, in the plane (degrees)")
    p.add_argument("--ring_radius", type=float, default=None, help="ring: radius of the circle (default: neighbours --spacing apart)")
    p.add_argument("--azimuth0", type=float, default=0.0, help="ring: angle of the first bead; columns: rotation of the square in the plane (degrees)")
    p.add_argument("--layer_spacing", type=float, default=None,
                   help="columns: distance between the layers of a column (default --spacing: a cube for N = 8); pyramid: distance between the layers "
                        "(default close packing, spacing * sqrt(2/3) for a triangular base, spacing / sqrt(2) for a square one)")
    p.add_argument("--base", choices=["triangle", "square"], default="triangle", help="pyramid: shape of the base (N = 1, 4, 10, 20, 35 ... or 1, 5, 14, 30, 55 ...)")
    p.add_argument("--hollow", action="store_true",
                   help="pyramid: only the beads of the surface (base layer and sloped faces), none inside; N is the number of beads placed "
                        "(triangle 1, 4, 10, 20, 34, 52, 74 ...; square 1, 5, 14, 29, 50, 77 ...)")
    p.add_argument("--r_inner", type=float, default=0.0, help="random: inner radius of the annulus")
    p.add_argument("--r_outer", type=float, default=None, help="random: outer radius of the annulus (default: a disc sized for N beads)")
    p.add_argument("--min_sep", type=float, default=None, help="random: smallest distance between two beads (default --spacing)")
    p.add_argument("--seed", type=int, default=0, help="random: seed of the positions")
    p.add_argument("--margin", type=float, default=9.0, help="radii between the outermost bead and the clamped edge (sets the automatic --rescale)")
    p.add_argument("--rescale", type=float, default=None, help="disk scale: the disk has radius 10*rescale (default: automatic from --margin, 1.0 for a single bead)")
    p.add_argument("--tension", choices=["surface", "excess"], default="surface",
                   help="surface: KA * area, excess: KA * (area - initial area) only while the area is above the initial one")
    p.add_argument("--bead_normal", action="store_true",
                   help="Constraint \"Normal\" on every bead: in BFGS-Normal (and its polish) the beads move only along the normal of the plane (x) and are variables of the "
                        "optimiser instead of staying frozen. It does NOT hold in the gradient descent and BFGS phases before: there the beads move freely, "
                        "and a column or a hexagon can slide apart sideways")
    p.add_argument("--shell_width", type=float, default=None, help="half width of the adhesion shell as a fraction of the bead radius (code default 0.25)")
    p.add_argument("--shell_power", type=float, default=None, help="exponent of the adhesion shell (code default 1)")
    p.add_argument("--shell_profile", choices=["cosine", "frenkel", "lj", "morse"], default=None, help="radial profile of the adhesion (code default cosine)")
    p.add_argument("--shell_alpha", type=float, default=None, help="range of the morse profile")
    p.add_argument("--size_min", type=float, default=0.005, help="remesher smallest edge")
    p.add_argument("--size_max", type=float, default=0.65, help="remesher largest edge")
    p.add_argument("--refine_angle", type=float, default=0.4, help="remesher refine_angle")
    p.add_argument("--saved_states", type=int, default=15, help="BFGS_saved_states")
    p.add_argument("--switch_bfgs", type=int, default=2000, help="step of the Gradient_descent -> BFGS switch")
    p.add_argument("--switch_normal", type=int, default=40000, help="step of the BFGS -> BFGS-Normal switch (the stopping block also switches on a plateau)")
    p.add_argument("--freeze_beads", type=int, default=80000, help="step at which the beads are frozen")
    p.add_argument("--timesteps", type=int, default=200000, help="step cap (the stopping block ends the run earlier)")
    p.add_argument("--save_interval", type=int, default=500)
    p.add_argument("--first_dir", default="../Results/Planar_beads/")
    p.add_argument("--subfolder", default="1/", help="results subfolder inside --first_dir (config key Subfolder, default 1/)")
    p.add_argument("--out_dir", default="../Config_files/")
    p.add_argument("--subjob", action="store_true", help="also write the sbatch subjob ../Subjobs/<name>_subjob")
    p.add_argument("--time", default="48:00:00", help="--subjob: time limit of the job")
    p.add_argument("--mem", default="16G", help="--subjob: memory of the job")
    a = vars(p.parse_args())
    a["n"] = a.pop("N")
    subjob, time_limit, mem = a.pop("subjob"), a.pop("time"), a.pop("mem")
    path = create_config(**a)
    print(path)
    if subjob:
        from subjob import write_subjob
        name = os.path.splitext(os.path.basename(path))[0]
        print(write_subjob(name + "_subjob", name + "_output.output", "../build/bin/main_cluster {}".format(path),
                           job_name=name, mem=mem, time=time_limit, copy_output_to=a["first_dir"]))
