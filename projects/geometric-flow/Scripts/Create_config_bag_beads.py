"""Config for a bag of beads: N beads at random inside a membrane that is left to shrink on them (a sack of grapes).

The membrane starts as a sphere (Simplesphere.obj scaled by R0) and has bending and surface tension but NO volume constraint, so the tension shrinks it until it is pressed on the
beads. The beads (radius a, free to move) repel each other: every pair is bonded with the repulsive core "Shifted_LJ" [eps, sigma] (sigma = 2 a, cut at 2^(1/6) sigma = 2.245 a).
The bead-membrane interaction is --mem_inter:
  frenkel    (default) the Frenkel interaction of the wrapping studies (Frenkel_Normal_nopush): zero at the bead surface, repulsive inside, a well of depth --inter_str just outside
             it, cut at the length rc = --rc radii (default 1.2, i.e. 0.2 radii of reach); the beads are inside the membrane ("outside": -1)
  adhesion   the Adhesion interaction with its soft shell (--shell_width, --shell_power, --shell_profile, --shell_alpha), beads inside ("outside": -1)
  repulsive  only the repulsive core "Shifted-LJ" (cut at 2^(1/6) a, strength --inter_str): no adhesion, the membrane is just pressed on the beads
--inter_str is 6 for frenkel and adhesion (the W of the one-bead vesicle runs) and 10 for repulsive unless given. Their start positions are random, inside the sphere, away from the
wall and at least --min_sep (default the cutoff of the bead-bead repulsion, so nothing pushes at t = 0) apart; --seed repeats them.

The start sphere is as big as needed for the beads to fit: R0 = a (N / fill)^(1/3) with --fill (default 0.10) the fraction of the sphere that the beads fill; --R0 sets it directly.

    python Create_config_bag_beads.py -N 10 --KA 1 --KB 2
    python Create_config_bag_beads.py -N 10 --KA 1 --KB 2 --inter_str 3 --rc 1.5
    python Create_config_bag_beads.py -N 10 --KA 1 --KB 2 --mem_inter adhesion --inter_str 2
    python Create_config_bag_beads.py -N 10 --KA 0.5 --KB 2 --seed 3 --subjob

Run from Scripts/ (like the Create_subjob_*.py scripts). The config goes to ../Config_files/ and its path is printed; --subjob also writes the sbatch file ../Subjobs/<name>_subjob.
"""
import argparse
import itertools
import json
import math
import os
import random

from jinja2 import Environment, FileSystemLoader


def random_positions(n, region_radius, min_sep, seed, tries=4000, restarts=200):
    """N points uniform in a ball of the given radius, every pair at least min_sep apart (random sequential addition), or None if they do not fit."""
    rng = random.Random(seed)
    for _ in range(restarts):
        pts = []
        for _ in range(n):
            for _ in range(tries):
                v = [rng.gauss(0.0, 1.0) for _ in range(3)]
                norm = math.sqrt(sum(c * c for c in v))
                rho = region_radius * rng.random() ** (1.0 / 3.0)
                p = [rho * c / norm for c in v]
                if all(math.dist(p, q) >= min_sep for q in pts):
                    pts.append(p)
                    break
            else:
                break  # stuck: start over
        if len(pts) == n:
            return pts
    return None


def create_config(n, KA=1.0, KB=2.0, radius=1.0, eps=10.0, sigma=None, mem_inter="frenkel", inter_str=None, rc=1.2, shell_width=None, shell_power=None,
                  shell_profile=None, shell_alpha=None, bending="Bending_tan", fill=0.10, R0=None, clearance=None,
                  min_sep=None, seed=0, size_min=0.02, size_max=0.5, refine_angle=0.3, saved_states=15, switch_bfgs=2000, switch_normal=80000,
                  timesteps=200000, save_interval=500, first_dir="../Results/Bag_beads/", subfolder="1/", out_dir="../Config_files/"):
    if sigma is None:
        sigma = 2.0 * radius
    r_wca = 2.0 ** (1.0 / 6.0)
    if min_sep is None:
        min_sep = max(r_wca * sigma, 2.0 * radius)  # outside the bead-bead repulsion: no force between the beads at t = 0
    if inter_str is None:
        inter_str = 10.0 if mem_inter == "repulsive" else 6.0
    if mem_inter == "frenkel" and rc <= 1.0:
        raise SystemExit("--rc {:g} must be larger than 1 (the cutoff in bead radii, the potential is zero at the bead surface)".format(rc))
    if clearance is None:
        # bead surface to the membrane at t = 0: outside the reach of the bead-membrane interaction with room to spare
        reach = {"repulsive": r_wca - 1.0, "frenkel": rc - 1.0, "adhesion": shell_width if shell_width is not None else 0.25}[mem_inter]
        clearance = (reach + 0.2) * radius
    if fill <= 0.0 or fill >= 0.5:
        raise SystemExit("--fill {:g} must be between 0 and 0.5".format(fill))

    # The sphere: given, or grown from the fill fraction until the beads fit
    auto = R0 is None
    R = R0 if R0 is not None else radius * (n / fill) ** (1.0 / 3.0)
    pts = None
    for _ in range(40):
        region = R - radius - clearance
        if region > 0.0:
            pts = random_positions(n, region, min_sep, seed)
        if pts is not None or not auto:
            break
        R *= 1.05
    if pts is None:
        raise SystemExit("Could not place {} beads of radius {:g} at least {:g} apart inside a sphere of radius {:g}: use a bigger --R0 (or a smaller --fill) or fewer beads"
                         .format(n, radius, min_sep, R))
    beads = [{"pos": [round(c, 6) for c in p], "partners": [j for j in range(n) if j != k]} for k, p in enumerate(pts)]

    env = Environment(loader=FileSystemLoader(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../Templates/")), trim_blocks=True, lstrip_blocks=True)
    text = env.get_template("Bag_beads.txt").render(
        beads=beads, radius=radius, eps=eps, sigma=sigma, inter_str=inter_str, mem_inter=mem_inter, rc_abs=round(rc * radius, 6),
        **{k: v for k, v in (("shell_width", shell_width), ("shell_power", shell_power), ("shell_profile", shell_profile), ("shell_alpha", shell_alpha)) if v is not None}, KA=KA, KB=KB, bending=bending, R0=round(R, 6), size_min=size_min, size_max=size_max,
        refine_angle=refine_angle, saved_states=saved_states, switch_bfgs=switch_bfgs, switch_normal=switch_normal, timesteps=timesteps,
        save_interval=save_interval, first_dir=first_dir, subfolder=subfolder)
    data = json.loads(text)  # fails here if the template renders invalid JSON

    os.makedirs(out_dir, exist_ok=True)
    name = "Bag_beads_N_{}_{}_str_{}{}_r_{}_KA_{}_KB_{}_seed_{}".format(n, mem_inter, inter_str, "_rc_{}".format(rc) if mem_inter == "frenkel" else "", radius, KA, KB, seed)
    path = os.path.join(out_dir, name + ".json")
    with open(path, "w") as f:
        json.dump(data, f, indent=4)

    dmin = min(math.dist(a, b) for a, b in itertools.combinations(pts, 2)) if n > 1 else float("inf")
    rmax = max(math.sqrt(sum(c * c for c in p)) for p in pts)
    bead_vol = n * 4.0 / 3.0 * math.pi * radius ** 3
    print("{}: {} beads of radius {:g} in a sphere of radius {:g} (they fill {:.0%}), closest beads {:.3f} apart, furthest centre at {:.3f}"
          .format(name, n, radius, R, bead_vol / (4.0 / 3.0 * math.pi * R ** 3), dmin, rmax))
    print("The beads alone, packed, need a sphere of about radius {:.2f}: the bag shrinks to roughly that".format(radius * (n / 0.64) ** (1.0 / 3.0)))
    return path


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-N", type=int, default=10, help="number of beads (default 10)")
    p.add_argument("--KA", type=float, default=1.0, help="surface tension (Surface_tension constant, default 1)")
    p.add_argument("--KB", type=float, default=2.0, help="bending constant (kappa = KB/2, default 2)")
    p.add_argument("--bending", choices=["Bending_tan", "Bending"], default="Bending_tan", help="bending energy")
    p.add_argument("--radius", type=float, default=1.0, help="bead radius")
    p.add_argument("--eps", type=float, default=10.0, help="bead-bead Shifted_LJ epsilon (repulsive core, all pairs)")
    p.add_argument("--sigma", type=float, default=None, help="bead-bead Shifted_LJ sigma (default 2 radii: the repulsion starts at 2^(1/6) sigma)")
    p.add_argument("--mem_inter", choices=["frenkel", "adhesion", "repulsive"], default="frenkel",
                   help="bead-membrane interaction: frenkel (default, Frenkel_Normal_nopush with cutoff --rc), adhesion (soft shell) or repulsive (Shifted-LJ core only)")
    p.add_argument("--inter_str", type=float, default=None,
                   help="bead-membrane strength: the depth of the Frenkel well, the adhesion energy per area, or the repulsive core strength (default 6; 10 for repulsive)")
    p.add_argument("--rc", type=float, default=1.2, help="frenkel: cutoff of the interaction in bead radii (default 1.2, above 1)")
    p.add_argument("--shell_width", type=float, default=None, help="adhesion: half width of the contact shell as a fraction of the radius (code default 0.25)")
    p.add_argument("--shell_power", type=float, default=None, help="adhesion: exponent of the shell (code default 1)")
    p.add_argument("--shell_profile", choices=["cosine", "frenkel", "lj", "morse"], default=None, help="adhesion: radial profile (code default cosine)")
    p.add_argument("--shell_alpha", type=float, default=None, help="adhesion: range of the morse profile")
    p.add_argument("--fill", type=float, default=0.10, help="fraction of the start sphere filled by the beads (sets R0 unless --R0 is given)")
    p.add_argument("--R0", type=float, default=None, help="radius of the start sphere (default from --fill)")
    p.add_argument("--clearance", type=float, default=None, help="smallest distance of a bead surface from the membrane at t = 0 (default a bit more than the 2^(1/6) a repulsion range)")
    p.add_argument("--min_sep", type=float, default=None, help="smallest distance between two bead centres at t = 0 (default the cutoff of the bead-bead repulsion)")
    p.add_argument("--seed", type=int, default=0, help="seed of the random positions")
    p.add_argument("--size_min", type=float, default=0.02, help="remesher smallest edge")
    p.add_argument("--size_max", type=float, default=0.5, help="remesher largest edge")
    p.add_argument("--refine_angle", type=float, default=0.3, help="remesher refine_angle")
    p.add_argument("--saved_states", type=int, default=15, help="BFGS_saved_states")
    p.add_argument("--switch_bfgs", type=int, default=2000, help="step of the Gradient_descent -> BFGS switch")
    p.add_argument("--switch_normal", type=int, default=80000, help="step of the BFGS -> BFGS-Normal switch (the stopping block also switches on a plateau); the beads stay put in BFGS-Normal")
    p.add_argument("--timesteps", type=int, default=200000, help="step cap (the stopping block ends the run earlier)")
    p.add_argument("--save_interval", type=int, default=500)
    p.add_argument("--first_dir", default="../Results/Bag_beads/")
    p.add_argument("--subfolder", default="1/", help="results subfolder inside --first_dir (config key Subfolder, default 1/)")
    p.add_argument("--out_dir", default="../Config_files/")
    p.add_argument("--subjob", action="store_true", help="also write the sbatch subjob ../Subjobs/<name>_subjob")
    p.add_argument("--time", default="48:00:00", help="--subjob: time limit of the job")
    p.add_argument("--mem", default="8G", help="--subjob: memory of the job")
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
