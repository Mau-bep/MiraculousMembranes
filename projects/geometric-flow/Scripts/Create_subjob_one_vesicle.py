import os
import json
import math
import argparse
from jinja2 import Environment, FileSystemLoader
from cli_args import number_str
from subjob import write_subjob

# One bead wrapped by a closed vesicle (or a flat disk), a chain of relaxations that scans ONE control parameter, each state warm-started from
# the previous one (the sister of Create_subjob_two_vesicle.py: there the two beads move, here the bead stays and a parameter changes).
# Writes the config of the cold start (from ../Templates/Wrapping_one_vesicle.txt), the chain file (the list of values that
# Run_one_vesicle_chain.py walks through) and the sbatch subjob that runs the chain. See ../Studies/one_bead_vesicle/README.md.

# variable: (what it needs, default list, default largest change between two relaxations (None: the list is walked as given), what it is)
SCAN = {
    "f": ("patch", "0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.98", 0.05, "wrapped area fraction of the patch, f = (1 - cos theta_Wr)/2"),
    "theta": ("patch", "60,75,90,105,120,135,150,165", 6.0, "polar angle theta_Wr of the edge of the patch, degrees"),
    "wt": ("adhesion", "2,2.6,3,3.3,3.5,3.7,3.8,3.9,4,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5,5.1,5.2,5.4,5.6,5.8,6,6.5,7,8,10,12", None,
           "dimensionless interaction strength w~ = 4 W a^2/KB = 2 W a^2/kappa of the Deserno mapping, W = the depth of the potential per area of membrane in contact "
           "(--wrap frenkel, adhesion or patch); the default list is 30 values for the closed vesicle at low tension: a few below the sharp-contact transition of Deserno (w~ = 3.6 for "
           "the Frenkel potential of rc = 1.2, whose contact radius is 1.055 a; 4 for the Adhesion shell), steps of 0.1 from 3.8 to 5.2 (a pilot of the cosine shell at kappa = 21, "
           "sigma~ = 0.02 went from 20 to 94 percent wrapped between 4 and 5; a single 6000-step cold start of the Frenkel potential at 4.5 was 56 percent wrapped and still going), "
           "then a few up to 12"),
    "W": ("adhesion", "1,2,3,4,5,6,8", None, "adhesion energy per unit area of the bead, in the energy unit of the config (kappa/a^2 when KB = 2)"),
    "st": ("tension", "0.5,0.3,0.2,0.1,0.05,0.02,0.01,0.005", None, "dimensionless tension sigma~ = sigma a^2/kappa = 2 sigma a^2/KB (--area tension or excess)"),
    "sigma": ("tension", "3,1,0.3,0.1,0.03,0.01", None, "tension of the membrane, in the energy unit of the config per a^2 (kappa/a^2 when KB = 2; --area tension or excess)"),
    "excess": ("constraint", "0.2,0.15,0.1,0.07,0.05", 0.02, "relative excess area of the fixed-area vesicle (--area constraint or excess)"),
    "volume": ("closed", "1.0,0.95,0.9,0.85,0.8", 0.05, "volume of the vesicle, fraction of the sphere of radius R_v"),
    "cov": ("coverage", "0.2,0.3,0.4,0.5,0.6,0.7,0.8", 0.05, "target of the Coverage energy (--wrap coverage)"),
}

parser = argparse.ArgumentParser(description="Writes the config, the chain file and the sbatch subjob of one scan of a single wrapped bead.")
# --- the chain
parser.add_argument("--scan", choices=list(SCAN), default="wt",
                    help="the control parameter the chain scans (default wt, the interaction strength): " + "; ".join("%s = %s" % (k, v[3]) for k, v in SCAN.items()))
parser.add_argument("--scan_list", type=str, default=None,
                    help="values of the parameter in the order they are run (default per parameter: wt 2 ... 12 (30 values), W 1 ... 8, st 0.5 ... 0.005, sigma 3 ... 0.01, "
                         "f 0.3 ... 0.98, theta 60 ... 165, excess 0.2 ... 0.05, volume 1.0 ... 0.8, cov 0.2 ... 0.8). A decreasing list starts from the end and unwraps (a fully wrapped cold start)")
parser.add_argument("--cold_value", type=number_str, default=None,
                    help="value of the parameter in the cold start (default: the first of the list; if it is another one the chain first walks to the list)")
parser.add_argument("--max_step", type=number_str, default=None,
                    help="largest change of the parameter between two relaxations, bigger changes of the list go through waypoints (default per parameter: f 0.05, theta 6, "
                         "excess 0.02, volume 0.05, cov 0.05, none for W and sigma: the list is walked as given; 0 turns the waypoints off). The contact line of a patch "
                         "advances by a * dtheta per relaxation, which is 0.1 a for 6 degrees")
parser.add_argument("--no_return", action="store_true", help="do not walk the list backwards after the last value (the way back gives the hysteresis and a second value at every point: "
                                                            "it doubles the time of the job)")
# --- geometry
parser.add_argument("--planar", action="store_true",
                    help="a flat disk with a clamped edge instead of the closed vesicle (input/Big_planar_mem.obj, radius --R_disk): no volume or area constraint, only --area tension "
                         "(default) or none. The bead starts tangent above the disk and is free to sink in the cold start, then frozen at that height")
parser.add_argument("--inside", action="store_true",
                    help="the bead is INSIDE the vesicle, adhered to its inner surface and wrapped by an outward bulge, not outside with the membrane invaginating around it (the default). "
                         "Needs the build with the bead key \"outside\": -1 for Adhesion, and --wrap patch or adhesion (not frenkel)")
parser.add_argument("--free_cold", action="store_true",
                    help="the bead is free (state default) in the cold start and frozen from the first link on (default for --planar and --inside; the patch axis stays fixed during the cold start)")
parser.add_argument("--R_disk", type=number_str, default="20.0", help="radius of the planar disk in bead radii (--planar)")
parser.add_argument("--Rv", type=number_str, default="6.25", help="radius of the vesicle in bead radii (the paper: D_v/D_p = 6.25)")
parser.add_argument("--radius", type=number_str, default="1.0", help="bead radius (the length unit)")
# --- how the wrapped area is set
parser.add_argument("--wrap", choices=["frenkel", "adhesion", "patch", "coverage"], default="frenkel",
                    help="frenkel (default) = the Frenkel_Normal_nopush interaction cut at --rc, no patch: the contact line is free; adhesion = the Adhesion interaction with its cosine shell "
                         "(--shell_width, no patch, free contact line); patch = Adhesion limited to an angular patch (pinned contact line, the paper's theta_Wr; needs the build with "
                         "patch_angle; the patch is off unless you ask for it); coverage = the Coverage energy (the wrapped fraction is a target)")
parser.add_argument("--rc", type=number_str, default=None,
                    help="--wrap frenkel: cutoff of the Frenkel potential in bead radii (default 1.2), V = eps alpha (a^2/r^2 - 1)(rc^2/r^2 - 1)^2 of the distance r from the bead centre to the membrane: zero at the bead "
                         "surface r = a, repulsive inside it, minimum -eps at r = a sqrt(3 rc^2/(a^2 + 2 rc^2)) (1.055 a for 1.2), zero value and slope at rc")
parser.add_argument("--patch_angle", type=number_str, default="143.13",
                    help="--wrap patch: theta_Wr in degrees when the chain does not scan it: 143.13 = 90 percent of the bead area, 165 = 98 percent, 101.5 = 60 percent, 90 = 50 percent; "
                         "the active area of a cap of half-angle theta is (1 - cos theta)/2")
parser.add_argument("--patch_width", type=number_str, default="5.0", help="--wrap patch: width in degrees of the smooth ramp inside the edge of the patch")
parser.add_argument("--free_axis", action="store_true",
                    help="--wrap patch: let the patch turn: its axis maximises the covered weight at every evaluation, every solve starting from the axis set for the link (needs the build with "
                         "commit 8ced2182 of the branch bead-patch-adhesion; finds the axis within about 25 degrees of the start)")
parser.add_argument("--W", type=number_str, default="6.0", help="depth of the potential per area of membrane in contact (eps of the Frenkel potential, W of Adhesion), in the energy unit of the config per a^2, "
                                                                 "when the chain does not scan it (w~ = 4 W a^2/KB: 12 for 6 and KB = 2)")
parser.add_argument("--w_tilde", type=number_str, default=None, help="the same as the dimensionless w~ = 4 W a^2/KB (overrides --W)")
parser.add_argument("--shell_width", type=number_str, default=None,
                    help="Adhesion only (--wrap adhesion or patch): half width of the shell, fraction of the radius (default 0.1 with --wrap patch, 0.25 with adhesion: the soft shell of the Deserno study)")
parser.add_argument("--shell_power", type=number_str, default=None, help="Adhesion only: exponent p of the cosine shell [1/2 (1 + cos(pi x))]^p (default 1; steeper for p > 1, the cosine shell only)")
parser.add_argument("--shell_profile", choices=["cosine", "frenkel", "lj", "morse"], default=None,
                    help="Adhesion only: radial profile of the adhesion, the cosine shell (default) or a pair potential with a repulsive core (needs the build with shell_profile). "
                         "This is a profile of the Adhesion interaction, not the Frenkel_Normal_nopush interaction of --wrap frenkel")
parser.add_argument("--shell_alpha", type=number_str, default=None, help="Adhesion only: range of the Morse profile (alpha a, template default 6)")
parser.add_argument("--cov", type=number_str, default="0.6", help="coverage target for --wrap coverage when the chain does not scan it")
parser.add_argument("--Kcov", type=number_str, default=None, help="stiffness of the Coverage energy for --wrap coverage (default 2e4 kappa)")
# --- the vesicle
parser.add_argument("--area", choices=["tension", "constraint", "excess", "none"], default="tension",
                    help="tension: Surface_tension sigma and a free area (default, floppy: sigma~ = sigma a^2/kappa = 0.02); constraint: Area_constraint, A_t = (1 + excess) 4 pi R_v^2 "
                         "(the tension is the multiplier); excess: Excess_tension, free area up to A_t then a tension; none: no area term")
parser.add_argument("--sigma", type=number_str, default="0.02", help="tension, in the energy unit of the config per a^2 (kappa/a^2 when KB = 2), when the chain does not scan it (tension, and above A_t for excess)")
parser.add_argument("--sigma_tilde", type=number_str, default=None, help="the same as the dimensionless sigma~ = sigma a^2/kappa = 2 sigma a^2/KB (overrides --sigma; the floppy vesicles of the experiment have sigma~ < 0.03)")
parser.add_argument("--excess", type=number_str, default="0.05", help="relative excess area when the chain does not scan it (constraint and excess), the paper: 0.05")
parser.add_argument("--volume", type=number_str, default="1.0", help="volume of the vesicle as a fraction of the sphere of radius R_v when the chain does not scan it")
parser.add_argument("--KB", type=number_str, default="2.0", help="bending constant KB = 2 kappa: 2 gives kappa = 1 (energies in units of kappa), 42 gives the kappa = 21 kBT of the paper (energies in kBT). "
                                                                 "The energy scale only enters through the constants below, which are scaled with kappa so that the problem stays the same")
parser.add_argument("--bending", choices=["Bending", "Bending_tan"], default="Bending", help="the bending energy: Bending (edge term = dihedral angle) or Bending_tan (2 tan(angle/2))")
parser.add_argument("--KV", type=number_str, default=None, help="stiffness of the volume constraint (default 2e5 kappa)")
parser.add_argument("--KA", type=number_str, default=None, help="stiffness of the area constraint (default 2e5 kappa)")
# --- mesh
parser.add_argument("--size_min", type=number_str, default=None, help="smallest edge length of the remesher (template 0.05)")
parser.add_argument("--size_max", type=number_str, default=None, help="largest edge length of the remesher (template 0.6)")
parser.add_argument("--refine_angle", type=number_str, default=None,
                    help="remesher refine_angle (template 0.3: edges of about 0.2 near a bead; 0.10 for --wrap frenkel unless given: the well of the potential is 0.2 a wide and a state that "
                         "wraps the bead needs edges of about 0.1 there, with 0.3 a 5000-step relaxation at w~ = 4.5 had 51 vertices in the well and 16 percent of the bead covered, with 0.1 793 and 56 percent)")
# --- how long each run is (the numbers of the stopping block of the config; the links copy the cold start unless a --link_ value is given)
parser.add_argument("--switch", type=str, default="100,1500", help="Switch_times of the warm-started links (BFGS, BFGS-Normal)")
parser.add_argument("--timesteps", type=number_str, default=None, help="largest number of steps of a run (template 200000, never reached in practice); the links copy it")
parser.add_argument("--window", type=number_str, default=None, help="stopping: steps of one window (template 50)")
parser.add_argument("--patience", type=number_str, default=None, help="stopping: windows in a row that must pass the test (template 3)")
parser.add_argument("--bfgs_tol", type=number_str, default=None, help="stopping: relative net energy change of a window below which BFGS hands over to BFGS-Normal (template 5e-05)")
parser.add_argument("--normal_tol", type=number_str, default=None, help="stopping: relative net energy change of a window below which BFGS-Normal ends the run (template 1e-08; 1e-07 "
                                                                         "shortens a link by 20-30 percent for 0.0004 kappa, 1e-06 by 40 percent for 0.001-0.008 kappa)")
parser.add_argument("--stall_steps", type=number_str, default=None, help="stopping: line searches of tiny steps in a row that reset or end a run (template 10)")
parser.add_argument("--link_timesteps", type=number_str, default=None, help="largest number of steps of the links only (the cold start keeps --timesteps)")
parser.add_argument("--link_bfgs_tol", type=number_str, default=None, help="bfgs_tol of the links only")
parser.add_argument("--link_normal_tol", type=number_str, default=None, help="normal_tol of the links only: a tight cold start and looser links")
parser.add_argument("--prewrap", choices=["coverage", "none"], default=None,
                    help="coverage: the bead is first wrapped quickly by the stiff Coverage energy (--prewrap_steps steps, target --prewrap_cov), then the cold start of the chain begins from that mesh, "
                         "the bead where it was (default for --wrap patch, and for any wrap mode when --prewrap_cov is given; the other modes start from the tangent contact, as with none). "
                         "none: start from the tangent contact. Beads outside a closed vesicle only")
parser.add_argument("--prewrap_cov", type=number_str, default=None, help="Coverage target of the pre-wrap (default: for --wrap patch the area fraction of the patch minus 0.2, between 0.1 and 0.5; 0.3 otherwise)")
parser.add_argument("--prewrap_steps", type=number_str, default="8000", help="steps of the pre-wrap run")
# --- the job
parser.add_argument("--replica", type=int, default=0, help="replica number: the vesicle starts 0.003 * replica bead radii off the symmetry axis, a negligible change that sends the (chaotic) mesh history elsewhere, for error bars")
parser.add_argument("--batch_tag", type=str, required=True, help="tag of the batch (also the results folder name)")
parser.add_argument("--unique_tag", type=str, required=True, help="unique tag of this scan (prefix of the config, subjob and output files)")
parser.add_argument("--time", type=str, default="24:00:00", help="time limit of the sbatch job (a chain is about 10-40 runs of a few minutes)")
parser.add_argument("--no_binary_check", action="store_true", help="do not look for the patch option in ../build/bin/main_cluster")
args = parser.parse_args()

BINARY = '../build/bin/main_cluster'
if args.wrap == "patch" and not args.no_binary_check:
    # an old build ignores the unknown "patch_angle" key and would run the bead as a full-sphere adhesion, silently
    try:
        with open(BINARY, 'rb') as f:
            content = f.read()
        has_patch = b'patch_angle' in content and (not args.free_axis or b'patch_free' in content)
    except IOError:
        has_patch = None
    if has_patch is False:
        raise SystemExit("{} does not know \"patch_angle\" (or \"patch_free\"): build main from the merge of bead-patch-adhesion (68af27c6) or later first "
                         "(make main_cluster in ../build), or use --wrap coverage".format(BINARY))
    if has_patch is None:
        print("Warning: {} not found, the patch option of the build was not checked".format(BINARY))
if args.shell_profile is not None and args.shell_profile != "cosine" and not args.no_binary_check:
    try:
        with open(BINARY, 'rb') as f:
            if b'shell_profile' not in f.read():
                raise SystemExit("{} does not know \"shell_profile\": build main from commit 7c4576c2 or later first".format(BINARY))
    except IOError:
        print("Warning: {} not found, the shell_profile option of the build was not checked".format(BINARY))
if args.inside and not args.no_binary_check:
    try:
        with open(BINARY, 'rb') as f:
            if b'inside the membrane' not in f.read():
                raise SystemExit("{} does not know beads inside the membrane (\"outside\": -1 of the Adhesion interaction): build main from the merge of bead-patch-adhesion (68af27c6) or later first".format(BINARY))
    except IOError:
        print("Warning: {} not found, the inside option of the build was not checked".format(BINARY))

# --- what the scan needs
need, default_list, default_step, _ = SCAN[args.scan]
problems = []
if need == "patch" and args.wrap != "patch":
    problems.append("--scan %s needs --wrap patch" % args.scan)
if need == "adhesion" and args.wrap == "coverage":
    problems.append("--scan %s needs --wrap frenkel, adhesion or patch" % args.scan)
if need == "tension" and args.area not in ("tension", "excess"):
    problems.append("--scan sigma needs --area tension or excess (with --area constraint the tension is the multiplier of the constraint)")
if need == "constraint" and args.area not in ("constraint", "excess"):
    problems.append("--scan excess needs --area constraint or excess")
if need == "closed" and args.planar:
    problems.append("--scan volume needs a closed vesicle, not --planar")
if need == "coverage" and args.wrap != "coverage":
    problems.append("--scan cov needs --wrap coverage")
if args.planar and args.inside:
    problems.append("--planar and --inside do not combine")
if args.planar and args.area in ("constraint", "excess"):
    problems.append("--planar has no volume to fix: use --area tension or none")
if args.inside and args.wrap == "coverage":
    problems.append("--inside needs --wrap patch or adhesion: the Coverage energy knows beads outside only")
if args.inside and args.wrap == "frenkel":
    problems.append("--inside needs --wrap patch or adhesion: the Frenkel interaction of this pipeline is set up for a bead outside the vesicle")
if args.wrap == "frenkel":
    adhesion_only = [o for o, v in (("--shell_width", args.shell_width), ("--shell_power", args.shell_power), ("--shell_profile", args.shell_profile), ("--shell_alpha", args.shell_alpha)) if v is not None]
    if adhesion_only:
        problems.append("%s %s to the Adhesion interaction (--wrap adhesion or patch); the Frenkel potential has --rc" % (" and ".join(adhesion_only), "belongs" if len(adhesion_only) == 1 else "belong"))
elif args.rc is not None:
    problems.append("--rc is the cutoff of the Frenkel potential (--wrap frenkel); the shell of Adhesion has --shell_width")
if args.free_axis and args.wrap != "patch":
    problems.append("--free_axis needs --wrap patch")
rc_rel = float(args.rc) if args.rc is not None else 1.2        # the cutoff of the Frenkel potential in bead radii
if args.wrap == "frenkel" and not (1.0 < rc_rel <= 4.0):
    problems.append("--rc must be in (1, 4] bead radii, not %g" % rc_rel)
if problems:
    raise SystemExit("; ".join(problems))

values = [float(x) for x in (args.scan_list if args.scan_list is not None else default_list).split(",")]
cold = float(args.cold_value) if args.cold_value is not None else values[0]
max_step = float(args.max_step) if args.max_step is not None else default_step
if max_step is not None and max_step <= 0:
    max_step = None
for v in values + [cold]:
    if args.scan == "f" and not (0.0 < v <= 1.0):
        raise SystemExit("f must be in (0, 1], not %g" % v)
    if args.scan == "theta" and not (0.0 < v <= 180.0):
        raise SystemExit("theta must be in (0, 180] degrees, not %g" % v)
    if args.scan == "cov" and not (0.0 < v <= 1.0):
        raise SystemExit("cov must be in (0, 1], not %g" % v)
    if args.scan in ("W", "wt", "volume") and not v > 0.0:
        raise SystemExit("%s must be positive, not %g" % (args.scan, v))
    if args.scan in ("sigma", "st", "excess") and v < 0.0:
        raise SystemExit("%s must not be negative, not %g" % (args.scan, v))


def theta_of_f(f):
    return math.degrees(math.acos(1.0 - 2.0 * f))


Rv = float(args.Rv)
a = float(args.radius)
KB = float(args.KB)
kappa = KB / 2.0                      # the code's bending energy is KB * H^2, H = (c1 + c2)/2
R_disk = float(args.R_disk)
Batch_dir = '../Results/{}/'.format(args.batch_tag)
V_ref = 4.0 / 3.0 * math.pi * Rv ** 3          # the volume of the ideal sphere, an absolute target for every link
A_sphere = 4.0 * math.pi * Rv ** 2
# the stiffnesses of the constraints scale with the energy: the problem at kappa = 21 is the one at kappa = 1 with every energy 21 times larger
KV = float(args.KV) if args.KV is not None else 2e5 * kappa
KA = float(args.KA) if args.KA is not None else 2e5 * kappa
Kcov = float(args.Kcov) if args.Kcov is not None else 2e4 * kappa
shell_width = args.shell_width if args.shell_width is not None else ("0.1" if args.wrap == "patch" else "0.25")
w_of_wt = lambda wt: wt * KB / (4.0 * a * a)       # W = w~ KB/(4 a^2): w~ = 4 KI r^2/KB of the Deserno mapping
s_of_st = lambda st: st * KB / (2.0 * a * a)       # sigma = sigma~ KB/(2 a^2): sigma~ = 2 KA r^2/KB

# the parameters of the cold start: the scanned one takes the cold value, the others the options
theta, W, sigma, excess, volume, cov = float(args.patch_angle), float(args.W), float(args.sigma), float(args.excess), float(args.volume), float(args.cov)
if args.w_tilde is not None:
    W = w_of_wt(float(args.w_tilde))
if args.sigma_tilde is not None:
    sigma = s_of_st(float(args.sigma_tilde))
if args.scan == "f":
    theta = theta_of_f(cold)
elif args.scan == "theta":
    theta = cold
elif args.scan == "W":
    W = cold
elif args.scan == "wt":
    W = w_of_wt(cold)
elif args.scan == "sigma":
    sigma = cold
elif args.scan == "st":
    sigma = s_of_st(cold)
elif args.scan == "excess":
    excess = cold
elif args.scan == "volume":
    volume = cold
elif args.scan == "cov":
    cov = cold
if args.wrap == "patch":
    thetas = [theta_of_f(v) if args.scan == "f" else v for v in values + [cold]] if args.scan in ("f", "theta") else [theta]
    if float(args.patch_width) >= min(thetas):
        raise SystemExit("--patch_width %s must be smaller than the smallest patch angle, %.2f degrees" % (args.patch_width, min(thetas)))


def create_json():
    os.makedirs("../Config_files/", exist_ok=True)
    env = Environment(loader=FileSystemLoader('../Templates/'))
    template = env.get_template('Wrapping_one_vesicle.txt')

    # The bead starts tangent to the vesicle at the origin; the vesicle is moved to it. Outside: the centre of the vesicle is Rv + a below the bead and the
    # patch axis points to it (the bead sinks in, the axis stays on the line); inside: the bead touches the wall from the inside, the centre is Rv - a below it and
    # the axis points up, to the membrane the bead wraps. Planar: the disk lies in the plane x = 0, the bead sits above it (centre at height a), the axis points down
    if args.planar:
        beads = [{"pos": [a, 0.0, 0.003 * args.replica], "axis": [-1.0, 0.0, 0.0]}]
        params = dict(first_dir=Batch_dir, beads=beads, displacement=None, Rv=R_disk, rescale=R_disk / 10.0, init_file="../../../input/Big_planar_mem.obj",
                      planar=True, bead_state="default", radius=a, KB=args.KB, wrap_mode=args.wrap, area_mode=args.area, outside=1,
                      size_max=1.0)  # the bead is free in the cold start of a planar run (it sinks to its height)
    else:
        if args.inside:
            beads = [{"pos": [0.0, 0.0, 0.0], "axis": [0.0, 1.0, 0.0]}]
            displacement = [0.0, -(Rv - a), 0.003 * args.replica]
        else:
            beads = [{"pos": [0.0, 0.0, 0.0], "axis": [0.0, -1.0, 0.0]}]
            displacement = [0.0, -(Rv + a), 0.003 * args.replica]
        params = dict(first_dir=Batch_dir, beads=beads, displacement=displacement, Vbar=V_ref * volume, Rv=Rv, rescale=Rv,
                      radius=a, KB=args.KB, KV=KV, wrap_mode=args.wrap, area_mode=args.area, outside=-1 if args.inside else 1,
                      bead_state="default" if (args.free_cold or args.inside) else "froze")
    params["bending"] = args.bending
    if args.wrap == "coverage":
        params.update(cov=cov, Kcov=Kcov)
    elif args.wrap == "frenkel":
        params.update(W=W, rc=rc_rel * a)       # "rc" of a Frenkel bead is an absolute length, its sigma is the bead radius
    else:
        params.update(W=W, shell_width=shell_width)
        for key in ("shell_power", "shell_profile", "shell_alpha"):
            if getattr(args, key) is not None:
                params[key] = getattr(args, key)
        if args.wrap == "patch":
            params.update(patch_angle=theta, patch_width=args.patch_width, patch_free=args.free_axis)
    if args.area == "tension":
        params.update(sigma=sigma)
    elif args.area == "constraint":
        params.update(KA=KA, Abar=A_sphere * (1.0 + excess), nu=-1.0)
    elif args.area == "excess":
        # Excess_tension takes the target as a fraction of the area of the initial mesh of each run: the cold start begins on the coarse Simplesphere,
        # which has 99.8 percent of the ideal area (the runner sets the factor of every link from the area of its own initial mesh)
        params.update(sigma=sigma, area_factor=(1.0 + excess) * 1.002)
    if args.wrap == "frenkel" and args.refine_angle is None:
        params["refine_angle"] = 0.10
    for key in ("size_min", "size_max", "refine_angle"):
        if getattr(args, key) is not None:
            params[key] = getattr(args, key)
    for key, option in (("timesteps", args.timesteps), ("stop_window", args.window), ("stop_patience", args.patience), ("bfgs_tol", args.bfgs_tol),
                        ("normal_tol", args.normal_tol), ("stall_steps", args.stall_steps)):
        if option is not None:
            params[key] = option

    data = json.loads(template.render(**params))
    Config_path = '../Config_files/{0}_ConfigFile.json'.format(args.unique_tag)
    with open(Config_path, 'w') as file:
        json.dump(data, file, indent=4)

    chain = {"config": Config_path, "binary": BINARY, "batch_dir": Batch_dir,
             "scan": {"var": args.scan, "list": values, "cold": cold, "max_step": max_step},
             "return": not args.no_return, "switch": [int(x) for x in args.switch.split(",")], "nice": 0, "area_mode": args.area, "excess": excess}
    if args.inside:
        chain["inside"] = True  # the runner points the patch axis away from the volume centre
    if args.planar:
        chain["planar"] = True  # the patch axis stays pointing down, the first link freezes the bead at the height it sank to
    link_stopping = {}
    if args.link_bfgs_tol is not None:
        link_stopping["bfgs_tol_E"] = float(args.link_bfgs_tol)
    if args.link_normal_tol is not None:
        link_stopping["normal_tol_E"] = float(args.link_normal_tol)
    if link_stopping:
        chain["link_stopping"] = link_stopping
    if args.link_timesteps is not None:
        chain["link_timesteps"] = int(float(args.link_timesteps))
    prewrap = args.prewrap if args.prewrap is not None else ("coverage" if (args.wrap == "patch" or args.prewrap_cov is not None) else "none")
    if prewrap == "coverage" and args.wrap != "coverage" and (args.inside or args.planar):
        print("note: the pre-wrap is for a bead outside a closed vesicle, ignored for --inside and --planar")
    if prewrap == "coverage" and args.wrap != "coverage" and not args.inside and not args.planar:
        if args.prewrap_cov is not None:
            cov_pre = float(args.prewrap_cov)
        elif args.wrap == "patch":
            f_patch = (1.0 - math.cos(math.radians(theta))) / 2.0
            cov_pre = round(max(0.1, min(0.5, f_patch - 0.2)), 3)
        else:
            cov_pre = 0.3
        pre = dict(params, wrap_mode="coverage", cov=cov_pre, Kcov=2e4 * kappa, timesteps=int(float(args.prewrap_steps)))
        for key in ("W", "rc", "shell_width", "patch_angle", "patch_width"):
            pre.pop(key, None)
        Pre_path = '../Config_files/{0}_prewrap.json'.format(args.unique_tag)
        with open(Pre_path, 'w') as file:
            json.dump(json.loads(template.render(**pre)), file, indent=4)
        chain["prewrap_config"] = Pre_path
    Chain_path = '../Config_files/{0}_chain.json'.format(args.unique_tag)
    with open(Chain_path, 'w') as file:
        json.dump(chain, file, indent=4)
    return Config_path, Chain_path


os.makedirs('../Subjobs/', exist_ok=True)
os.makedirs('../Outputs/', exist_ok=True)

Config_path, Chain_path = create_json()

Output_name = '{0}_output.output'.format(args.unique_tag)

write_subjob('{0}_subjob'.format(args.unique_tag),
             Output_name,
             'python3 Run_one_vesicle_chain.py {}'.format(Chain_path),
             job_name=args.batch_tag, mem='16G', copy_output_to=Batch_dir, time=args.time)

print("chain {}: scan {} over {} (cold start at {:g}{}), return leg {}, results in {}".format(
    args.unique_tag, args.scan, ",".join("%g" % v for v in values), cold, ", at most %g between two relaxations" % max_step if max_step else "",
    "off" if args.no_return else "on", Batch_dir))
