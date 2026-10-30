import os
import json
import math
import argparse
from jinja2 import Environment, FileSystemLoader
from cli_args import number_str
from subjob import write_subjob

# Two beads wrapped by a closed vesicle, separation scan of the pair interaction (van der Wel et al., Sci. Rep. 6, 32825, Fig. 3e / 4a).
# Writes the config of the cold start (from ../Templates/Wrapping_two_vesicle.txt), the chain file (the list of separations that
# Run_two_vesicle_chain.py walks through, warm-starting each one from the previous) and the sbatch subjob that runs the chain.
# See ../Studies/two_beads_vesicle/README.md for the physics and the defaults.

parser = argparse.ArgumentParser(description="Writes the config, the chain file and the sbatch subjob of one separation scan.")
parser.add_argument("--d_cold", type=number_str, default=None,
                    help="distance between the bead centres at the cold start (default 6.0 for beads outside, 9.0 inside, 10.0 on the disk; the beads start tangent to the vesicle and sink in "
                         "(outside) or out (inside) while they wrap: the arc s of the wrapped pair is about 5.5 D_p for 6.0 outside with a floppy vesicle (deep pockets), less "
                         "with the area constraint)")
parser.add_argument("--planar", action="store_true",
                    help="a flat disk with a clamped edge instead of the closed vesicle (input/Big_planar_mem.obj, radius --R_disk): the large-vesicle limit of the experiment "
                         "(D_v/D_p >= 14), no volume or area constraint, only --area tension (default) or none. The beads start tangent above the disk and are free to sink "
                         "in the cold start, then frozen at that height; s is the distance between the bead centres. Reynwar et al. computed an attraction at contact with a barrier at longer range "
                         "for the flat membrane; this setup (a pinned contact line) gives a repulsion, see Studies/two_beads_vesicle/README.md")
parser.add_argument("--free_cold", action="store_true",
                    help="the beads are free (state default) in the cold start and frozen from the first link on: the beads follow the membrane that wraps them instead of the whole vesicle "
                         "having to move relative to them (an inside cold start with frozen beads needs 60000-100000 steps because the vesicle has to slide away from the beads and "
                         "retract the lobe left between them; with free beads 5000 steps got further than 60000). The patch axes stay fixed during the cold start. Default for --planar and --inside, "
                         "except with --wrap coverage: its stiff force pulls a free bead through the membrane (a planar cold start broke within 100 steps) and a free bead inside drifts "
                         "0.9 D_p from where it started, so there the beads are frozen as in the outside chain unless you give --free_cold")
parser.add_argument("--R_disk", type=number_str, default="20.0", help="radius of the planar disk in length units (--planar; Big_planar_mem.obj has radius 10, the mesh is rescaled by R_disk / 10): bead radii for the default --radius 1; "
                                                     "the deformation of a tension sigma~ = 0.02 decays over sqrt(kappa/sigma) = 7 a")
parser.add_argument("--inside", action="store_true",
                    help="the beads are INSIDE the vesicle, adhered to its inner surface and wrapped by outward bulges (the snapshots of the paper's Fig. 4b,c look like that), "
                         "not outside with the membrane invaginating around them (the experiment, the default). Needs the build with the bead key \"outside\": -1 "
                         "(Adhesion; for --wrap coverage the Coverage energy counts the faces that face away from the bead, any build that loads the config has it)")
parser.add_argument("--s_list", type=str, default=None,
                    help="separation s between the beads (--distance: chord or arc, R_v phi), in D_p, of the links after the cold start, in the order they are run "
                         "(default: the far reference at 5.8, then in to 1.45: 5.8,5.0,4.3,3.8,3.4,3.0,2.7,2.4,2.2,2.05,1.9,1.78,1.66,1.55,1.45; "
                         "the cold start itself is not a point of the curve). Below s ~ 1.8 D_p the beads of a floppy 90 %% wrapped vesicle overlap "
                         "(they sit deep, no excluded volume). With --planar s is the distance between the bead centres in D_p, the default list "
                         "5.0,4.2,3.6,3.2,2.9,2.6,2.4,2.2,2.0,1.8,1.65,1.5,1.4,1.3,1.2,1.1 (the beads touch at 1). With --inside the list goes on to 0.9: the beads sit "
                         "1.3 R_v out from the vesicle centre, d = 2 R_c sin(phi/2), they touch at s ~ 0.8 D_p")
parser.add_argument("--distance", choices=["chord", "arc"], default="chord",
                    help="the separation s of --s_list: chord (default) = the straight distance between the bead centres in D_p, the beads touch at s = 1 and never overlap "
                         "(values below 1 are refused); arc = R_v phi, the arc of the sphere (the older chains; beads in deep pockets overlap at s < 1.8 D_p)")
parser.add_argument("--max_step", type=number_str, default="0.5", help="largest change of s (D_p) between two relaxations, bigger moves of the list go through waypoints")
parser.add_argument("--no_return", action="store_true", help="do not walk the list backwards after the last separation (the way back gives the hysteresis, 0.03-0.2 kappa, "
                                                            "and a second value at every s: it doubles the time of the job)")
parser.add_argument("--wrap", choices=["patch", "coverage", "adhesion"], default="adhesion",
                    help="how the wrapped area is set: adhesion (default) = plain Adhesion, a shell around the bead, the contact line is free and the wrapping is "
                         "whatever the energy gives (no patch logic, works with any build, --planar and --inside); patch = Adhesion limited to an angular patch "
                         "(pinned contact line, the paper's theta_Wr; needs the build with patch_angle, and --patch_angle/--patch_width/--free_axis only act with it); "
                         "coverage = the Coverage energy of Wrapping_two_cov")
parser.add_argument("--patch_angle", type=number_str, default=None,
                    help="only with --wrap patch (default 143.13). theta_Wr in degrees, the polar angle of the edge of the patch: 143.13 = 90 percent of the bead area (the paper's text), "
                         "165 = 98 percent (the example of the paper's Fig. S5: an inactive cone of 30 degrees), 101.5 = 60 percent, 90 = 50 percent; "
                         "the active area of a cap of half-angle theta is (1 - cos theta)/2")
parser.add_argument("--patch_width", type=number_str, default=None, help="only with --wrap patch (default 5.0): width in degrees of the smooth ramp inside the edge of the patch")
parser.add_argument("--replica", type=int, default=0, help="replica number: the vesicle starts 0.003 * replica bead radii off the symmetry plane, a negligible change that sends the (chaotic) mesh history of the run elsewhere, for error bars from several replicas")
parser.add_argument("--free_axis", action="store_true", help="only with --wrap patch: let each patch turn: its axis maximises the covered weight at every evaluation, every solve starting from the axis set for the link (default: fixed in space, pointed to the centre of the vesicle at every separation); needs the build with commit 8ced2182 of the branch bead-patch-adhesion (the merge 68af27c6 alone has a warm-started solve that ends on a false maximum, 2 kappa too high); finds the axis within ~25 degrees of the start")
parser.add_argument("--W", type=number_str, default="6.0", help="adhesion energy per area in kappa/a^2 (w~ = 2 W a^2/kappa = 12 for 6: far above 4)")
parser.add_argument("--shell_width", type=number_str, default="0.1", help="half width of the adhesion shell, fraction of the radius")
parser.add_argument("--cov", type=number_str, default="0.6", help="coverage target for --wrap coverage")
parser.add_argument("--Kcov", type=number_str, default="20000.0", help="stiffness of the Coverage energy for --wrap coverage")
parser.add_argument("--area", choices=["tension", "constraint", "excess", "none"], default="tension",
                    help="tension: Surface_tension sigma and a free area (default, floppy: sigma~ = sigma a^2/kappa = 0.02); "
                         "constraint: the paper's Area_constraint, A_t = (1 + excess) 4 pi R_v^2; excess: Excess_tension, free area up to A_t "
                         "then a tension; none: no area term")
parser.add_argument("--sigma", type=number_str, default="0.02", help="tension in kappa/a^2 (tension, and above A_t for excess)")
parser.add_argument("--excess", type=number_str, default="0.05", help="relative excess area (constraint and excess), the paper: 0.05")
parser.add_argument("--KB", type=number_str, default="2.0", help="bending constant KB = 2 kappa, 2 gives kappa = 1: energies in units of kappa")
parser.add_argument("--bending", choices=["Bending", "Bending_tan"], default="Bending",
                    help="the bending energy: Bending (default, edge term = dihedral angle) or Bending_tan (2 tan(angle/2), the energy of all runs before 2026-10-12); "
                         "same KB, same kappa = KB/2, they differ by a nearly constant 0.09 kappa on the wrapped shapes (the dihedral angles are below 0.4 rad)")
parser.add_argument("--KV", type=number_str, default="200000.0", help="stiffness of the volume constraint")
parser.add_argument("--KA", type=number_str, default="200000.0", help="stiffness of the area constraint")
parser.add_argument("--Rv", type=number_str, default="6.25", help="radius of the vesicle in bead radii (the paper: D_v/D_p = 6.25)")
parser.add_argument("--radius", type=number_str, default="1.0", help="bead radius (the length unit)")
parser.add_argument("--size_min", type=number_str, default=None, help="smallest edge length of the remesher (template 0.05)")
parser.add_argument("--size_max", type=number_str, default=None, help="largest edge length of the remesher (template 0.6)")
parser.add_argument("--refine_angle", type=number_str, default=None, help="remesher refine_angle (template 0.3: edges of about 0.2 near a bead)")
parser.add_argument("--switch", type=str, default="100,1500", help="Switch_times of the warm-started links (BFGS, BFGS-Normal)")
parser.add_argument("--prewrap", choices=["coverage", "none"], default="coverage",
                    help="coverage (default, used with --wrap patch only): the beads are first wrapped quickly by the stiff Coverage energy (--prewrap_steps steps, "
                         "target --prewrap_cov), then the patch run starts from that mesh. A patch run started from the tangent contact needs "
                         "30000 steps just to wrap; the pre-wrapped state is already close to its minimum. none: start from the tangent contact "
                         "(the patch axes then point where the vesicle centre is expected, which is only right for shallow pockets). With --inside there is no pre-wrap "
                         "(not made for beads inside): the patch run starts from the tangent contact, 30000-60000 steps (about half an hour) just to wrap")
parser.add_argument("--prewrap_cov", type=number_str, default=None, help="Coverage target of the pre-wrap (default: 0.5 or the area fraction of the patch minus 0.2, whichever is smaller: forcing 0.8 from the tangent contact crumples the membrane)")
parser.add_argument("--prewrap_steps", type=number_str, default="30000", help="steps of the pre-wrap run")
parser.add_argument("--batch_tag", type=str, required=True, help="tag of the batch (also the results folder name)")
parser.add_argument("--unique_tag", type=str, required=True, help="unique tag of this scan (prefix of the config, subjob and output files)")
parser.add_argument("--time", type=str, default="48:00:00", help="time limit of the sbatch job (default 48:00:00: a scan is a chain of about 20 runs of minutes to an hour)")
parser.add_argument("--no_binary_check", action="store_true", help="do not look for the patch option in ../build/bin/main_cluster")
args = parser.parse_args()

if args.wrap != "patch":
    # the patch keys would be dropped without a word: a launcher line written for the old default (--wrap patch) must not silently run another model
    given = [name for name, value in (("--patch_angle", args.patch_angle), ("--patch_width", args.patch_width), ("--free_axis", args.free_axis or None)) if value is not None]
    if given:
        raise SystemExit("{} only act with --wrap patch (the default is now --wrap {}, no patch): add --wrap patch or drop them".format(", ".join(given), args.wrap))
if args.patch_angle is None:
    args.patch_angle = "143.13"
if args.patch_width is None:
    args.patch_width = "5.0"

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
        raise SystemExit("{} does not know \"patch_angle\": build main from the merge of bead-patch-adhesion (68af27c6) or later first (make main_cluster in ../build), "
                         "or use --wrap coverage".format(BINARY))
    if has_patch is None:
        print("Warning: {} not found, the patch option of the build was not checked".format(BINARY))

if args.planar:
    if args.inside:
        raise SystemExit("--planar and --inside do not combine")
    if args.area in ("constraint", "excess"):
        raise SystemExit("--planar has no volume to fix: use --area tension or none")
if args.inside:
    # --wrap coverage needs no check: a build that does not know beads inside refuses the config when it loads ("outside": -1 with Coverage)
    if args.wrap != "coverage" and not args.no_binary_check:
        try:
            with open(BINARY, 'rb') as f:
                content = f.read()
            if b'inside the membrane' not in content:
                raise SystemExit("{} does not know beads inside the membrane (\"outside\": -1 of the Adhesion interaction): build main from the merge of bead-patch-adhesion (68af27c6) or later first".format(BINARY))
        except IOError:
            print("Warning: {} not found, the inside option of the build was not checked".format(BINARY))

# the beads move in the cold start (state "default") for --free_cold and, with Adhesion, for --planar and --inside; with the Coverage energy they are frozen unless --free_cold
free_cold = args.free_cold or ((args.planar or args.inside) and args.wrap != "coverage")

Rv = float(args.Rv)
a = float(args.radius)
d_cold = float(args.d_cold) if args.d_cold is not None else (10.0 if args.planar else (9.0 if args.inside else 6.0))
R_disk = float(args.R_disk)
if args.s_list is None and args.distance == "chord" and not args.planar:
    # touching beads (s = 1) are the last link; the chord is the straight distance, so the arc list's tail below 1.45 is replaced by a finer approach to contact
    args.s_list = "5.8,5.0,4.3,3.8,3.4,3.0,2.7,2.4,2.2,2.05,1.9,1.78,1.66,1.55,1.45,1.35,1.25,1.15,1.07,1.0"
if args.s_list is None:
    args.s_list = ("5.0,4.2,3.6,3.2,2.9,2.6,2.4,2.2,2.0,1.8,1.65,1.5,1.4,1.3,1.2,1.1" if args.planar else
                   "5.8,5.0,4.3,3.8,3.4,3.0,2.7,2.4,2.2,2.05,1.9,1.78,1.66,1.55,1.45,1.3,1.2,1.1,1.0,0.9" if args.inside else
                   "5.8,5.0,4.3,3.8,3.4,3.0,2.7,2.4,2.2,2.05,1.9,1.78,1.66,1.55,1.45")
Batch_dir = '../Results/{}/'.format(args.batch_tag)
V_bar = 4.0 / 3.0 * math.pi * Rv ** 3          # the volume of the ideal sphere, an absolute target for every link
A_sphere = 4.0 * math.pi * Rv ** 2


def create_json():
    os.makedirs("../Config_files/", exist_ok=True)
    env = Environment(loader=FileSystemLoader('../Templates/'))
    template = env.get_template('Wrapping_two_vesicle.txt')

    # The beads start tangent to the vesicle, side by side along x at height 0; the centre of the vesicle is below, on the bisector.
    # Outside (default): each patch axis points from its bead to where the centre of the vesicle will be once the beads are wrapped: the beads
    # sink into the vesicle, their centres end up 1.1 a inside its surface (measured: 5.3 for Rv = 6.25), at a distance Rc from the centre.
    # Inside: the beads start tangent from the inside, 1 a below the surface; the wall wraps them with an outward bulge, the axis points
    # from the centre of the vesicle to the bead (to the membrane); the bead ends up beyond the surface of the main sphere
    if args.planar:
        pass
    elif args.inside:
        yc = -math.sqrt((Rv - a) ** 2 - (d_cold / 2.0) ** 2)
        Rc = Rv + 0.5 * a
        y_wrapped = -math.sqrt(max(Rc ** 2 - (d_cold / 2.0) ** 2, 0.25))
    else:
        yc = -math.sqrt((Rv + a) ** 2 - (d_cold / 2.0) ** 2)
        Rc = Rv - 1.1 * a
        y_wrapped = -math.sqrt(max(Rc ** 2 - (d_cold / 2.0) ** 2, 0.25))
    beads = []
    if args.planar:
        # the disk lies in the plane x = 0 (normal +x), the beads sit above it side by side along y, tangent (centres at height a) and free to sink
        # in the cold start; each patch axis points down to the membrane. The replica shifts the beads by 0.003 a in z
        for y in (d_cold / 2.0, -d_cold / 2.0):
            beads.append({"pos": [a, y, 0.003 * args.replica], "axis": [-1.0, 0.0, 0.0]})
        params = dict(first_dir=Batch_dir, dist=d_cold, beads=beads, displacement=None, Rv=R_disk, rescale=R_disk / 10.0, init_file="../../../input/Big_planar_mem.obj",
                      planar=True, bead_state="default" if free_cold else "froze", radius=a, KB=args.KB, wrap_mode=args.wrap, area_mode=args.area, outside=1,
                      size_max=1.0)  # the beads are free in the cold start of a planar run (they sink to their height)
    else:
        for x in (d_cold / 2.0, -d_cold / 2.0):
            axis = (x, -y_wrapped, 0.0) if args.inside else (-x, y_wrapped, 0.0)
            n = math.hypot(axis[0], axis[1])
            beads.append({"pos": [x, 0.0, 0.0], "axis": [axis[0] / n, axis[1] / n, 0.0]})
        params = dict(first_dir=Batch_dir, dist=d_cold, beads=beads, displacement=[0.0, yc, 0.003 * args.replica], Vbar=V_bar, Rv=Rv, rescale=Rv,
                      radius=a, KB=args.KB, KV=args.KV, wrap_mode=args.wrap, area_mode=args.area, outside=-1 if args.inside else 1,
                      bead_state="default" if free_cold else "froze")
    params["bending"] = args.bending
    if args.wrap == "coverage":
        params.update(cov=args.cov, Kcov=args.Kcov)
    else:
        params.update(W=args.W, shell_width=args.shell_width)
        if args.wrap == "patch":
            params.update(patch_angle=args.patch_angle, patch_width=args.patch_width, patch_free=args.free_axis)
    if args.area == "tension":
        params.update(sigma=args.sigma)
    elif args.area == "constraint":
        params.update(KA=args.KA, Abar=A_sphere * (1.0 + float(args.excess)), nu=-1.0)
    elif args.area == "excess":
        # Excess_tension takes the target as a fraction of the area of the initial (coarse) mesh: use the ideal sphere's area
        # relative to it (the coarse Simplesphere has 99.8 percent of the ideal area)
        params.update(sigma=args.sigma, area_factor=(1.0 + float(args.excess)) * 1.002)
    for key in ("size_min", "size_max", "refine_angle"):
        if getattr(args, key) is not None:
            params[key] = getattr(args, key)

    data = json.loads(template.render(**params))
    Config_path = '../Config_files/{0}_ConfigFile.json'.format(args.unique_tag)
    with open(Config_path, 'w') as file:
        json.dump(data, file, indent=4)

    chain = {"config": Config_path, "binary": BINARY, "batch_dir": Batch_dir, "s_list": [float(x) for x in args.s_list.split(",")],
             "return": not args.no_return, "switch": [int(x) for x in args.switch.split(",")], "nice": 0, "max_step": float(args.max_step)}
    if not args.planar:
        chain["distance"] = args.distance
        if args.distance == "chord" and min(chain["s_list"]) < 1.0:
            raise SystemExit("--distance chord: the beads would overlap at s = %g D_p < 1" % min(chain["s_list"]))
    if args.inside:
        chain["inside"] = True  # the runner points each patch axis away from the volume centre
    if args.planar:
        chain["planar"] = True  # the runner translates the beads along their line, the first link freezes them at the height they sank to
    if args.wrap == "patch" and args.prewrap == "coverage" and not args.inside and not args.planar:
        f_patch = (1.0 - math.cos(math.radians(float(args.patch_angle)))) / 2.0
        pre = dict(params, wrap_mode="coverage", cov=args.prewrap_cov if args.prewrap_cov is not None else round(min(0.5, f_patch - 0.2), 3),
                   Kcov=20000.0, timesteps=int(args.prewrap_steps))
        for key in ("W", "shell_width", "patch_angle", "patch_width"):
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
             'python3 Run_two_vesicle_chain.py {}'.format(Chain_path),
             job_name=args.batch_tag, mem='16G', copy_output_to=Batch_dir, time=args.time)
