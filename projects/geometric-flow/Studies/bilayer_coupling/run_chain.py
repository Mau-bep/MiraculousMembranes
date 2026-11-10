"""A chain of Bilayer runs: Delta_a0 is stepped through the given values, every run restarting from the Final_state.obj of the
previous one, at a fixed area and reduced volume (the area is the one of the first mesh).

    python run_chain.py --init ../../../../input/Sphere_453.obj --out ~/chains/cup --da0 0.95 0.90 0.85 0.80 0.75 0.70

    # finer mesh near the closing neck of the stomatocyte: remeshing on, continue from a state of the chain above
    python run_chain.py --init ~/chains/cup/da0_0.700/1/Final_state.obj --out ~/chains/fine --da0 0.68 0.66 0.64 0.62 --remesh

    # discocytes: start from the oblate (the mesh stretched along z to the reduced volume), no remeshing so that the mirror
    # symmetry z -> -z of the icosphere is kept
    python run_chain.py --init ../../../../input/Icosphere.obj --oblate --out ~/chains/oblate --da0 1.03 1.0 0.95 0.9

Each run is <out>/da0_<value>/1/ (Output_data.txt with the Delta_a column, membrane_*.obj, Final_state.obj); the config of
the run is <out>/da0_<value>.json. The template is Bilayer_example.json (Bending 1, K_A = K_V = 1e4, BFGS-Normal, the mesh is remeshed at the start only).
The area target is the area of --init, the volume target is v (4 pi / 3) R0^3 with R0 = sqrt(A / 4 pi).
"""
import argparse
import json
import math
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def load_obj(path):
    V, F = [], []
    for line in open(path):
        t = line.split()
        if t and t[0] == "v":
            V.append([float(x) for x in t[1:4]])
        elif t and t[0] == "f":
            F.append([int(x.split("/")[0]) - 1 for x in t[1:4]])
    return np.array(V), np.array(F)


def area_volume(V, F):
    P = V[F]
    A = 0.5 * np.linalg.norm(np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0]), axis=1).sum()
    return A, np.einsum("ij,ij->", P[:, 0], np.cross(P[:, 1], P[:, 2])) / 6


def oblate_factor(v):
    """The factor c of z for which the unit sphere becomes the oblate spheroid (1, 1, c) of reduced volume v."""
    def reduced(c):
        e = math.sqrt(1 - c * c)
        area = 2 * math.pi * (1 + c * c / e * math.atanh(e))
        return c / (area / (4 * math.pi)) ** 1.5
    lo, hi = 1e-3, 1 - 1e-9
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = (lo, mid) if reduced(mid) > v else (mid, hi)  # the reduced volume grows with c
    return 0.5 * (lo + hi)


def write_oblate(src, dst, v):
    c = oblate_factor(v)
    with open(dst, "w") as out:
        for line in open(src):
            t = line.split()
            if t and t[0] == "v":
                line = "v %.17g %.17g %.17g\n" % (float(t[1]), float(t[2]), float(t[3]) * c)
            out.write(line)
    print("oblate: z scaled by %.6f" % c)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--init", required=True, help="first mesh (OBJ)")
    ap.add_argument("--out", required=True, help="folder of the chain")
    ap.add_argument("--da0", type=float, nargs="+", required=True, help="the values of Delta_a0, in order")
    ap.add_argument("--K", type=float, default=2000.0, help="K of the Bilayer term (default 2000)")
    ap.add_argument("--KAV", type=float, default=1e4, help="K_A and K_V of the area and volume penalties (default 1e4)")
    ap.add_argument("--v", type=float, default=0.8, help="reduced volume (default 0.8)")
    ap.add_argument("--oblate", action="store_true", help="stretch --init along z to the oblate of reduced volume v first")
    ap.add_argument("--remesh", action="store_true", help="remeshing on throughout (finer mesh, breaks mirror symmetry)")
    ap.add_argument("--exe", default=os.path.join(HERE, "../../build/bin/main_cluster"), help="main_cluster")
    args = ap.parse_args()

    out = os.path.abspath(args.out)
    os.makedirs(out, exist_ok=True)
    init = os.path.abspath(args.init)
    if args.oblate:
        init_oblate = os.path.join(out, "oblate_start.obj")
        write_oblate(init, init_oblate, args.v)
        init = init_oblate
    A, _ = area_volume(*load_obj(init))
    V_target = args.v * 4 * math.pi / 3 * (A / (4 * math.pi)) ** 1.5
    print("area target %.6f, volume target %.6f (v = %.3f)" % (A, V_target, args.v))

    template = json.load(open(os.path.join(HERE, "Bilayer_example.json")))
    for da0 in args.da0:
        name = "da0_%.3f" % da0
        cfg = json.loads(json.dumps(template))
        cfg["init_file"] = init
        cfg["first_dir"] = os.path.join(out, name) + "/"
        cfg["Subfolder"] = "1/"
        cfg["timesteps"] = 20000
        if args.remesh:
            cfg["Switches"] = ["BFGS-Normal", "Save_all"]
            cfg["Switch_times"] = [1, -1]
            cfg["remesh_every"] = 10
            cfg["remeshing"] = True
            cfg["remesher"] = {"aspect_min": 0.2, "refine_angle": 0.4, "refine_compression": 0.0001, "refine_velocity": 1.0,
                               "size_max": 0.06, "size_min": 0.01}
        elif args.oblate:
            cfg["remeshing"] = False  # never remesh: the mirror symmetry of the mesh is what keeps the discocytes
        # else: "remeshing" is not given, as in Bilayer_example.json (the mesh is remeshed at the start, then No_remesh)
        for e in cfg["Energies"]:
            if e["Name"] == "Bilayer":
                e["constants"] = [args.K, da0]
            elif e["Name"] == "Area_constraint":
                e["constants"] = [args.KAV, A]
            elif e["Name"] == "Volume_constraint":
                e["constants"] = [args.KAV, V_target]
        path = os.path.join(out, name + ".json")
        json.dump(cfg, open(path, "w"), indent=1)

        result = subprocess.run([args.exe, path], capture_output=True, text=True)
        log = open(os.path.join(out, name + ".log"), "w")
        log.write(result.stdout + result.stderr)
        log.close()
        final = os.path.join(out, name, "1", "Final_state.obj")
        status = [l for l in result.stdout.splitlines() if "converged" in l or "broke" in l]
        da = ""
        output = os.path.join(out, name, "1", "Output_data.txt")
        if os.path.exists(output):
            da = open(output).read().splitlines()[-1].split()[-1]
        print("da0 = %.3f -> Delta_a = %s | %s" % (da0, da, (status[0][:90] if status else "no convergence message")), flush=True)
        if not os.path.exists(final) or any("broke" in l for l in status):
            sys.exit("the run did not finish properly, the chain stops here (see %s.log)" % name)
        init = final


if __name__ == "__main__":
    main()
