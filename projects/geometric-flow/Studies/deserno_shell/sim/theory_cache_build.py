"""Pre-solve Deserno's free-membrane problem at the sigma~ values of the planar runs (sigma~ = 2 k/15, k = 1..15)
and store F(z), psi_dot_0(z) on a fine z grid in ../data/theory_cache.npz (used by extract_lines.py)."""
import os, sys, time
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[v] = "1"
import numpy as np
from pathlib import Path
HERE = Path(__file__).resolve().parent
P = HERE.parents[2]
sys.path.insert(0, str(P / "Scripts"))
import deserno_theory as dt

out = HERE.parent / "data" / "theory_cache.npz"
z = np.linspace(0.0, 2.0, 1001)
res = {}
if out.exists():
    old = np.load(out)
    res = {k: old[k] for k in old.files}
res["z"] = z
for k in range(1, 16):
    s = 2.0 * k / 15.0
    key = f"s{k:02d}"
    if key + "_F" in res:
        continue
    t0 = time.time()
    F = np.asarray(dt.free_energy(z, s), float)
    pd = np.asarray(dt.contact_curvature(z, s), float)
    res[key + "_F"] = F
    res[key + "_pd"] = pd
    res[key + "_sigma"] = np.array(s)
    np.savez_compressed(out, **res)
    print(key, s, "%.1f s" % (time.time() - t0), flush=True)
print("done")
