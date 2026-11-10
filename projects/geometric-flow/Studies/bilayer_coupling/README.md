# Bilayer coupling (`Bilayer` energy)

The bilayer-coupling model of Seifert, Berndl and Lipowsky, *Shape transformations of vesicles: Phase diagram for
spontaneous-curvature and bilayer-coupling models*, Phys. Rev. A **44**, 1182 (1991)
(`~/Documents/DDG/Membrane_shape_analitical.pdf`). The two monolayers do not exchange lipids, so for a given area `A`,
volume `V` the Helfrich energy is minimised at a fixed area difference of the monolayers
`ΔA = A_ex − A_in ≈ 2 D M` (eq. 2.6), `M = ∮ H dA = ½ ∮ (C1 + C2) dA` the integrated mean curvature and `D` the distance of
the monolayers. The shape depends on `v = V / (4π/3 R0³)` and on the reduced area difference

    Δa = M / (4π R0),    R0 = sqrt(A / 4π)        (eq. 4.3)

(`Δa = 1` for a sphere, `Δa ≈ 1.03` for the oblate at `v = 0.8`, `Δa ≈ 1.127` for the prolate; Fig. 3 of the paper is the
phase diagram in `(v, Δa)`).

## The energy term

```json
{"Name": "Bilayer", "constants": [K, da0]}
```

    E = K/2 (Δa − Δa0)²

`Δa` is evaluated at the current mesh: `M = ½ Σ_e l_e θ_e` over the interior edges (dihedral angle `θ_e`, positive for a convex
edge with outward normals; this is the sum over the vertices of the `mean_curv` of the Bending terms) and `A` the current
total area. It is scale invariant, as the bending energy. The force is `−dE/dx` with
`dΔa = (4πA)^(-1/2) (dM − M/(2A) dA)` (`F_Bilayer`).

* Constants: `K` (energy) and `da0` (both required; `K = 0` switches the term off). With a large `K` the term is the hard
  constraint of the paper; with a finite one it is the area-difference-elasticity (ADE) energy
  `k_r π/(2 A D²) (ΔA − ΔA0)²`, for which `K = 16 π² k_r` (`k_r` the non-local bending modulus, same units as `K_B`).
* The `Bending` constants are `KB H² A`, `H = (C1 + C2)/2`, so `K_B = 2 κ` of the paper (`G_b = 8πκ` for a sphere). Energies in
  units of `8πκ = 4π K_B` compare directly with the figures.
* Closed membranes with outward normals. Boundary edges do not count and boundary vertices get no force (as in the other
  terms). A reversed mesh has `Δa = −Δa` of the original.
* Penalty like `Area_constraint` and `Volume_constraint`: it stays an energy in every method (it is not replaced by an
  Lagrange multiplier). It has **no Hessian**: with `Newton` the term is missing in the Hessian (the gradient is exact); the
  quasi-Newton methods (`BFGS`, `BFGS-Normal`) are not affected.
* Keep the area and the volume penalties on (`Area_constraint`, `Volume_constraint`): `Δa` and `v` are only meaningful at fixed
  area, and `Bilayer` alone does not fix them.

`Bilayer_example.json` is a vesicle at `v = 0.8` (area target = the area of the mesh, volume target 0.8 of the sphere's) with
`Δa0 = 0.9`; the mesh is `input/Sphere_453.obj`.

## Using it

* **Stiffness.** The residual `Δa − Δa0` of a converged state is (the restoring derivative of the rest of the energy)`/K`.
  Measured on `Sphere_453` at `v = 0.8`, `K_B = 1`: `K = 2000` leaves 0.004 at `Δa0 = 0.7` (0.014 at 0.95, the first step of the
  chain below), so that derivative is about 8. A weak `K` gives the compromise state: `K = 200` settled at `Δa = 1.02` for
  `Δa0 = 0.7` (the bending cost of the cup is larger than the penalty gain).
* **Large jumps do not converge.** Starting the sphere at `Δa0 = 0.7` with `K = 5000` and area/volume penalties of 2e4 stalls the
  line search (the initial penalty, 225, is 20 times the bending energy). Change `Δa0` in steps, restarting from the previous
  `Final_state.obj` (the chain protocol of the other studies). Steps of 0.05 from 0.95 to 0.70 with `K = 2000`,
  `K_A = K_V = 1e4` on `Sphere_453` converged at every step, with `Δa` = 0.964, 0.910, 0.857, 0.805, 0.754, 0.704 and
  `v = 0.800 ± 0.004` and the area within 0.2 % of its target.
* **Output.** With a `Bilayer` term `Output_data.txt` gets a last column `Delta_a` (the reduced area difference of every saved
  row; without the term the file is unchanged).
* `area_difference.py` computes `A`, `v` and `Δa` of OBJ files without the C++ code.

## Checks

`build/bin/Test_bilayer [mesh [workdir]]` (target `Test_bilayer`): `Δa` of the sphere (1.0013 for `Sphere_453`) and independent of
its size, `Δa` of stretched spheres against the quadrature of `∫H dA` of the spheroid, `Δa → −Δa` for the reversed mesh,
the constants read from the input file, and `F_Bilayer` against central finite differences of the energy on a perturbed sphere
and spheroid (relative error 1e-7), and `Calculate_gradient` adds it.

## Snapshots

`plot_slices.py` draws the cross-sections through the axis of symmetry with `Δa` and `v` of each mesh.
`chain_stomatocytes.png` is the chain `Δa0 = 0.95 … 0.60` (`K = 2000`, `v = 0.80`) and `chain_prolates.png` the chain
`Δa0 = 1.05 … 1.125`, both on `Sphere_453`. They reproduce the shapes of Fig. 7b (stomatocytes: 0.7, 0.8, 0.9 look like the paper's
cups) and Fig. 5a (prolate to dumbbell). Below `Δa ≈ 0.7` the neck of the cup closes into the cavity of the limit shape
`L^sto` (the 0.62 of Fig. 7b); the 453-vertex mesh cannot resolve it (the normal force density grows to 5 at 0.655 and 1761 at
0.612), so those two states are only qualitative: use a finer mesh.

## Reproducing the chains

`run_chain.py` runs a chain (every run restarts from the `Final_state.obj` of the previous one; area and reduced volume are fixed
by the first mesh; it stops at the first run that breaks). From this folder, with `build/bin/main_cluster` built:

```bash
# stomatocytes, chain_stomatocytes.png (the runs below 0.70 are the ones the 453-vertex mesh cannot resolve)
python3 run_chain.py --init ../../../../input/Sphere_453.obj --out ~/chains/cup --da0 0.95 0.90 0.85 0.80 0.75 0.70 0.65 0.60

# prolate to dumbbell, chain_prolates.png
python3 run_chain.py --init ../../../../input/Sphere_453.obj --out ~/chains/pro --da0 1.05 1.08 1.10 1.12 1.125

# the 0.62 stomatocyte (da = 0.627, 8000 vertices): remeshing on, from the 0.70 state of the first chain
python3 run_chain.py --init ~/chains/cup/da0_0.700/1/Final_state.obj --out ~/chains/fine --remesh --da0 0.68 0.66 0.64 0.62

# oblate and discocyte (da = 1.033, 1.029, 0.989): the icosphere squashed to v = 0.8, no remeshing
python3 run_chain.py --init ../../../../input/Icosphere.obj --oblate --out ~/chains/oblate --da0 1.03 1.0 0.95

# the cross-sections
python3 plot_slices.py cup.png ~/chains/cup/da0_0.900/1/Final_state.obj ~/chains/cup/da0_0.700/1/Final_state.obj
```

(`plot_slices.py` needs matplotlib: `micromamba run -n mir_mem python plot_slices.py ...`.) Each run keeps its config
(`da0_<value>.json`) and log (`da0_<value>.log`) next to its folder; `Output_data.txt` has the `Delta_a` column. Options:
`--K` (2000), `--KAV` (1e4, K_A and K_V), `--v` (0.8). Checked: the first chain gives `Δa` = 0.964306, 0.909618 for 0.95, 0.90
and the oblate one 1.0333, 1.0293, 0.98896 (the numbers above); the remeshed `--remesh` chain was run with the same settings by an
earlier copy of this script and has not been rerun with this one. The 0.90 oblate step does not converge (the squashed mesh
degenerates and the shape starts to leave the symmetric branch), so the discocyte chain ends at 0.95.
