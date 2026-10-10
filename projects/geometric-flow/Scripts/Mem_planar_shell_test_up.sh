#!/bin/tcsh
#$ -cwd
#$ -j y


# Shell test at zero surface tension (sigma~ = 2 KA r^2/KB = 0): how the wrapping transition of one bead on a planar
# membrane changes with the shape of the adhesion shell, W = [1/2 (1 + cos(pi x))]^p, x = (r - a)/(s a).
#   exponent  p = 1, 2, 4      (p = 1 is the shell of the previous runs; p = 2, 4 are cos^4, cos^8 of pi x/2)
#   half width s = 0.25, 0.1, 0.05   (fraction of the bead radius)
# For every (p, s) the interaction strength KI is scanned from "not wrapped" to "fully wrapped": KI = 0.5 ... 1.8 in steps of
# 0.05, w~ = 4 KI r^2/KB = 2 ... 7.2. The previous scan (Wrapping_planar_excess, KA = 0, p = 1, s = 0.25, KI step 0.1) has the
# coverage z = 2 CoverageUnion going from 0.42 at KI = 0.5 over 0.96 at 1.0 and 1.50 at 1.2 to 1.78 at 1.4 and > 1.9 from 1.7 on;
# Deserno's transition at sigma~ = 0 is at w~ = 4 (KI = 1), the narrower the shell the closer the simulated one should come to it.
# Flat start only (Big_planar_mem / Planar_mem with the bead above the plane, no wrapping at the beginning): at zero tension the
# transition is a smooth crossover for the wide shells, if the narrow ones bring back a jump add the full start (the lines of
# the second scan, see Mem_planar_shell_series_up.sh and Studies/deserno_shell/README.md).
#
# Mesh: the shell must be resolved. The simulations follow the continuum soft-shell model (Studies/deserno_shell) to ~0.006 pi kappa
# when the edges near the bead are h <~ 0.5 ell, ell = a sqrt(s_eff) w~^(-1/4), s_eff = s (1, 0.74, 0.54 for p = 1, 2, 4). The
# remesher refine_angle sets the edge near the bead (refine 0.15 -> h ~ 0.10 with size_max 1.0), so it is lowered for the narrower
# shells (more faces: the runs get slower, ask for more time). Optimiser: BFGS_saved_states 60, switch to BFGS-Normal at the given
# step (it also switches when the energy has plateaued).
# Edge bounds and rescale: the old runs used a disk of radius 10 (rescale 1.0 with Big_planar_mem.obj). If Create_subjob_planar.py
# has FLAT_MESH = Planar_mem.obj (a disk of radius 2) the same disk needs rescale 5.0: it is detected below.
# Excess_tension with KA = 0 is no tension at all (as in the KA = 0 row of the old data). The results go to
# ../Results/Shell_test_p<p>_s<s>/
set Nsim=1
set KB = 1
set radius = 1.0
set KA = 0

set Rescale = 1.0
grep -q "^FLAT_MESH = '../../../input/Planar_mem.obj'" Create_subjob_planar.py
if ( $status == 0 ) set Rescale = 5.0
echo "rescale ${Rescale}"

foreach Pw ( 1 2 4 )
foreach Sw ( 0.25 0.1 0.05 )

# refine_angle (edge length near the bead), BFGS -> BFGS-Normal step and slurm time limit for this shell
if ( ${Sw} == "0.25" ) then
set Refine = 0.15
set Switch = 60000
set Time = 12:00:00
endif
if ( ${Sw} == "0.1" ) then
if ( ${Pw} == "1" ) set Refine = 0.15
if ( ${Pw} == "2" ) set Refine = 0.14
if ( ${Pw} == "4" ) set Refine = 0.12
set Switch = 100000
set Time = 24:00:00
endif
if ( ${Sw} == "0.05" ) then
if ( ${Pw} == "1" ) set Refine = 0.11
if ( ${Pw} == "2" ) set Refine = 0.10
if ( ${Pw} == "4" ) set Refine = 0.08
set Switch = 120000
set Time = 48:00:00
endif

set Batch_Tag = "Shell_test_p${Pw}_s${Sw}"

# 27 values: 0.5 ... 1.8 in steps of 0.05
foreach Strg ( 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0 1.05 1.1 1.15 1.2 1.25 1.3 1.35 1.4 1.45 1.5 1.55 1.6 1.65 1.7 1.75 1.8 )

set Unique_Tag = "${Batch_Tag}_Strg_${Strg}_r_${radius}_KA_${KA}_KB_${KB}_Nsim_${Nsim}"

python3 Create_subjob_planar.py --inter_str ${Strg} --radius ${radius} --KA ${KA} --KB ${KB} --tension excess --init flat --rescale ${Rescale} --shell_width ${Sw} --shell_power ${Pw} --size_min 0.01 --size_max 1.0 --refine_angle ${Refine} --saved_states 60 --switch_normal ${Switch} --time ${Time} --batch_tag ${Batch_Tag} --unique_tag ${Unique_Tag}
sbatch ../Subjobs/${Unique_Tag}_subjob

end
end
end
