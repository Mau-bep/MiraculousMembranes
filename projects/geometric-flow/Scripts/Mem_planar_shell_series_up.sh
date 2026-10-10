#!/bin/tcsh
#$ -cwd
#$ -j y


# Shell series of the Deserno study (Studies/deserno_shell): how the simulated phase lines approach those of
# Deserno (arXiv:cond-mat/0303656) when the adhesion shell gets narrower and steeper, with the mesh scaled to the
# shell so that it stays resolved.
#
# Why: the finite range of the adhesion shell acts like a negative contact-line tension of the order
#   tau = - T_p w sqrt(s a / Dc),   Dc = sqrt(w~)/a  (curvature jump at the contact line),
# which moves E, S1, S2 to larger w~ like sqrt(s). With the shell power p, W = [1/2 (1 + cos(pi x))]^p, a steeper shell
# is equivalent to a narrower one, s_eff = s * (0.74, 0.54, 0.39 for p = 2, 4, 8). The simulations follow the
# continuum soft-shell model (Studies/deserno_shell/theory/softshell_direct.py) to ~0.006 pi kappa in energy when the
# edges near the bead are h <~ 0.5 ell, ell = sqrt(s_eff) w~^(-1/4) a. Matching Deserno to 0.1 needs s_eff ~ 0.003, so
# the way to the zero-range limit is a series s = 0.25, 0.15, 0.10, 0.07, 0.05 (p = 8) and an extrapolation of the lines
# quadratically in sqrt(s) (Studies/deserno_shell/theory/line_tension_analysis.md, section on the extrapolation).
#
# Mesh: remesher refine_angle sets the edge length near the bead (target ~ refine_angle/curvature, clamped to
# [size_min, size_max]); refine_angle 0.15 gives edges of ~0.10 near the bead (ell = 0.20 for p = 8, s = 0.25). It is
# scaled with ell, so the number of faces grows like 1/ell^2 and the runs get much slower (~1 h for s = 0.25 on one
# core, several hours for s = 0.05; ask for enough time).
# ALWAYS --rescale 1.0: the template default (3.0) scales the mesh (radius 30 disk, a wrapped mesh around a bead that is not
# rescaled).
#
# Two surface tensions, sigma~ = 2 KA r^2/KB = 0.53 and 1.07 (KA = 0.2667, 0.5333), KB = r = 1, both starts, Excess_tension.
# KI windows (steps of 0.05, w~ steps of 0.2) cover [E_Deserno - 0.4, E_Deserno + 1.6] in w~ (E_Deserno = 5.13, 6.29).
# Compare the flat start (jump S1, partially wrapped energies) and the full start (unwrapping S2, enveloped energies):
# the E line is where Total_E of the two starts crosses. Results in ../Results/Shell_series_<cond>_<init>/.
set Nsim=1
set KB = 1
set radius = 1.0
set Power = 8

# Conditions (name, half width s, refine_angle, step of the BFGS -> BFGS-Normal switch, slurm time). Comment out what
# you do not need; the ones with the finest meshes are the slowest.
foreach Cond ( s25 s15 s10 s07 s05 )
if ( ${Cond} == "s25" ) set Cfg = ( --shell_width 0.25 --refine_angle 0.15 --switch_normal 60000 --time 12:00:00 )
if ( ${Cond} == "s15" ) set Cfg = ( --shell_width 0.15 --refine_angle 0.116 --switch_normal 80000 --time 24:00:00 )
if ( ${Cond} == "s10" ) set Cfg = ( --shell_width 0.10 --refine_angle 0.095 --switch_normal 100000 --time 48:00:00 )
if ( ${Cond} == "s07" ) set Cfg = ( --shell_width 0.07 --refine_angle 0.079 --switch_normal 120000 --time 72:00:00 )
if ( ${Cond} == "s05" ) set Cfg = ( --shell_width 0.05 --refine_angle 0.067 --switch_normal 150000 --time 96:00:00 )

foreach Init ( flat full )

set Batch_Tag = "Shell_series_${Cond}_${Init}"

foreach KA ( 0.2667 0.5333 )
if ( ${KA} == "0.2667" ) set KIs = ( 1.15 1.2 1.25 1.3 1.35 1.4 1.45 1.5 1.55 1.6 1.65 1.7 )
if ( ${KA} == "0.5333" ) set KIs = ( 1.45 1.5 1.55 1.6 1.65 1.7 1.75 1.8 1.85 1.9 1.95 2.0 )

foreach Strg ( ${KIs} )

set Unique_Tag = "${Batch_Tag}_Strg_${Strg}_r_${radius}_KA_${KA}_KB_${KB}_Nsim_${Nsim}"

python3 Create_subjob_planar.py --inter_str ${Strg} --radius ${radius} --KA ${KA} --KB ${KB} --tension excess --init ${Init} --rescale 1.0 --shell_power ${Power} --size_min 0.01 --size_max 1.0 --saved_states 60 ${Cfg} --batch_tag ${Batch_Tag} --unique_tag ${Unique_Tag}
sbatch ../Subjobs/${Unique_Tag}_subjob

end
end
end
end
