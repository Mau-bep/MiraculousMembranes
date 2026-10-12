#!/bin/tcsh
#$ -cwd
#$ -j y


# Planar wrapping with other radial profiles of the adhesion weight (Frenkel, Lennard-Jones, Morse), the same
# comparison with Deserno (arXiv:cond-mat/0303656) as the shell exponent series (Mem_planar_shell_series_up.sh), but at
# zero surface tension: KA = 0 (sigma~ = 2 KA r^2/KB = 0), KB = r = 1, so only bending against adhesion.
#
# The profile is the potential turned into a weight, w(r) = -V(r)/depth (bead_geometry::ShellProfile): 1 at r = radius
# (the minimum of the well) for every potential, so KI is still the contact energy per area (w~ = 4 KI r^2/KB).
# Unlike the cosine shell they keep the repulsive core (w < 0 for r below ~ radius) and are cut at rc = radius (1 + s),
# s = --shell_width (up to 4):
#   frenkel  V = alpha (s_F^2/r^2 - 1)(rc^2/r^2 - 1)^2, zero value and slope at rc, s_F chosen so that the minimum is at r = radius
#   lj       V = (r0/r)^12 - 2 (r0/r)^6, shifted to zero at rc (the force jumps there)
#   morse    V = (1 - exp(-alpha (r - r0)/r0))^2 - 1, shifted to zero at rc; alpha = --shell_alpha (default 6, about LJ)
# The depth is renormalised to 1 after the shift, so a small s makes the wall of lj/morse relatively stiff (LJ: depth
# 1 + V(rc) = 0.54 for s = 0.25, 0.064 for s = 0.05). The mesh must resolve the range of the attraction: edges <~ s * radius / 2.
# Tests: build/bin/Test_shell_profile.
#
# Mesh as the s = 0.25 row of the shell series (refine_angle 0.15, edges ~0.1 near the bead). ALWAYS --rescale 1.0.
# Compare with the cosine shell of the same width (Cond cos below) and the p-series; flat and full starts, the E line is
# where Total_E of the two starts crosses. Results in ../Results/Potentials_<cond>_<init>/.
set Nsim=1
set KB = 1
set radius = 1.0
set KA = 0
# 17 values, w~ = 4 KI from 3.2 to 9.6 in steps of 0.4. At sigma~ = 0 Deserno's wrapping transition E is at w~ = 4 (KI = 1.0,
# deserno_table.csv); the finite range of the adhesion moves it to larger KI (the cosine shell scans put it near w~ = 5 and above), so the
# window reaches w~ = 9.6; widen it if the lines are not bracketed
set KIs = ( 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3 2.4 )

# Conditions: name, profile and range. Comment out what you do not need.
# LJ: rc = 2.5 and 2 sigma_LJ, sigma_LJ = radius / 2^(1/6) (the minimum is at the radius), i.e. rc = 2.2272 and 1.7818 radius,
#     shell_width = rc/radius - 1 = 1.2272468 and 0.7817974.
# Morse: a = alpha/radius = 2 and 10 (radius = 1), cut at 2.5 sigma_LJ like lj25. a = 2 is a soft wide well (w = 0 at 0.68 radius,
#     attraction down to half at 1.52 radius), a = 10 a short stiff one (w = 0 at 0.93 radius, w(0.8 radius) = -40, half at 1.12 radius).
#     (a = 0.25 and 0.1 were rejected: almost no repulsive core, the membrane could sink into the bead.)
foreach Cond ( cos frenkel lj25 lj20 morse2 morse10 )
if ( ${Cond} == "cos" )      set Cfg = ( --shell_profile cosine --shell_width 0.25 )
if ( ${Cond} == "frenkel" )  set Cfg = ( --shell_profile frenkel --shell_width 0.25 )
if ( ${Cond} == "lj25" )     set Cfg = ( --shell_profile lj --shell_width 1.2272468 )
if ( ${Cond} == "lj20" )     set Cfg = ( --shell_profile lj --shell_width 0.7817974 )
if ( ${Cond} == "morse2" )   set Cfg = ( --shell_profile morse --shell_width 1.2272468 --shell_alpha 2 )
if ( ${Cond} == "morse10" )  set Cfg = ( --shell_profile morse --shell_width 1.2272468 --shell_alpha 10 )

foreach Init ( flat full )

set Batch_Tag = "Potentials_${Cond}_${Init}"

foreach Strg ( ${KIs} )

set Unique_Tag = "${Batch_Tag}_Strg_${Strg}_r_${radius}_KA_${KA}_KB_${KB}_Nsim_${Nsim}"

python3 Create_subjob_planar.py --inter_str ${Strg} --radius ${radius} --KA ${KA} --KB ${KB} --tension excess --init ${Init} --rescale 1.0 --size_min 0.01 --size_max 1.0 --refine_angle 0.15 --switch_normal 60000 --saved_states 60 ${Cfg} --batch_tag ${Batch_Tag} --unique_tag ${Unique_Tag} --time 12:00:00
sbatch ../Subjobs/${Unique_Tag}_subjob

end
end
end
