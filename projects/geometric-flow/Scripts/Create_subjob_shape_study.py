import argparse
from cli_args import number_str
from subjob import write_subjob

# The study of the remesher settings on vesicle shapes (Study_shape_mesh.py) as ONE sbatch job: the runs of every case and setting (about 50 relaxations, `--workers` at a time on the cpus of the job),
# then the remesh stress test (the three stress meshes of ../Studies/shape_chains/stress_meshes are made if missing, 60 runs of 150 steps). Results: ../Results/Shape_mesh_study/ (a chain folder per run,
# mesh_study.csv, stress_robust.csv). The figures and tables (`Study_shape_mesh.py plot`, `stress_plot`) are made where matplotlib is: python3 Study_shape_mesh.py collect ; plot --out ../Studies/shape_chains/mesh
# See ../Studies/shape_chains/README.md for what each case and group is.

parser = argparse.ArgumentParser(description="Writes the sbatch subjob of the study of the remesher settings on vesicle shapes.")
parser.add_argument("--unique_tag", type=str, default="Shape_mesh_study", help="prefix of the subjob and output files")
parser.add_argument("--workers", type=int, default=8, help="runs at the same time (and cpus of the job): every run is one core")
parser.add_argument("--cap", type=number_str, default="30000", help="largest number of steps of a run (a run that reaches it is marked capped)")
parser.add_argument("--cases", nargs="*", default=None, help="only these cases (sphere, prolate1, prolate2, pear, pear3)")
parser.add_argument("--groups", nargs="*", default=None, help="only these groups (angle, size_max, size_min, aspect, integrator)")
parser.add_argument("--binary", type=str, default="../build/bin/main_cluster", help="the main_cluster that runs the relaxations")
parser.add_argument("--no_stress", action="store_true", help="leave the remesh stress test out")
parser.add_argument("--time", type=str, default="48:00:00", help="time limit of the sbatch job")
parser.add_argument("--mem", type=str, default="32G", help="memory of the sbatch job")
args = parser.parse_args()

run = "python3 Study_shape_mesh.py run --workers {} --cap {} --binary {}".format(args.workers, args.cap, args.binary)
if args.cases:
    run += " --cases " + " ".join(args.cases)
if args.groups:
    run += " --groups " + " ".join(args.groups)
steps = [run]
if not args.no_stress:
    steps += ["python3 Study_shape_mesh.py stress_make --binary {}".format(args.binary),
              "python3 Study_shape_mesh.py stress_robust --workers {} --binary {} --out ../Studies/shape_chains/mesh".format(args.workers, args.binary)]
steps += ["python3 Study_shape_mesh.py collect"]
command = 'bash -c "{}"'.format(" && ".join(steps))

write_subjob('{}_subjob'.format(args.unique_tag), '{}_output.output'.format(args.unique_tag), command, job_name='Shape_mesh_study', mem=args.mem, cpus=args.workers, time=args.time,
             copy_output_to='../Results/Shape_mesh_study')
print("study job {}: {} workers, cap {}, {}, time limit {}, results in ../Results/Shape_mesh_study/".format(args.unique_tag, args.workers, args.cap, "stress test included" if not args.no_stress else "no stress test", args.time))
