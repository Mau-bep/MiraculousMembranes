"""Shared sbatch subjob writer for the Create_subjob_*.py scripts.

The subjob text lives in ../Templates/Subjob.txt (Jinja); write_subjob() only fills it in.
Paths are relative to the Scripts folder, like in the Create_subjob_*.py scripts.
"""
import os
from pathlib import Path

from jinja2 import Environment, FileSystemLoader, StrictUndefined

TEMPLATE_DIR = Path(__file__).resolve().parent.parent / 'Templates'

_env = Environment(loader=FileSystemLoader(str(TEMPLATE_DIR)), undefined=StrictUndefined,
                   trim_blocks=True, lstrip_blocks=True, keep_trailing_newline=True)


def write_subjob(subjob_name, output_name, command, job_name="Mem3DGpa", time="10:00:00",
                 mem="3G", cpus=1, conda=True, copy_output_to=None, mail_user=None,
                 error=None, subjob_dir='../Subjobs/', output_dir='../Outputs/'):
    """Write the sbatch file subjob_dir/subjob_name and return its path.

    command         what srun runs, e.g. '../build/bin/main_cluster ../Config_files/x.json'
    output_name     file name (in output_dir) that slurm writes stdout/stderr to
    cpus            used for both --cpus-per-task and OMP_NUM_THREADS
    conda           load conda and activate mir_membranes before running
    copy_output_to  if given, the output file is copied there (as output_name) when the run ends
    mail_user       if given, slurm mails this address on start / end / failure
    error           if given, a separate --error file (slurm patterns like %x_%j.err work)
    """
    os.makedirs(subjob_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)

    text = _env.get_template('Subjob.txt').render(
        job_name=job_name, output_name=output_name, output_path=output_dir + output_name,
        command=command, time=time, mem=mem, cpus=cpus, conda=conda,
        copy_output_to=copy_output_to, mail_user=mail_user, error=error)

    subjob_path = subjob_dir + subjob_name
    with open(subjob_path, 'w') as f:
        f.write(text)
    return subjob_path
