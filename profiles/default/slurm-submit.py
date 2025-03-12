import os
import sys
from snakemake.utils import read_job_properties

jobscript = sys.argv[1]
job_properties = read_job_properties(jobscript)

# Extract properties
time = job_properties.get("time", "01:00:00")
partition = job_properties.get("partition", "short")
mem_mb = job_properties.get("mem_mb", 4096)
cpus_per_task = job_properties.get("threads", 1)

# Create SLURM command
cmd = (
    f"sbatch --time={time} --partition={partition} --mem={mem_mb} "
    f"--cpus-per-task={cpus_per_task} {jobscript}"
)

# Submit job
os.system(cmd)
