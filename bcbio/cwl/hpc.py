"""Setup configurations for running on HPC clusters with CWL.

Contains support for setting up configuration inputs for Cromwell.
"""
import os

def create_cromwell_config(args, work_dir):
    """Prepare a cromwell configuration within the current working directory.
    """
    out_file = os.path.join(work_dir, "bcbio-cromwell.conf")
    cl_args, conf_args, scheduler = _args_to_cromwell(args)
    main_config = {"hpc": (HPC_CONFIGS[scheduler] % conf_args) if scheduler else ""}
    with open(out_file, "w") as out_handle:
        out_handle.write(CROMWELL_CONFIG % main_config)
    return out_file

def args_to_cromwell_cl(args):
    """Convert input bcbio arguments into cromwell command line arguments.
    """
    cl_args, conf_args, scheduler = _args_to_cromwell(args)
    return cl_args

def _args_to_cromwell(args):
    """Convert input arguments into cromwell inputs for config and command line.
    """
    scheduler_map = {"sge": "SGE", "slurm": "SLURM", "lsf": "LSF"}
    default_config = {"slurm": {"timelimit": "1-00:00"}, "sge": {}, "lsf": {}}
    # HPC scheduling
    if args.scheduler:
        if args.scheduler not in scheduler_map:
            raise ValueError("Scheduler not yet supported by Cromwell: %s" % args.scheduler)
        if not args.queue:
            raise ValueError("Need to set queue (-q) for running with an HPC scheduler")
        cl = []
        config = default_config[args.scheduler]
        cl.append("-Dbackend.default=%s" % scheduler_map[args.scheduler])
        config["queue"] = args.queue
        for r in args.resources:
            parts = r.split("=")
            if len(parts) == 2:
                key, val = parts
                config[key] = val
        return cl, config, args.scheduler
    # Local multicore runs
    # Avoid overscheduling jobs for local runs by limiting concurrent jobs
    # Longer term would like to keep these within defined core window
    else:
        return ["-Dbackend.providers.Local.config.concurrent-job-limit=1",
                "-Dload-control.memory-threshold-in-mb=1"], {}, None

CROMWELL_CONFIG = """
include required(classpath("application"))

backend {
  providers {
%(hpc)s
  }
}
"""

HPC_CONFIGS = {
"slurm": """
    SLURM {
      actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
      config {
        runtime-attributes = \"\"\"
        Int cores = 1
        Int memory_mb = 2048
        String timelimit = "%(timelimit)s"
        String queue = "%(queue)s"
        \"\"\"
        submit = \"\"\"
            sbatch -J ${job_name} -D ${cwd} -o ${out} -e ${err} -t ${timelimit} -p ${queue} \
            ${"--cpus-per-task=" + cores} \
            --mem=${memory_mb} \
            --wrap "/usr/bin/env bash ${script}"
        \"\"\"
        kill = "scancel ${job_id}"
        check-alive = "squeue -j ${job_id}"
        job-id-regex = "Submitted batch job (\\\\d+).*"
      }
    }

"""
}
