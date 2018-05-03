"""Setup configurations for running on HPC clusters with CWL.

Contains support for setting up configuration inputs for Cromwell.
"""
import os

def create_cromwell_config(args, work_dir):
    """Prepare a cromwell configuration within the current working directory.
    """
    out_file = os.path.join(work_dir, "bcbio-cromwell.conf")
    cl_args, conf_args, scheduler = _args_to_cromwell(args)
    main_config = {"hpc": (HPC_CONFIGS[scheduler] % conf_args) if scheduler else "",
                   "work_dir": work_dir}
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
    default_config = {"slurm": {"timelimit": "1-00:00", "account": ""},
                      "sge": {"memtype": "mem_type", "pename": "smp"},
                      "lsf": {},
                      "torque": {"walltime": "1-00:00", "account": ""},
                      "pbspro": {"walltime": "1-00:00", "account": ""}}
    prefixes = {("account", "slurm"): "-A ", ("account", "pbspro"): "-A "}
    cl = ["-Dload-control.memory-threshold-in-mb=1"]
    # HPC scheduling
    if args.scheduler:
        if args.scheduler not in default_config:
            raise ValueError("Scheduler not yet supported by Cromwell: %s" % args.scheduler)
        if not args.queue:
            raise ValueError("Need to set queue (-q) for running with an HPC scheduler")
        config = default_config[args.scheduler]
        cl.append("-Dbackend.default=%s" % args.scheduler.upper())
        config["queue"] = args.queue
        for r in args.resources:
            parts = r.split("=")
            if len(parts) == 2:
                key, val = parts
                config[key] = prefixes.get((key, args.scheduler), "") + val
        return cl, config, args.scheduler
    # Local multicore runs
    # Avoid overscheduling jobs for local runs by limiting concurrent jobs
    # Longer term would like to keep these within defined core window
    else:
        return ["-Dbackend.providers.Local.config.concurrent-job-limit=1"] + cl, {}, None

CROMWELL_CONFIG = """
include required(classpath("application"))

system {
  workflow-restart = true
}
call-caching {
  enabled = true
}

database {
  profile = "slick.jdbc.HsqldbProfile$"
  db {
    driver = "org.hsqldb.jdbcDriver"
    url = "jdbc:hsqldb:file:%(work_dir)s/persist/metadata;shutdown=false;hsqldb.tx=mvcc"
    connectionTimeout = 20000
  }
}

backend {
  providers {
    Local {
      config {
        runtime-attributes = \"\"\"
        String? docker
        String? docker_user
        Int? cpu
        Int? memory_mb
        Int? cpuMin
        Int? cpuMax
        Int? memoryMin
        Int? memoryMax
        String? outDirMin
        String? outDirMax
        String? tmpDirMin
        String? tmpDirMax
      \"\"\"
      }
    }
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
        Int cpu = 1
        Int memory_mb = 2048
        String queue = "%(queue)s"
        String timelimit = "%(timelimit)s"
        String account = "%(account)s"
        \"\"\"
        submit = \"\"\"
            sbatch -J ${job_name} -D ${cwd} -o ${out} -e ${err} -t ${timelimit} -p ${queue} \
            ${"--cpus-per-task=" + cpu} --mem=${memory_mb} ${account} \
            --wrap "/usr/bin/env bash ${script}"
        \"\"\"
        kill = "scancel ${job_id}"
        check-alive = "squeue -j ${job_id}"
        job-id-regex = "Submitted batch job (\\\\d+).*"
      }
    }
""",
"sge": """
    SGE {
      actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
      config {
        runtime-attributes = \"\"\"
        Int cpu = 1
        Int memory_mb = 2048
        String queue = "%(queue)s"
        String pename = "%(pename}s"
        String memtype = "%(memtype)s"
        \"\"\"
        submit = \"\"\"
        qsub -V -w w -j y -N ${job_name} -wd ${cwd} \
        -o ${out} -e ${err} -q ${queue} \
        -pe ${pename} ${cpu} ${"-l " + mem_type + "=" + memory_mb + "m"} \
        /usr/bin/env bash ${script}
        \"\"\"
        kill = "qdel ${job_id}"
        check-alive = "qstat -j ${job_id}"
        job-id-regex = "(\\\\d+)"
      }
    }
""",
"pbspro": """
    PBSPRO {
      actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
      config {
        runtime-attributes = \"\"\"
        Int cpu = 1
        Int memory_mb = 2048
        String queue = "%(queue)s"
        String account = "%(account)s"
        String walltime = "%(walltime)s"
        \"\"\"
        submit = \"\"\"
        qsub -V -N ${job_name} -W sandbox=${cwd} \
        -o ${out} -e ${err} -q ${queue} \
        -l select=1:ncpus=${cpu}:mem=${memory_mb}mb -l walltime=${walltime} \
        /usr/bin/env bash ${script}
        \"\"\"
        kill = "qdel ${job_id}"
        check-alive = "qstat -j ${job_id}"
        job-id-regex = "(\\\\d+)"
      }
    }

""",
"torque": """
    TORQUE {
      actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
      config {
        runtime-attributes = \"\"\"
        Int cpu = 1
        Int memory_mb = 2048
        String queue = "%(queue)s"
        String account = "%(account)s"
        String walltime = "%(walltime)s"
        \"\"\"
        submit = \"\"\"
        qsub -V -N ${job_name} -W sandbox=${cwd} \
        -o ${out} -e ${err} -q ${queue} \
        -l nodes=1:ppn=${cpu} -l mem=${memory_mb}mb -l walltime=${walltime} \
        /usr/bin/env bash ${script}
        \"\"\"
        kill = "qdel ${job_id}"
        check-alive = "qstat -j ${job_id}"
        job-id-regex = "(\\\\d+)"
      }
    }
"""
}
