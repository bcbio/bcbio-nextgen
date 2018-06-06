"""Setup configurations for running on HPC clusters with CWL.

Contains support for setting up configuration inputs for Cromwell.
"""
import os

def create_cromwell_config(args, work_dir):
    """Prepare a cromwell configuration within the current working directory.
    """
    docker_attrs = ["String? docker", "String? docker_user"]
    cwl_attrs = ["Int? cpuMin", "Int? cpuMax", "Int? memoryMin", "Int? memoryMax", "String? outDirMin",
                 "String? outDirMax", "String? tmpDirMin", "String? tmpDirMax"]
    out_file = os.path.join(work_dir, "bcbio-cromwell.conf")
    run_config = _load_custom_config(args.runconfig) if args.runconfig else {}
    std_args = {"docker_attrs": "" if args.no_container else "\n        ".join(docker_attrs),
                "submit_docker": 'submit-docker: ""' if args.no_container else "",
                "cwl_attrs": "\n        ".join(cwl_attrs),
                "filesystem": FILESYSTEM_CONFIG,
                "database": run_config.get("database", DATABASE_CONFIG % {"work_dir": work_dir})}
    cl_args, conf_args, scheduler = _args_to_cromwell(args)
    conf_args.update(std_args)
    main_config = {"hpc": (HPC_CONFIGS[scheduler] % conf_args) if scheduler else "",
                   "work_dir": work_dir}
    main_config.update(std_args)
    # Local run always seems to need docker set because of submit-docker in default configuration
    # Can we unset submit-docker based on configuration so it doesn't inherit?
    #main_config["docker_attrs"] = "\n        ".join(docker_attrs)
    with open(out_file, "w") as out_handle:
        out_handle.write(CROMWELL_CONFIG % main_config)
    return out_file

def _load_custom_config(run_config):
    """Load custom configuration input HOCON file for cromwell.
    """
    from pyhocon import ConfigFactory, HOCONConverter, ConfigTree
    conf = ConfigFactory.parse_file(run_config)
    out = {}
    if "database" in conf:
        out["database"] = HOCONConverter.to_hocon(ConfigTree({"database": conf.get_config("database")}))
    return out

def args_to_cromwell_cl(args):
    """Convert input bcbio arguments into cromwell command line arguments.
    """
    cl_args, conf_args, scheduler = _args_to_cromwell(args)
    return cl_args

def _args_to_cromwell(args):
    """Convert input arguments into cromwell inputs for config and command line.
    """
    default_config = {"slurm": {"timelimit": "1-00:00", "account": ""},
                      "sge": {"memtype": "mem_type", "pename": "smp"},
                      "lsf": {},
                      "htcondor": {},
                      "torque": {"walltime": "24:00:00", "account": ""},
                      "pbspro": {"walltime": "24:00:00", "account": "",
                                 "cpu_and_mem": "-l select=1:ncpus=${cpu}:mem=${memory_mb}mb"}}
    prefixes = {("account", "slurm"): "-A ", ("account", "pbspro"): "-A "}
    custom = {("noselect", "pbspro"): ("cpu_and_mem", "-l ncpus=${cpu} -l mem=${memory_mb}mb")}
    cl = []
    config = {}
    # HPC scheduling
    if args.scheduler:
        if args.scheduler not in default_config:
            raise ValueError("Scheduler not yet supported by Cromwell: %s" % args.scheduler)
        if not args.queue and args.scheduler not in ["htcondor"]:
            raise ValueError("Need to set queue (-q) for running with an HPC scheduler")
        config = default_config[args.scheduler]
        cl.append("-Dbackend.default=%s" % args.scheduler.upper())
        config["queue"] = args.queue
        for rs in args.resources:
            for r in rs.split(";"):
                parts = r.split("=")
                if len(parts) == 2:
                    key, val = parts
                    config[key] = prefixes.get((key, args.scheduler), "") + val
                elif len(parts) == 1 and (parts[0], args.scheduler) in custom:
                    key, val = custom[(parts[0], args.scheduler)]
                    config[key] = val
        return cl, config, args.scheduler
    return cl, config, args.scheduler


FILESYSTEM_CONFIG = """
      filesystems {
        local {
          localization: ["soft-link"]
          caching {
            duplication-strategy: ["soft-link"]
            hashing-strategy: "path"
          }
        }
      }
"""

DATABASE_CONFIG = """
database {
  profile = "slick.jdbc.HsqldbProfile$"
  db {
    driver = "org.hsqldb.jdbcDriver"
    url = "jdbc:hsqldb:file:%(work_dir)s/persist/metadata;shutdown=false;hsqldb.tx=mvcc"
    connectionTimeout = 100000
  }
}
"""

CROMWELL_CONFIG = """
include required(classpath("application"))

system {
  workflow-restart = true
}
call-caching {
  enabled = true
}
load-control {
  # Avoid watching memory, since the load-controller stops jobs on local runs
  memory-threshold-in-mb = 1
}

%(database)s

backend {
  providers {
    Local {
      config {
        # Avoid overscheduling jobs for local runs by limiting concurrent jobs
        # Longer term would like to keep these within defined core window
        concurrent-job-limit = 1
        runtime-attributes = \"\"\"
        Int? cpu
        Int? memory_mb
        %(docker_attrs)s
        %(cwl_attrs)s
        \"\"\"
        %(submit_docker)s
        %(filesystem)s
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
        %(docker_attrs)s
        %(cwl_attrs)s
        \"\"\"
        submit = \"\"\"
            sbatch -J ${job_name} -D ${cwd} -o ${out} -e ${err} -t ${timelimit} -p ${queue} \
            ${"--cpus-per-task=" + cpu} --mem=${memory_mb} ${account} \
            --wrap "/usr/bin/env bash ${script}"
        \"\"\"
        kill = "scancel ${job_id}"
        check-alive = "squeue -j ${job_id}"
        job-id-regex = "Submitted batch job (\\\\d+).*"
        %(filesystem)s
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
        %(docker_attrs)s
        %(cwl_attrs)s
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
        %(filesystem)s
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
        %(docker_attrs)s
        %(cwl_attrs)s
        \"\"\"
        submit = \"\"\"
        qsub -V -l wd -N ${job_name} -o ${out} -e ${err} -q ${queue} -l walltime=${walltime} \
        %(cpu_and_mem)s \
        -- /usr/bin/env bash ${script}
        \"\"\"
        kill = "qdel ${job_id}"
        check-alive = "qstat -j ${job_id}"
        job-id-regex = "(\\\\d+).*"
        %(filesystem)s
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
        %(docker_attrs)s
        %(cwl_attrs)s
        \"\"\"
        submit = \"\"\"
        qsub -V -d ${cwd} -N ${job_name} -o ${out} -e ${err} -q ${queue} \
        -l nodes=1:ppn=${cpu} -l mem=${memory_mb}mb -l walltime=${walltime} \
        ${script}
        \"\"\"
        kill = "qdel ${job_id}"
        check-alive = "qstat ${job_id}"
        job-id-regex = "(\\\\d+).*"
        %(filesystem)s
      }
    }
""",
"htcondor": """
    HTCONDOR {
      actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
      config {
        runtime-attributes = \"\"\"
          Int cpu = 1
          Float memory_mb = 512.0
          Float disk_kb = 256000.0
          String? nativeSpecs
          %(docker_attrs)s
          %(cwl_attrs)s
        \"\"\"
        submit = \"\"\"
          chmod 755 ${script}
          cat > ${cwd}/execution/submitFile <<EOF
          Iwd=${cwd}/execution
          requirements=${nativeSpecs}
          leave_in_queue=true
          request_memory=${memory_mb}
          request_disk=${disk_kb}
          error=${err}
          output=${out}
          log_xml=true
          request_cpus=${cpu}
          executable=${script}
          log=${cwd}/execution/execution.log
          description=${job_name}
          getenv=true
          queue
          EOF
          condor_submit ${cwd}/execution/submitFile
        \"\"\"
        # submit-docker = \"\"\"
        #   chmod 755 ${script}
        #   cat > ${cwd}/execution/dockerScript <<EOF
        #   #!/bin/bash
        #   docker run --rm -i -v ${cwd}:${docker_cwd} ${docker} /bin/bash ${script}
        #   EOF
        #   chmod 755 ${cwd}/execution/dockerScript
        #   cat > ${cwd}/execution/submitFile <<EOF
        #   Iwd=${cwd}/execution
        #   requirements=${nativeSpecs}
        #   leave_in_queue=true
        #   request_memory=${memory_mb}
        #   request_disk=${disk_kb}
        #   error=${cwd}/execution/stderr
        #   output=${cwd}/execution/stdout
        #   log_xml=true
        #   request_cpus=${cpu}
        #   executable=${cwd}/execution/dockerScript
        #   log=${cwd}/execution/execution.log
        #   queue
        #   EOF
        #   condor_submit ${cwd}/execution/submitFile
        # \"\"\"
        kill = "condor_rm ${job_id}"
        check-alive = "condor_q ${job_id}"
        job-id-regex = "(?sm).*cluster (\\\\d+)..*"
        %(filesystem)s
      }
    }
"""
}
