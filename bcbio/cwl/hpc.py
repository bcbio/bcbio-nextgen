"""Setup configurations for running on HPC clusters with CWL.

Contains support for setting up configuration inputs for Cromwell.
"""
import json
import os

def create_cromwell_config(args, work_dir, sample_file):
    """Prepare a cromwell configuration within the current working directory.
    """
    docker_attrs = ["String? docker", "String? docker_user"]
    cwl_attrs = ["Int? cpuMin", "Int? cpuMax", "Int? memoryMin", "Int? memoryMax", "String? outDirMin",
                 "String? outDirMax", "String? tmpDirMin", "String? tmpDirMax"]
    out_file = os.path.join(work_dir, "bcbio-cromwell.conf")
    run_config = _load_custom_config(args.runconfig) if args.runconfig else {}
    # Avoid overscheduling jobs for local runs by limiting concurrent jobs
    # Longer term would like to keep these within defined core window
    joblimit = args.joblimit
    if joblimit == 0 and not args.scheduler:
        joblimit = 1
    file_types = _get_filesystem_types(args, sample_file)
    std_args = {"docker_attrs": "" if args.no_container else "\n        ".join(docker_attrs),
                "submit_docker": 'submit-docker: ""' if args.no_container else "",
                "joblimit": "concurrent-job-limit = %s" % (joblimit) if joblimit > 0 else "",
                "cwl_attrs": "\n        ".join(cwl_attrs),
                "filesystem": _get_filesystem_config(file_types),
                "database": run_config.get("database", DATABASE_CONFIG % {"work_dir": work_dir})}
    cl_args, conf_args, scheduler, cloud_type = _args_to_cromwell(args)
    std_args["engine"] = _get_engine_filesystem_config(file_types, args, conf_args)
    conf_args.update(std_args)
    main_config = {"hpc": (HPC_CONFIGS[scheduler] % conf_args) if scheduler else "",
                   "cloud": (CLOUD_CONFIGS[cloud_type] % conf_args) if cloud_type else "",
                   "work_dir": work_dir}
    main_config.update(std_args)
    # Local run always seems to need docker set because of submit-docker in default configuration
    # Can we unset submit-docker based on configuration so it doesn't inherit?
    # main_config["docker_attrs"] = "\n        ".join(docker_attrs)
    with open(out_file, "w") as out_handle:
        out_handle.write(CROMWELL_CONFIG % main_config)
    return out_file

def _get_file_paths(cur):
    """Retrieve a list of file paths, recursively traversing the
    """
    out = []
    if isinstance(cur, (list, tuple)):
        for x in cur:
            new = _get_file_paths(x)
            if new:
                out.extend(new)
    elif isinstance(cur, dict):
        if "class" in cur:
            out.append(cur["path"])
        else:
            for k, v in cur.items():
                new = _get_file_paths(v)
                if new:
                    out.extend(new)
    return out

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
    cl_args, conf_args, scheduler, cloud = _args_to_cromwell(args)
    return cl_args

def _args_to_cromwell(args):
    """Convert input arguments into cromwell inputs for config and command line.
    """
    default_config = {"slurm": {"timelimit": "1-00:00", "account": ""},
                      "sge": {"memtype": "mem_free", "pename": "smp"},
                      "lsf": {"walltime": "24:00", "account": ""},
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
    cloud_type = None
    if args.cloud_project:
        if args.cloud_root and args.cloud_root.startswith("gs:"):
            cloud_type = "PAPI"
            cloud_root = args.cloud_root
            cloud_region = None
        elif ((args.cloud_root and args.cloud_root.startswith("s3:")) or
              (args.cloud_project and args.cloud_project.startswith("arn:"))):
            cloud_type = "AWSBATCH"
            cloud_root = args.cloud_root
            if not cloud_root.startswith("s3://"):
                cloud_root = "s3://%s" % cloud_root
            # split region from input Amazon Resource Name, ie arn:aws:batch:us-east-1:
            cloud_region = args.cloud_project.split(":")[3]
        else:
            raise ValueError("Unexpected inputs for Cromwell Cloud support: %s %s" %
                             (args.cloud_project, args.cloud_root))
        config = {"cloud_project": args.cloud_project, "cloud_root": cloud_root, "cloud_region": cloud_region}
        cl.append("-Dbackend.default=%s" % cloud_type)
    return cl, config, args.scheduler, cloud_type

def _get_filesystem_types(args, sample_file):
    """Retrieve the types of inputs and staging based on sample JSON and arguments.
    """
    out = set([])
    ext = "" if args.no_container else "_container"
    with open(sample_file) as in_handle:
        for f in _get_file_paths(json.load(in_handle)):
            if f.startswith("gs:"):
                out.add("gcp%s" % ext)
            elif f.startswith("s3:"):
                out.add("s3%s" % ext)
            elif f.startswith(("https:", "http:")):
                out.add("http%s" % ext)
            else:
                out.add("local%s" % ext)
    return out

def _get_filesystem_config(file_types):
    """Retrieve filesystem configuration, including support for specified file types.
    """
    out = "     filesystems {\n"
    for file_type in sorted(list(file_types)):
        if file_type in _FILESYSTEM_CONFIG:
            out += _FILESYSTEM_CONFIG[file_type]
    out += "      }\n"
    return out


_FILESYSTEM_CONFIG = {
  "gcp": """
        gcs {
          auth = "gcp-auth"
          caching {
            duplication-strategy = "reference"
          }
        }
  """,
  "gcp_container": """
        gcs {
          auth = "gcp-auth"
          caching {
            duplication-strategy = "copy"
          }
        }
  """,
  "http": """
        http { }
  """,
  "http_container": """
        http { }
  """,
  "local": """
        local {
          localization: ["soft-link"]
          caching {
            duplication-strategy: ["soft-link"]
            hashing-strategy: "path"
          }
        }
""",
  "local_container": """
        local {
          localization: ["hard-link", "copy"]
          caching {
            duplication-strategy: ["hard-link", "copy"]
            hashing-strategy: "file"
          }
        }
"""
}

DATABASE_CONFIG = """
database {
  profile = "slick.jdbc.HsqldbProfile$"
  db {
    driver = "org.hsqldb.jdbcDriver"
    url = "jdbc:hsqldb:file:%(work_dir)s/persist/metadata;shutdown=false;hsqldb.tx=mvcc"
    connectionTimeout = 200000
  }
}
"""

def _get_engine_filesystem_config(file_types, args, conf_args):
    """Retriever authorization and engine filesystem configuration.
    """
    file_types = [x.replace("_container", "") for x in list(file_types)]
    out = ""
    if "gcp" in file_types:
        out += _AUTH_CONFIG_GOOGLE
    if "s3" in file_types:
        out += _AUTH_CONFIG_AWS % conf_args["cloud_region"]
    if "gcp" in file_types or "http" in file_types or "s3" in file_types:
        out += "engine {\n"
        out += "  filesystems {\n"
        if "gcp" in file_types:
            out += '    gcs {\n'
            out += '      auth = "gcp-auth"\n'
            if args.cloud_project:
                out += '      project = "%s"\n' % args.cloud_project
            out += '    }\n'
        if "http" in file_types:
            out += '    http {}\n'
        if "s3" in file_types:
            out += '    s3 { auth = "default" }'
        out += "  }\n"
        out += "}\n"

    return out

_AUTH_CONFIG_AWS = """
aws {
  application-name = "cromwell"
  auths = [{
      name = "default"
      scheme = "default"
  }]
  region = "%s"
}
"""

_AUTH_CONFIG_GOOGLE = """
google {
  application-name = "cromwell"
  auths = [
    {
      name = "gcp-auth"
      scheme = "service_account"
      json-file = ${?GOOGLE_APPLICATION_CREDENTIALS}
    }
  ]
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

cwltool-runner {
  # Use external cwltool to avoid slow runtimes with java embedded pre-processing
  class = "cwl.CwltoolProcess"
}

%(database)s

%(engine)s

backend {
  providers {
    Local {
      config {
        %(joblimit)s
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
%(cloud)s
  }
}
"""

HPC_CONFIGS = {
"slurm": """
    SLURM {
      actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
      config {
        %(joblimit)s
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
        %(joblimit)s
        runtime-attributes = \"\"\"
        Int cpu = 1
        Int memory_mb = 2048
        String queue = "%(queue)s"
        String pename = "%(pename)s"
        String memtype = "%(memtype)s"
        %(docker_attrs)s
        %(cwl_attrs)s
        \"\"\"
        submit = \"\"\"
        qsub -V -w w -j y -N ${job_name} -wd ${cwd} \
        -o ${out} -e ${err} -q ${queue} \
        -pe ${pename} ${cpu} -l ${memtype}=${memory_mb}m \
        ${script}
        \"\"\"
        kill = "qdel ${job_id}"
        check-alive = "qstat -j ${job_id}"
        job-id-regex = "(\\\\d+)"
        %(filesystem)s
      }
    }
""",
"lsf": """
    LSF {
      actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
      config {
        %(joblimit)s
        runtime-attributes = \"\"\"
        Int cpu = 1
        Int memory_gb = 2
        String queue = "%(queue)s"
        String account = "%(account)s"
        String walltime = "%(walltime)s"
        %(docker_attrs)s
        %(cwl_attrs)s
        \"\"\"
        submit = \"\"\"
        bsub -J ${job_name} -cwd ${cwd} -o ${out} -e ${err} -q ${queue} -W ${walltime} \
        -R rusage[mem=${memory_gb}] -n ${cpu} \
        /usr/bin/env bash ${script}
        \"\"\"
        kill = "bkill ${job_id}"
        check-alive = "bjobs ${job_id}"
        job-id-regex = "Job <(\\\\d+)>.*"
        %(filesystem)s
      }
    }

""",
"pbspro": """
    PBSPRO {
      actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
      config {
        %(joblimit)s
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
        %(joblimit)s
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
        %(joblimit)s
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

CLOUD_CONFIGS = {
"PAPI": """
    PAPI {
      actor-factory = "cromwell.backend.google.pipelines.v2alpha1.PipelinesApiLifecycleActorFactory"
      config {
        project = "%(cloud_project)s"
        root = "%(cloud_root)s"

        genomics {
          auth = "gcp-auth"
          endpoint-url = "https://genomics.googleapis.com/"
        }

        filesystems {
          gcs {
            auth = "gcp-auth"
            project = "%(cloud_project)s"
          }
        }
      }
    }
""",
"AWSBATCH": """
    AWSBATCH {
          actor-factory = "cromwell.backend.impl.aws.AwsBatchBackendLifecycleActorFactory"
          config {
            root = "%(cloud_root)s/cromwell-execution"
            auth = "default"

            numSubmitAttempts = 3
            numCreateDefinitionAttempts = 3

            default-runtime-attributes {
              queueArn: "%(cloud_project)s"
            }
            filesystems {
              s3 {
                auth = "default"
              }
            }
          }
        }
"""
}
