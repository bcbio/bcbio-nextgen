# Parallel execution
The pipeline runs in parallel in two different ways:
* multiple cores -- Analyses will run in parallel using multiple cores on
a single machine. This requires only the `multiprocessing` Python library,
included by default with most Python installations.
* parallel messaging -- This allows scaling beyond the cores available on a single machine, and requires multiple machines with a shared filesystem like standard cluster environments. Machine to machine communication occurs via messaging, using the [IPython parallel](https://ipython.readthedocs.io/en/stable/) framework.

## Tuning core and memory usage
bcbio has two ways to specify core usage, helping provide options for parallelizing different types of processes:
* Total available cores: specified with `-n` on the commandline, this tells bcbio how many total cores to use. This applies either to a local multicore run or a distributed job.
* Maximum cores to use for multicore processing of individual jobs. You specify this in the `resource` section of either a sample YAML file or `bcbio_system.yaml`. Ideally you specify this in the `default` section (along with memory usage). For example, this would specify that processes using multiple cores can get up to 16 cores with 2GB of memory per core:
```yaml
resources:
  default:
    memory: 2G
    cores: 16
    jvm_opts: ["-Xms750m", "-Xmx2000m"]
```
bcbio uses these settings, along with memory requests, to determine how to partition jobs. For example, if you had `-n 32` and `cores: 16` for a run on a single 32 core machine, this would run two simultaneous bwa mapping jobs using 16 cores each.

Memory specifications (both in `memory` and `jvm_opts`) are per-core. bcbio takes care of adjusting this memory to match the cores used. In the example above, if bcbio was running a 16 core java process, it would use 32GB of memory for the JVM, adjusting `Xmx` and `Xms` to match cores used. Internally bcbio looks at the memory and CPU usage on a machine and matches your configuration options to the available system resources. It will scale down core requests if memory is limiting, avoiding over-scheduling resources during the run. You ideally want to set both `memory` and `jvm_opts` to match the average memory per core on the run machine and adjust upwards if this does not provide enough memory for some processes during the run.

For single machine runs with a small number of samples, you generally want to set `cores` close to or equal the number of total cores you're allocating to the job with `-n`. This will allow individual samples to process as fast as possible and take advantage of multicore software.

For distributed jobs, you want to set `cores` to match the available cores on a single node in your cluster, then use `-n` as a multiple of this to determine how many nodes to spin up. For example, `cores: 16` and `-n 64` would try to make four 16 core machines available for analysis.

## Multiple cores

Running using multiple cores only requires setting the `-n` command line flag:
```shell
bcbio_nextgen.py bcbio_sample.yaml -t local -n 12
```

## IPython parallel

[IPython parallel](https://ipython.readthedocs.io/en/stable/) provides a distributed framework for performing parallel computation in standard cluster environments. The bcbio-nextgen setup script installs both IPython and [pyzmq](https://github.com/zeromq/pyzmq), which provides Python bindings for the [ZeroMQ](https://zeromq.org/) messaging library. The only additional requirement is that the work directory where you run the analysis is accessible to all processing nodes. This is typically accomplished with a distributed file system like [NFS](https://en.wikipedia.org/wiki/Network_File_System), [Gluster](https://www.gluster.org/) or [Lustre](http://wiki.lustre.org/Main_Page).

Run an analysis using ipython for parallel execution:
```shell
bcbio_nextgen.py bcbio_sample.yaml -t ipython -n 12 -s lsf -q queue
```
The `-s` flag specifies a type of scheduler to use `(lsf, sge, torque, slurm, pbspro)`.

The `-q` flag specifies the queue to submit jobs to.

The `-n` flag defines the total number of cores to use on the cluster during processing. The framework will select the appropriate number of cores and type of cluster (single core versus multi-core) to use based on the pipeline stage (see the `Parallel` section in the internals documentation for more details). For multiple core steps, the number of cores to use for programs like `bwa`, `novoalign` and `gatk` comes from the `Resources` section of the configuration. Ensure the `cores` specification matches the physical cores available on machines in your cluster, and the pipeline will divide the total cores specified by `-n` into the appropriate number of multicore jobs to run.

The pipeline default parameters assume a system with minimal time to obtain processing cores and consistent file system accessibility. These defaults allow the system to fail fast in the case of cluster issues which need diagnosis. For running on shared systems with high resource usage and potential failures due to intermittent cluster issues, there are turning parameters that increase resiliency. The `--timeout` flag specifies the numbers of minutes to wait for a cluster to start up before timing out. This defaults to 15 minutes. The `--retries` flag specify the number of times to retry a job on failure. In systems with transient distributed file system hiccups like lock errors or disk availability, this will provide recoverability at the cost of resubmitting jobs that may have failed for reproducible reasons.

Finally, the `-r resources` flag specifies resource options to pass along to the underlying queue scheduler. This currently supports SGE's `-l` parameter, Torque's `-l` parameter and LSF and SLURM native flags. This allows specification or resources to the scheduler (see the [qsub man page](http://gridscheduler.sourceforge.net/htmlman/htmlman1/qsub.html)). You may specify multiple resources, so `-r mem=4g -r ct=01:40:00` translates to `-l mem=4g -l ct=01:40:00` when passed to `qsub` or `-r "account=a2010002" -r "timelimit=04:00:00"` when using SLURM, for instance. SLURM and Torque support specification of an account parameter with `-r account=your_name`, which IPython transfers into `-A`.

SGE supports special parameters passed using resources to help handle the heterogeneity of possible setups.

Specify an [SGE parallel environment](https://docs.oracle.com/cd/E19957-01/820-0698/6ncdvjcmd/index.html) that supports using multiple cores on a single node with `-r pename=your_pe`. Since this setup is system specific it is hard to write general code for find a suitable environment. Specifically, when there are multiple usable parallel environments, it will select the first one which may not be correct. Manually specifying it with a `pename=` flag to resources will ensure correct selection of the right environment. If you're administering a grid engine cluster and not sure how to set this up you'd typically want a `smp` queue using `allocation_rule: $pe_slots` like in this [example pename configuration](https://github.com/WGLab/biocluster/blob/431a05f6dfd532205aacfc7477ac740b0e7b2a0a/03%20System%20customization.md#setting-up-parallel-environment) or [smp template](https://gist.github.com/dan-blanchard/6586533#file-smp_template).

SGE has other specific flags you may want to tune, depending on your setup. To specify an advanced reservation with the `-ar` flag, use `-r ar=ar_id`. To specify an alternative memory management model instead of `mem_free` use `-r memtype=approach`. It is further recommended to configure `mem_free` (or any other chosen memory management model) as a consumable, requestable resource in SGE to prevent overfilling hosts that do not have sufficient memory per slot. This can be done in two steps. First, launch `qmon` as an admin, select `Complex Configuration` in qmon, click on `mem_free`, under the `Consumable` dialog select `JOB` (instead of `YES` or `NO`) and finally click `Modify` for the changes to take effect. Secondly, for each host in the queue, configure`mem_free` as a complex value. If a host called `myngshost` has 128GB of RAM, the corresponding command would be `qconf -mattr exechost complex_values mem_free=128G myngshost`.

There are also special `-r`resources parameters to support pipeline configuration:
* `-r conmem=4` -- Specify the memory for the controller process, in GB. This currently applies to SLURM processing and defaults to 4GB.
* `-r minconcores=2` -- The minimum number of cores to use for the controller process. The controller one works on a single core but this can help in queues where you can only specify multicore jobs.
* `-r mincores=16` -- Specify the minimum number of cores to batch together for parallel single core processes like variant calling. This will run multiple processes together under a single submission to allow sharing of resources like memory, which is helpful when a small percentage of the time a process like variant calling will use a lot of memory. By default, bcbio will calculate `mincores` based on specifications for multicore calling so this doesn't normally require a user to set.

## Troubleshooting

### Diagnosing job failures

Parallel jobs can often terminate with rather generic failures like any of the following:
* `joblib/parallel.py, ... TypeError: init() takes at least 3 arguments (2 given)`
* `Multiprocessing exception:`
* `CalledProcessError: Command '<command line that failed>'`

These errors unfortunately don't help diagnose the problem, and you'll likely see the actual error triggering this generic exception earlier in the run. This error can often be hard to find due to parallelization. If you run into a confusing failure like this, the best approach is to re-run with a single core:
```shell
bcbio_nextgen.py your_input.yaml -n 1
```
which should produce a more helpful debug message right above the failure. It's also worth re-trying the failed command line outside of bcbio to look for errors. You can find the failing command by cross-referencing the error message with command lines in `log/bcbio-nextgen-commands.log`. You may have to change temporary directories (`tx/tmp**`) in some of the job outputs. Reproducing the error outside of bcbio is a good first step to diagnosing and fixing the underlying issue.

### No parallelization where expected

This may occur if the current execution is a re-run of a previous project:
* Files in `checkpoints_parallel/*.done` tell bcbio not to parallelize already executed pipeline tasks. This makes restarts faster by avoiding re-starting a cluster (when using distributed runs) for finished stages. If that behaviour is not desired for a task, removing the checkpoint file will get things parallelizing again.
* If the processing of a task is nearly finished the last jobs of this task will be running and bcbio will wait for those to finish.

### IPython parallelization problems

Networking problems on clusters can prevent the IPython parallelization framework from working properly. Be sure that the compute nodes on your cluster are aware of IP addresses that they can use to communicate with each other (usually these will be local IP addresses). Running:
```shell
python -c 'import socket; print socket.gethostbyname(socket.gethostname())'
```
Should return such an IP address (as opposed to localhost). This can be fixed by adding an entry to the hosts file.
The line:
```
host-ip hostname
```
where `host-ip` is replaced by the actual IP address of the machine and `hostname` by the machine's own hostname, should be aded to `/etc/hosts` on each compute node. This will probably involve contacting your local cluster administrator.

## Memory management

The memory information specified in the system configuration `Resources` enables scheduling of memory intensive processes. The values are specified on a *memory-per-core* basis and thus bcbio-nextgen handles memory scheduling by:
* [Determining available cores and memory per machine](#determining-available-cores-and-memory-per-machine)
* Calculating the memory and core usage. The system configuration `Resources` contains the expected core and memory usage of external programs.
* Adjusting the specified number of total cores to avoid over-scheduling memory. This allows running programs with more than the available memory per core without getting out of memory system   errors.
* Passing total memory usage along to schedulers. The SLURM, SGE, Torque and PBSPro schedulers use this information to allocate memory to processes, avoiding issues with other scheduled programs using available memory on a shared machine.

As a result of these calculations, the cores used during processing will not always correspond to the maximum cores provided in the input `-n` parameter. The goal is rather to intelligently maximize cores and memory while staying within system resources. Note that memory specifications are for a single core, and the pipeline takes care of adjusting this to actual cores used during processing.

## Determining available cores and memory per machine

bcbio automatically tries to determine the total available memory and cores per machine for balancing resource usage. For multicore runs, it retrieves total memory from the current machine. For parallel runs, it spawns a job on the queue and extracts the system information from that machine. This expects a homogeneous set of machines within a cluster queue. You can see the determined cores and total memory in `provenance/system-ipython-queue.yaml`.

For heterogeneous clusters or other cases where bcbio does not correctly identify available system resources, you can manually set the machine cores and total memory in the`resource` section of either a sample YAML file or `bcbio_system.yaml`:
```yaml
resources:
  machine:
    memory: 48.0
    cores: 16
```
The memory usage is total available on the machine in GB, so this specifies that individual machines have 48GB of total memory and 16 cores.

## Tuning systems for scale

bcbio-nextgen scales out on clusters including hundreds of cores and is stress tested on systems with 1000 simultaneous processes. Scaling up often requires system specific tuning to handle simultaneous processes. This section collects useful tips and tricks for managing scaling issues.

### Open file handles

A common failure mode is having too many open file handles. This error report can come from the IPython infrastructure logs as ZeroMQ attempts to open sockets, or from the processing logs as third party software gets file handles. You can check your available file handles with`ulimit -a | grep open`. Setting open file handle limits is open system and cluster specific and below are tips for specific setups.

In addition to open file handle limits (`ulimit -n`) large processes may also run into issues with available max user processes (`ulimit -u`). Some systems set a low soft limit (`ulimit -Su`) like 1024 but a higher hard limit (`ulimit -Hu`), allowing adjustment without root privileges. The IPython controllers and engines do this automatically, but the main `bcbio_nextgen.py` driver process cannot. If this scheduler puts this process on the same node as worker processes, you may run into open file handle limits due to work happening on the workers. To fix this, manually set `ulimit
-u a_high_number` as part of the submission process for the main process.

For a Ubuntu system, edit `/etc/security/limits.conf` to set the soft and hard `nofile`descriptors, and edit `/etc/pam.d/common-session` to add `pam_limits.so`. See [this blog post](https://viewsby.wordpress.com/2013/01/29/ubuntu-increase-number-of-open-files/) for more details.

For CentOS/RedHat systems, edit `/etc/security/limits.conf` and `/etc/security/limits.d/90-nproc.conf` to [increase maximum open files and user limits](https://ithubinfo.blogspot.com/2013/07/how-to-increase-ulimit-open-file-and.html).

SGE needs configuration at the qmaster level. Invoke `qconf -mconf` from a host with admin privileges, and edit `execd_params`:
```
execd_params                 S_DESCRIPTORS=20000
```

### IO and Network File Systems

bcbio-nextgen makes use of distributed network file systems to manage sharing large files between compute nodes. While we strive to minimize disk-based processing by making use of pipes, the pipeline still has a major IO component. To help manage IO and network bottlenecks, this section contains pointers on deployments and benchmarking. Please contribute your tips and thoughts.

Harvard and Dell: See the 'Distributed File Systems' section of our [post on scaling bcbio-nextgen](https://bcb.io/2013/05/22/scaling-variant-detection-pipelines-for-whole-genome-sequencing-analysis/) for details about the setup within [Harvard FAS Research Computing](https://www.rc.fas.harvard.edu/) and thoughts on scaling and hardware. We also collaborate with Dell to test the pipeline on Dell's Active Infrastructure for Life Sciences. We found the biggest initial factor limiting scaling was network bandwidth between compute and storage nodes.

### Spark
Some GATK tools like recalibration use Apache Spark for parallelization. By default bcbio runs these with multicore parallelization on a single node, to fit in standard cluster and local compute environments. If you have a custom Spark cluster on your system you can use that for GATK by setting up the appropriate configuration in your `Sample or run specific resources`:
```yaml
resources:
  gatk-spark:
    options: [--spark-master, 'spark://your-spark-cluster:6311']
```

## Profiling

Profiling (tracking CPU, memory, IO usage) could help to optimize resource usage of bcbio, especially when running on a server or AWS instance. Sometimes running a bcbio project with 32 cores is just 10% more efficient than with 16 cores, because a particular configuration might have memory or IO related bottlenecks. IO bottlenecks are when you see low CPU utilization and high read/write values.

1. [Install and start sysstat deamon](http://www.leonardoborda.com/blog/how-to-configure-sysstatsar-on-ubuntudebian/).
1. Create a cron job to gather system statistics every minute or two.
1. Before bcbio start, drop system memory caches. Otherwise memory usage statistic might be misleading:
    ```shell
    # become root
    sudo su
    echo 1 > /proc/sys/vm/drop_caches
    ```
1. Record bcbio project start and stop time (_date_)
1. Collect usage statistics:
    ```shell
    # CPU load
    sar -q -s $start -e $end |  awk '{print $1","$4}'  | sed 1d | sed 1d > cpu.csv
    # memory
    sar -r -s $start -e $end | awk '{print $5}' | sed 1d | sed 1d > mem.csv
    # IO
    sar -b -s $start -e $end | awk '{print $5","$6}' | sed 1d | sed 1d > io.csv
    paste -d "," cpu.csv mem.csv io.csv > usage.csv
    ```
    ```
    # Example of usage.csv

    23:07:01,ldavg-1,%memused,bread/s,bwrtn/s
    23:08:01,1.77,3.10,23238.66,20204.10
    23:09:01,4.34,24.28,208650.45,26270.58
    23:10:01,11.09,25.67,0.13,15.46
    23:11:01,13.56,27.00,4.27,21.99
    23:12:01,15.44,29.63,26.52,2749.22
    23:13:01,15.16,29.75,42.93,27.06
    23:14:01,16.94,30.54,205.26,2740.95
    23:15:01,15.76,30.57,28.92,2751.62
    23:16:01,15.77,30.88,6.13,33.59
    ```
    See [sar man page](https://linux.die.net/man/1/sar) for more fields and field definitions.
1. Overlap profiling results with bcbio-nextgen-commands.log to investigate the performance of particular steps.
