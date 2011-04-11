"""Script to check for finalized illumina runs and report to messaging server.

This is meant to be run via a cron job on a regular basis, and looks for newly
dumped output directories that are finished and need to be processed.

Usage:
    illumina_finished_msg.py <Galaxy config> <YAML local config>

The Galaxy config needs to have information on the messaging server and queues.
The local config should have the following information:

    msg_process_tag, msg_store_tag: tag names to send messages for processing and
                                    storage
    dump_directories: directories to check for machine output
    msg_db: flat file of output directories that have been reported
"""
import os
import sys
import json
import operator
import ConfigParser
import socket
import glob
import getpass
import subprocess
from optparse import OptionParser

import yaml
from amqplib import client_0_8 as amqp
import logbook

from bcbio.solexa import samplesheet
from bcbio.log import create_log_handler
from bcbio import utils
from bcbio.solexa.flowcell import (get_flowcell_info, get_fastq_dir, get_qseq_dir)

LOG_NAME = os.path.splitext(os.path.basename(__file__))[0]
log = logbook.Logger(LOG_NAME)

def main(galaxy_config, local_config, process_msg=True, store_msg=True,
         qseq=True, fastq=True):
    amqp_config = _read_amqp_config(galaxy_config)
    with open(local_config) as in_handle:
        config = yaml.load(in_handle)
    log_handler = create_log_handler(config, LOG_NAME)
    with log_handler.applicationbound():
        search_for_new(config, amqp_config, process_msg, store_msg, qseq, fastq)

def search_for_new(config, amqp_config, process_msg, store_msg, qseq, fastq):
    """Search for any new directories that have not been reported.
    """
    reported = _read_reported(config["msg_db"])
    for dname in _get_directories(config):
        if os.path.isdir(dname) and dname not in reported:
            if _is_finished_dumping(dname):
                log.info("The instrument has finished dumping on directory %s" % dname)
                _update_reported(config["msg_db"], dname)

                ss_file = samplesheet.run_has_samplesheet(dname, config)
                if ss_file:
                    out_file = os.path.join(dname, "run_info.yaml")
                    log.info("CSV Samplesheet %s found, converting to %s" %
                             (ss_file, out_file))
                    samplesheet.csv2yaml(ss_file, out_file)
                if qseq:
                    log.info("Generating qseq files for %s" % dname)
                    _generate_qseq(get_qseq_dir(dname), config)
                if fastq:
                    log.info("Generating fastq files for %s" % dname)
                    _generate_fastq(dname, config)

                store_files, process_files = _files_to_copy(dname)

                if process_msg:
                    finished_message(config["msg_process_tag"], dname,
                                     process_files, amqp_config)
                if store_msg:
                    finished_message(config["msg_store_tag"], dname,
                                     store_files, amqp_config)

def _generate_fastq(fc_dir, config):
    """Generate fastq files for the current flowcell.
    """
    fc_name, fc_date = get_flowcell_info(fc_dir)
    short_fc_name = "%s_%s" % (fc_date, fc_name)
    fastq_dir = get_fastq_dir(fc_dir)
    basecall_dir = os.path.split(fastq_dir)[0]
    if not fastq_dir == fc_dir and not os.path.exists(fastq_dir):
        log.info("Generating fastq files for %s" % fc_dir)
        with utils.chdir(basecall_dir):
            lanes = sorted(list(set([f.split("_")[1] for f in
                glob.glob("*qseq.txt")])))
            cl = ["solexa_qseq_to_fastq.py", short_fc_name,
                    ",".join(lanes)]
            log.info("Converting qseq to fastq on all lanes.")
            subprocess.check_call(cl)
            log.info("Qseq to fastq conversion completed.")
    return fastq_dir

def _generate_qseq(bc_dir, config):
    """Generate qseq files from illumina bcl files if not present.

    More recent Illumina updates do not produce qseq files. These can be
    generated from bcl, intensity and filter files with tools from
    the offline base caller OLB.
    """
    if not os.path.exists(os.path.join(bc_dir, "finished.txt")):
        log.info("Generating qseq files at %s" % bc_dir)
        bcl2qseq_log = os.path.join(config["log_dir"], "setupBclToQseq.log")
        cmd = os.path.join(config["program"]["olb"], "bin", "setupBclToQseq.py")
        cl = [cmd, "-L", bcl2qseq_log,"-o", bc_dir, "-P", "_pos.txt", "--in-place", "--overwrite"]
        # in OLB version 1.9, the -i flag changed to intensities instead of input
        version_cl = [cmd, "-v"]
        p = subprocess.Popen(version_cl, stdout=subprocess.PIPE)
        (out, _) = p.communicate()
        olb_version = float(out.strip().split()[-1].rsplit(".", 1)[0])
        if olb_version > 1.8:
            cl += ["-b", bc_dir]
        else:
            cl += ["-i", bc_dir, "-p", os.path.split(bc_dir)[0]]
        subprocess.check_call(cl)
        log.info("Qseq files generated.")
        with utils.chdir(bc_dir):
            try:
                processors = config["algorithm"]["num_cores"]
            except KeyError:
                processors = 8
            cl = config["program"].get("olb_make", "make").split() + ["-j", str(processors)]
            subprocess.check_call(cl)

def _is_finished_dumping(directory):
    """Determine if the sequencing directory has all files.

    The final checkpoint file will differ depending if we are a
    single or paired end run.
    """
    to_check = ["Basecalling_Netcopy_complete_SINGLEREAD.txt",
                "Basecalling_Netcopy_complete_READ2.txt",
                "Basecalling_Netcopy_complete_Read3.txt"]

    return reduce(operator.or_,
            [os.path.exists(os.path.join(directory, f)) for f in to_check])

def _files_to_copy(directory):
    """Retrieve files that should be remotely copied.
    """
    with utils.chdir(directory):
        image_redo_files = reduce(operator.add,
                                  [glob.glob("*.params"),
                                   glob.glob("Images/L*/C*"),
                                   ["RunInfo.xml", "runParameters.xml"]])
        qseqs = reduce(operator.add,
                     [glob.glob("Data/Intensities/*.xml"),
                      glob.glob("Data/Intensities/BaseCalls/*qseq.txt"),
                      ])
        reports = reduce(operator.add,
                     [glob.glob("*.xml"),
                      glob.glob("Data/Intensities/BaseCalls/*.xml"),
                      glob.glob("Data/Intensities/BaseCalls/*.xsl"),
                      glob.glob("Data/Intensities/BaseCalls/*.htm"),
                      ["Data/Intensities/BaseCalls/Plots", "Data/reports", 
                       "Data/Status.htm", "Data/Status_Files", "InterOp"]])
        
        logs = reduce(operator.add, [["Logs", "Recipe", "Diag", "Data/RTALogs", "Data/Log.txt"]])
        run_info = glob.glob("run_info.yaml")
        fastq = ["Data/Intensities/BaseCalls/fastq"]
        
    return sorted(image_redo_files + logs + reports + run_info), sorted(reports + fastq + run_info)

def _read_reported(msg_db):
    """Retrieve a list of directories previous reported.
    """
    reported = []
    if os.path.exists(msg_db):
        with open(msg_db) as in_handle:
            for line in in_handle:
                reported.append(line.strip())
    return reported

def _get_directories(config):
    for directory in config["dump_directories"]:
        for dname in sorted(glob.glob(os.path.join(directory, "*A?XX"))):
             if os.path.isdir(dname):
                 yield dname

def _update_reported(msg_db, new_dname):
    """Add a new directory to the database of reported messages.
    """
    with open(msg_db, "a") as out_handle:
        out_handle.write("%s\n" % new_dname)

def finished_message(tag_name, directory, files_to_copy, config):
    """Wait for messages with the give tag, passing on to the supplied handler.
    """
    log.info("Sending finished message to: %s" % tag_name)
    user = getpass.getuser()
    hostname = socket.gethostbyaddr(socket.gethostname())[0]
    data = dict(
            machine_type='illumina',
            hostname=hostname,
            user=user,
            directory=directory,
            to_copy=files_to_copy
            )
    conn = amqp.Connection(host=config['host'] + ":" + config['port'],
                           userid=config['userid'], password=config['password'],
                           virtual_host=config['virtual_host'], insist=False)
    chan = conn.channel()
    chan.queue_declare(queue=tag_name, exclusive=False, auto_delete=False,
            durable=True)
    try:
        chan.exchange_declare(exchange=config['exchange'], type="fanout", durable=True,
                auto_delete=False)
    except amqp.exceptions.AMQPChannelException:
        chan.exchange_delete(exchange=config['exchange'])
        chan.exchange_declare(exchange=config['exchange'], type="fanout", durable=True,
                auto_delete=False)
    msg = amqp.Message(json.dumps(data),
                       content_type='application/json',
                       application_headers={'msg_type': tag_name})
    msg.properties["delivery_mode"] = 2
    chan.basic_publish(msg, exchange=config['exchange'],
                       routing_key=config['routing_key'])
    chan.close()
    conn.close()

def _read_amqp_config(galaxy_config):
    """Read connection information on the RabbitMQ server from Galaxy config.
    """
    config = ConfigParser.ConfigParser()
    config.read(galaxy_config)
    amqp_config = {}
    for option in config.options("galaxy_amqp"):
        amqp_config[option] = config.get("galaxy_amqp", option)
    return amqp_config

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-p", "--noprocess", dest="process_msg",
            action="store_false", default=True)
    parser.add_option("-s", "--nostore", dest="store_msg",
            action="store_false", default=True)
    parser.add_option("-f", "--nofastq", dest="fastq",
            action="store_false", default=True)
    parser.add_option("-q", "--noqseq", dest="qseq",
            action="store_false", default=True)

    (options, args) = parser.parse_args()
    kwargs = dict(process_msg=options.process_msg, store_msg=options.store_msg,
                  fastq=options.fastq, qseq=options.qseq)
    main(*args, **kwargs)
