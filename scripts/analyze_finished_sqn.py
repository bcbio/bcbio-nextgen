"""Server which listens for finished NGS runs, processing for upload to galaxy.

Usage:
    analyze_finished_sqn.py <Galaxy config file> <Post-processing config file>

Need to configure the RabbitMQ server with:

    rabbitmqctl add_user galaxy password
    rabbitmqctl add_vhost galaxy_messaging_engine
    rabbitmqctl set_permissions -p galaxy_messaging_engine galaxy '.*' '.*' '.*'
"""
import os
import sys
import ConfigParser
import json
import subprocess
import contextlib

import yaml
from amqplib import client_0_8 as amqp

def main(galaxy_config, processing_config):
    amqp_config = _read_amqp_config(galaxy_config)
    with open(processing_config) as in_handle:
        config = yaml.load(in_handle)
    process_tag = config["msg_process_tag"]
    handlers = [(process_tag,
        analysis_handler(config, process_tag, processing_config))]
    message_reader(handlers, amqp_config)

def copy_and_analyze(remote_info, config, config_file):
    """Remote copy an output directory, process it, and upload to Galaxy.
    """
    print remote_info
    fc_dir = _remote_copy(remote_info, config["local_sqn_dir"])
    print fc_dir
    analysis_dir = os.path.join(config["analysis"]["base_dir"],
                                os.path.basename(remote_info["directory"]))
    if not config_file.startswith("/"):
        config_file = os.path.join(os.getcwd(), config_file)
    with _make_and_chdir(analysis_dir):
        cl = [config["analysis"]["process_program"], config_file, fc_dir]
        subprocess.check_call(cl)
    cl = [config["analysis"]["upload_program"], config_file, fc_dir, analysis_dir]
    subprocess.check_call(cl)

def _remote_copy(remote_info, local_sqn_dir):
    """Securely copy files from remote directory to the processing server.

    This requires ssh public keys to be setup so that no password entry
    is necessary.
    """
    fc_dir = os.path.join(local_sqn_dir,
            os.path.basename(remote_info['directory']))
    if not os.path.exists(fc_dir):
        os.makedirs(fc_dir)
    for fcopy in remote_info['to_copy']:
        target_loc = os.path.join(fc_dir, fcopy)
        if not os.path.exists(target_loc):
            target_dir = os.path.dirname(target_loc)
            if not os.path.exists(target_dir):
                os.makedirs(target_dir)
            cl = ["scp", "-r", "%s@%s:%s/%s" % (remote_info["user"],
                      remote_info["hostname"], remote_info["directory"], fcopy),
                  target_loc]
            subprocess.check_call(cl)
    return fc_dir

def analysis_handler(processing_config, tag_name, config_file):
    def receive_msg(msg):
        if msg.properties['application_headers'].get('msg_type') == tag_name:
            copy_and_analyze(json.loads(msg.body), processing_config,
                config_file)
    return receive_msg

def message_reader(handlers, config):
    """Wait for messages with the give tag, passing on to the supplied handler.
    """
    conn = amqp.Connection(host=config['host'] + ":" + config['port'],
                           userid=config['userid'], password=config['password'],
                           virtual_host=config['virtual_host'], insist=False)
    chan = conn.channel()
    for tag_name, handler in handlers:
        chan.queue_declare(queue=tag_name, exclusive=False, auto_delete=False,
                durable=True)
        chan.exchange_declare(exchange=config['exchange'], type="fanout", durable=True,
                auto_delete=False)
        chan.queue_bind(queue=tag_name, exchange=config['exchange'],
                        routing_key=config['routing_key'])
        chan.basic_consume(queue=tag_name, no_ack=True,
                           callback=handler, consumer_tag=tag_name)
    while True:
        chan.wait()
    for (tag_name, _) in handlers:
        chan.basic_cancel(tag_name)
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

@contextlib.contextmanager
def _make_and_chdir(new_dir):
    if not os.path.exists(new_dir):
        os.makedirs(new_dir)
    cur_dir = os.getcwd()
    try:
        os.chdir(new_dir)
        yield
    finally:
        os.chdir(cur_dir)

if __name__ == "__main__":
    main(*sys.argv[1:])
