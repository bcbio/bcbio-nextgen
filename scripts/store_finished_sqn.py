"""Transfer raw files from finished NGS runs for backup and storage.

Usage:
    store_finished_sqn.py <Galaxy config file> <Post-processing config file>
"""
import os
import sys
import ConfigParser
import json
import contextlib

import yaml
from amqplib import client_0_8 as amqp
import fabric.api as fabric
import fabric.contrib.files as fabric_files

def main(galaxy_config, processing_config):
    amqp_config = _read_amqp_config(galaxy_config)
    with open(processing_config) as in_handle:
        config = yaml.load(in_handle)
    store_tag = config["msg_store_tag"]
    handlers = [(store_tag, store_handler(config, store_tag))]
    message_reader(handlers, amqp_config)

def copy_for_storage(remote_info, config):
    """Securely copy files from remote directory to the storage server.

    This requires ssh public keys to be setup so that no password entry
    is necessary, Fabric is used to manage setting up copies on the remote
    storage server.
    """
    print remote_info
    base_dir = config["store_dir"]
    fabric.env.host_string = "%s@%s" % (config["store_user"], config["store_host"])
    fc_dir = os.path.join(base_dir, os.path.basename(remote_info['directory']))
    if not fabric_files.exists(fc_dir):
        fabric.run("mkdir %s" % fc_dir)
    for fcopy in remote_info['to_copy']:
        target_loc = os.path.join(fc_dir, fcopy)
        if not fabric_files.exists(target_loc):
            target_dir = os.path.dirname(target_loc)
            if not fabric_files.exists(target_dir):
                fabric.run("mkdir -p %s" % target_dir)
            cl = ["scp", "-r", "%s@%s:%s/%s" % (
                  remote_info["user"], remote_info["hostname"], remote_info["directory"],
                  fcopy), target_loc]
            fabric.run(" ".join(cl))

def store_handler(config, tag_name):
    def receive_msg(msg):
        if msg.properties['application_headers'].get('msg_type') == tag_name:
            copy_for_storage(json.loads(msg.body), config)
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

if __name__ == "__main__":
    main(*sys.argv[1:])
