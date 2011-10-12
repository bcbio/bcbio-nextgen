"""Loads configurations from .yaml files and expands
environment variables

The configuration yaml has the structure

galaxy_config: 
program:
	program1:
	program2:
algorithm:
	setting1:
	setting2:
analysis:
	towig_script:
distributed:
	rabbitmq_vhost:
custom_algorithms:
	setting1:
	setting2:

galaxy_config, program and analysis supports 
environment variables.
"""

import os
import yaml

def load_config(config_file):
	with open(config_file) as in_handle:
		config = yaml.load(in_handle)

	config['galaxy_config'] = os.path.expandvars(config['galaxy_config'])
	for program, setting in config['program'].items():
		config['program'][program] = os.path.expandvars(setting)

	config['analysis']['towig_script'] = os.path.expandvars(config['analysis']['towig_script'])

	return config