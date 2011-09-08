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
	
	for field, setting in config.items():
		config[field] = os.path.expandvars(setting)

		for sub_field, sub_setting in config[field].items():
			config[field][sub_field] = os.path.expandvars(sub_setting)	

	return config
