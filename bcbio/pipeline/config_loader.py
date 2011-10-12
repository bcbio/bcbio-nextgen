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
log_dir:
store_dir:
store_host:
analysis:
	config_file:
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
	""" Goes through two levels of the given yaml  post processing config file, 
	performs replacement of environment variables and returns a dictionary 
	representing the yaml file.
	"""
	with open(config_file) as in_handle:
		config = yaml.load(in_handle)
	
	for field, setting in config.items():
		try:
			config[field] = expand_path(setting)
		except AttributeError:
			pass

		try:
			for sub_field, sub_setting in config[field].items():
				config[field][sub_field] = expand_path(sub_setting)
		except AttributeError:
			pass

	return config

def expand_path(path):
	""" Combines os.path.expandvars with replacing ~ with $HOME.
	"""
	try:
		new_path = os.path.expandvars(path.replace("~", "$HOME"))
		return new_path
	except AttributeError:
		raise AttributeError("Not a path string")