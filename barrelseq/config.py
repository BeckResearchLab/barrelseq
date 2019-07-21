import argparse
import os
import types

import yaml

from barrelseq import sample

CONFIG_ARGS_KEY_EXCLUDED = ['config_file', 'command', 'config_command', 'func']
REQUIRED_CONFIG_OPTS = ['project_dir', 'project_name', 'reference_fasta_path',
        'reference_gff_path', 'reference_gb_path', 'reference_name', 'bwa_path', 
        'samtools_path', 'htseq_count_path', 'R_path', 'pair_ended'
        ]
WORKSPACE_DIR = 'workspace'

def create(args):
    """Take ``argparse.Namespace``, clean it up and save it as YAML.

    This function removes some keys from the argparse Namespace object
    returned by the ``parser.parse_args`` function and saves the cleaned
    dictionary to a YAML file.  The list of keys removed from the dict
    is a constant to this file and includes argparse state switches and
    open file streams, like the config file stream itself.

    Args:
        ``argparse.Namespace``` args - result from ``argparse.parse_args``

    Side-effects:
        Creates or overwrites a file on disk, the configuration file.
    """
    # start by converting the argparse Namespace object to a dictionary
    config_dict = vars(args).copy()
    # remove entries for non-configuration items
    for key in CONFIG_ARGS_KEY_EXCLUDED:
        config_dict.pop(key, None)

    save(args, config_dict)
    return


def save(args, config_dict):
    # rewind to the top of the config file
    args.config_file.seek(0)
    yaml.dump(config_dict, args.config_file, default_flow_style=False)
    args.config_file.truncate()
    args.config_file.close()
    return


def edit(args):
    """Updates the configuration file from ``argparse.Namespace``.

    This takes an ``argparse.Namespace`` object from parsing the
    configuration create and edit sub commands and cleans it up
    and then uses it to update the config read from a file.  Keys
    for entries from argparse that don't contain true config info
    are removed and followed by entries with a value of None. The
    remaining keys are updated in the config dict read from the
    existing config file.  The results are saved back to the confi
    file.

    Args:
        ``argparse.Namespace``` args - result from ``argparse.parse_args``

    Side-effects:
        Creates or overwrites a file on disk, the configuration file.
    """
    # convert Namespace to dict and remove non-config items
    edited_config_dict = vars(args).copy()
    for key in CONFIG_ARGS_KEY_EXCLUDED:
        edited_config_dict.pop(key, None)
    # drop elements with None values
    edited_config_dict = {k: v for k, v in edited_config_dict.items() if v}
    config = load(args)
    config_dict = vars(config)
    config_dict.update(edited_config_dict)
    save(args, config_dict)
    return


def validate(args):
    """Configuration file validation and load.
    """
    config = load(args)
    config_dict = vars(config)
    for option in REQUIRED_CONFIG_OPTS:
        if option not in config_dict:
            raise RuntimeError('required config option {} is missing from '
                    'the config file {}'.format(option, args.config_file))
    if not os.path.isdir(config.project_dir):
        raise RuntimeError('project directory {} does not'
                'exist'.format(config.project_dir))
    config_dict['workspace_dir'] = os.path.join(config.project_dir, 
                                                WORKSPACE_DIR)
    if not os.path.isdir(config.workspace_dir):
        try:
            print('WARNING: workspace directory {} does not exist\n\t'
                    'it will be created'.format(config.workspace_dir))
            os.mkdir(config.workspace_dir)
        except OSError as e:
            print('ERROR: unable to make workspace '
                    'directory {}'.format(config.workspace_dir))
            raise(e)
        if not os.path.isdir(config.workspace_dir):
            print('ERROR: workspace directory {} is not a '
                'directory'.format(config.workspace_dir))
            raise RuntimeError('unable to open workspace directory')
    else:
        pass
    sample_file = os.path.join(config.project_dir, sample.SAMPLE_INFO_FILE)
    if not os.path.isfile(sample_file):
        print('WARNING: unable to find existing sample info file\n\t'
                'if you are just getting started you can ignore this')
    for path2file in [x for x in REQUIRED_CONFIG_OPTS if x.endswith('_path')]:
        if not os.path.isfile(config_dict[path2file]):
            raise RuntimeError('file {} specified by the config option {} '
                'in the config file does not '
                'exist'.format(config_dict[path2file], path2file))
    print('config file validation succeeded')
    return config


def load(args):
    # as a defensive measure, rewind to start of file before loading
    args.config_file.seek(0)
    config_dict = yaml.safe_load(args.config_file)
    config = types.SimpleNamespace(**config_dict)
    return config
