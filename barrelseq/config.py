import argparse
import types

import yaml


CONFIG_ARGS_KEY_EXCLUDED = ['config_file', 'command', 'config_command', 'func']

def create(args):
    """Take an argsparse Namespace object, clean it up and save it as YAML.

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
    print(args)
    return


def load(args):
    # as a defensive measure, rewind to start of file before loading
    args.config_file.seek(0)
    config_dict = yaml.safe_load(args.config_file)
    config = types.SimpleNamespace(**config_dict)
    return config
