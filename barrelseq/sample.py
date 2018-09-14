import argparse
import os
import sys

import pandas as pd

from barrelseq import config


SAMPLE_INFO_FILE = 'sample_info.tsv'
SAMPLE_INFO_COLUMNS = ['name', 'group', 'fastq_files', 'description',
            'created', 'last_run', 'last_analyzed', 'md5sums'
        ]


def samples_init_df():
    """Creates a skeleton data frame with correct columns for sample info.
    """
    return pd.DataFrame(index=['name'], columns=SAMPLE_INFO_COLUMNS)


def load(cfg):
    """Load the sample metadata given a config object.

    The filename of the sample info file is assembled using the config
    object. It is attempted to be opened.  On the first try, if it fails
    it will try to create a skeleton file with only column headings.
    On a subsequent failure it prints an error message and exits"

    Args:
        cfg - ``object`` whose attributes are the configuration params

    """
    filename = os.path.join(cfg.project_dir, SAMPLE_INFO_FILE)
    while True:
        try:
            attempts = 0
            samples = pd.read_csv(filename, sep='\t')
            break
        except IOError as e:
            attempts = attempts + 1
            if attempts > 1:
                print('ERROR: unable to find the sample info file that '
                    'should be in {0}'.format(filename))
                sys.exit(2)
            else:
                save(cfg, samples_init_df())

    return samples


def save(cfg, samples):
    """Save a dataframe, assumed to be sample info, to the sample info file.

    This tries to save the sample info dataframe to a file.  On an error, it
    prints a message and exists.  It does no validate of the dataframe to
    actually check it is a valid sample info dataset.
    """
    filename = os.path.join(cfg.project_dir, SAMPLE_INFO_FILE)
    try:
        samples.to_csv(filename, sep='\t')
    except IOError as e:
        print('ERROR: unable to create sample info file as {0}'.format(filename))
        sys.exit(2)
    return


def add(args):
    cfg = config.validate(args)
    samples = load(cfg)
    # make sure the sample doesn't already exist
    # make sure the fastq files are present
    # compute md5sums of the fastq files
    # update created date
    save(cfg, samples)
    return


def remove(args):
    cfg = config.validate(args)
    samples = load(cfg)
    print(samples)
    save(cfg, samples)
    return
