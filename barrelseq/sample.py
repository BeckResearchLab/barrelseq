import argparse
import hashlib
import json
import os
import sys

import pandas as pd

from barrelseq import config


SAMPLE_INFO_FILE = 'sample_info.tsv'
SAMPLE_INFO_COLUMNS = ['name', 'group', 'fastq_files', 'description',
            'added', 'last_run', 'last_analyzed', 'md5sums', 'log'
        ]
CHUNK_SIZE = 1024 * 1024 * 16


def samples_init_df():
    """Creates a skeleton data frame with correct columns for sample info.
    """
    return pd.DataFrame(columns=SAMPLE_INFO_COLUMNS)


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
                    'should be in {}'.format(filename))
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
        samples.to_csv(filename, sep='\t', index=False)
    except IOError as e:
        print('ERROR: unable to create sample info file as {}'.format(filename))
        sys.exit(2)
    return


def add(args):
    cfg = config.validate(args)
    samples = validate(cfg)
    errors = False
    # make sure the sample doesn't already exist
    if args.name in samples.name.values:
        print('ERROR: sample {} already exists'.format(args.name))
        errors = True
    # make sure the fastq files are present
    for filename in args.fastq_files:
        if not os.path.exists(filename):
            print('ERROR: fastq file {} does not exist'.format(filename))
            errors = True
    fastq_str = ','.join(args.fastq_files)
    # compute md5sums of the fastq files
    md5sums = []
    if not errors:
        for filename in args.fastq_files:
            print('computing md5sum for file {}'.format(filename))
            m = hashlib.md5()
            with open(filename, 'rb') as f:
                while True:
                    chunk = f.read(CHUNK_SIZE)
                    if chunk == b'':
                        break
                    m.update(chunk)
            md5sums.append(m.hexdigest())
        md5sums_str = ','.join(md5sums)
        # update created date
        added = pd.Timestamp.now()
        sample_log = [{'ts':added.strftime('%c'),
                'msg':f'sample added {args.name}, '
                f'{args.group}, {fastq_str}, {md5sums_str}, {args.description}'}]
        samples.loc[len(samples)] = [ args.name, args.group, fastq_str,
                args.description, added, None, None, 
                md5sums_str, json.dumps(sample_log)
            ]
    if errors:
        raise RuntimeError('at least one failure occured while '
                'adding the sample')
    save(cfg, samples)
    return


def log_append(log_text, msg):
    log = json.loads(log_text)
    log.append({'ts':pd.Timestamp.now().strftime('%c'),
            'msg':msg})
    return json.dumps(log)


def view(args):
    cfg = config.validate(args)
    samples = validate(cfg)
    sample_info = samples.loc[samples['name'] == args.name]
    print(sample_info)
    return


def edit(args):
    cfg = config.validate(args)
    samples = validate(cfg)
    # make changes based on args
    print('WARNING: not implemented')
    save(cfg, samples)
    return


def remove(args):
    cfg = config.validate(args)
    samples = validate(cfg)
    if args.name in samples.name.values:
        samples = samples[samples.name != args.name]
    else:
        raise RuntimeError('sample {} not found in sample info '
                'file'.format(args.name))
    save(cfg, samples)
    return


def validate(cfg):
    samples = load(cfg)
    for sample in samples.name.values:
        workspace_dir = os.path.join(cfg.workspace_dir, sample)
        if not os.path.isdir(workspace_dir):
            try:
                print('WARNING: workspace directory {} does not exist for '
                        'sample {}\n\tit will be '
                        'created'.format(workspace_dir, sample))
                os.mkdir(workspace_dir)
            except OSError as e:
                print('ERROR: unable to create workspace directory {} for '
                        'sample {}'.format(workspace_dir, sample))
                raise(e)
            if not os.path.isdir(workspace_dir):
                print('ERROR: workspace directory path for sample {} '
                    '({}) is not a directory'.format(sample, workspace_dir))
        else:
            pass
    print('sample file validation succeeded')
    return samples
