import argparse
import multiprocessing as mp
from subprocess import Popen

from barrelseq import config
from barrelseq import sample


def make_bwa_cmd(sample_row, cfg):
    return "{} mem -M -t 1 {} {} > {}.sam".format(cfg.bwa_path, cfg.reference_fasta_path,
                                                              sample_row['fastq_files'], sample_row.name)


def make_view_cmd(name, cfg):
    return "{0} view -bt {1} -o {2}.bam {2}.sam".format(cfg.samtools_path, cfg.reference_fasta_path, name)


def make_view_cmd(name, cfg):
    return "{0} sort {1}.bam {1}.sorted".format(cfg.samtools_path, name)


def run(args):
    cfg = config.validate(args)
    samples = sample.load(cfg)
    # now the system is ready to go with a populated dataframe with sample info
    # and all of the system configuration data on the cfg object
    # the prints below demonstrate what attributes are on cfg and what
    # the schema of the sample_info table is

    # If save_as_scripts is true, don't run anything, but put it all in a bash file
    # Update run date?
    bwa_cmd = {}
    view_command = {}
    sort_command = {}
    index_command = {}
    htseq_command = {}

    samples['bwa_cmd'] = samples.apply(lambda x: make_bwa_cmd(x, cfg))
    samples['view_cmd'] = samples['name'].map(lambda x: make_view_cmd(x, cfg))
    samples['sort_cmd'] = samples['name'].map(lambda x: make_sort_cmd(x, cfg))



    print(args)
    print('\n\n')
    samples.to_csv("tmp.txt", sep="\t")
    return
