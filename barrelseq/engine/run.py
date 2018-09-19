import argparse
import multiprocessing as mp
from subprocess import Popen

from barrelseq import config
from barrelseq import sample


def make_bwa_cmd(sample_row, cfg):
    return "{} mem -M -t 1 {} {} > {}.sam".format(cfg.bwa_path, cfg.reference_fasta_path,
                                                  " ".join(sample_row['fastq_files'].split(",")), sample_row['name'])


def make_view_cmd(name, cfg):
    return "{0} view -bt {1} -o {2}.bam {2}.sam".format(cfg.samtools_path, cfg.reference_fasta_path, name)


def make_sort_cmd(name, cfg):
    return "{0} sort {1}.bam {1}.sorted".format(cfg.samtools_path, name)


def make_index_cmd(name, cfg):
    return "{} index {}.sorted.bam".format(cfg.samtools_path, name)


def make_htseq_cmd(name, cfg):
    return "python {0} -f bam -m intersection-nonempty -s no -t gene -i ID " \
           "{1}.sorted.bam {2} > {1}.summary.dat".format(cfg.htseq_count_path, name, cfg.reference_gff_path)


def run(args):
    cfg = config.validate(args)
    samples = sample.load(cfg)
    # now the system is ready to go with a populated dataframe with sample info
    # and all of the system configuration data on the cfg object
    # the prints below demonstrate what attributes are on cfg and what
    # the schema of the sample_info table is

    # If save_as_scripts is true, don't run anything, but put it all in a bash file
    # Update run date?

    samples['bwa_cmd'] = samples.apply(lambda x: make_bwa_cmd(x, cfg), axis=1)
    samples['view_cmd'] = samples['name'].map(lambda x: make_view_cmd(x, cfg))
    samples['sort_cmd'] = samples['name'].map(lambda x: make_sort_cmd(x, cfg))
    samples['index_cmd'] = samples['name'].map(lambda x: make_index_cmd(x, cfg))
    samples['htseq_cmd'] = samples['name'].map(lambda x: make_htseq_cmd(x, cfg))

    print(args)
    print('\n\n')
    samples.to_csv("tmp.txt", sep="\t")
    return
