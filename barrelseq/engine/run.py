import argparse
import multiprocessing as mp
import subprocess

from barrelseq import config
from barrelseq import sample


def make_bwa_cmd(sample_row, cfg):
    if cfg.opts_bwa_mem is not "null":
        optional_args = cfg.opts_bwa_mem
    else:
        optional_args = ""
    return "{} mem -M -t 1 {} {} {} > {}.sam".format(cfg.bwa_path, optional_args, cfg.reference_fasta_path,
                                                  " ".join(sample_row['fastq_files'].split(",")), sample_row['name'])


def make_view_cmd(name, cfg):
    if cfg.opts_samtools_sam2bam is not "null":
        optional_args = cfg.opts_bwa_mem
    else:
        optional_args = ""
    return "{0} view {1} -bt {2} -o {3}.bam {3}.sam".format(cfg.samtools_path,
                                                            optional_args, cfg.reference_fasta_path, name)


def make_sort_cmd(name, cfg):
    if cfg.opts_samtools_sort is not "null":
        optional_args = cfg.opts_bwa_mem
    else:
        optional_args = ""
    return "{0} sort {1} {2}.bam {2}.sorted".format(cfg.samtools_path, optional_args, name)


def make_index_cmd(name, cfg):
    if cfg.opts_index_mem is not "null":
        optional_args = cfg.opts_bwa_mem
    else:
        optional_args = ""
    return "{} index {} {}.sorted.bam".format(cfg.samtools_path, optional_args, name)


def make_htseq_cmd(name, cfg):
    if cfg.opts_htseq_count is not "null":
        optional_args = cfg.opts_bwa_mem
    else:
        optional_args = ""
    return "python {0} {1} -f bam -m intersection-nonempty -s no -t gene -i ID {2}.sorted.bam {3} > " \
           "{2}.summary.dat".format(cfg.htseq_count_path, optional_args, name, cfg.reference_gff_path)


def run_cmd(cmd):
    err = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
    print(err)
    return err


def run(args):
    cfg = config.validate(args)
    samples = sample.load(cfg)
    # now the system is ready to go with a populated dataframe with sample info
    # and all of the system configuration data on the cfg object
    # the prints below demonstrate what attributes are on cfg and what
    # the schema of the sample_info table is

    # If save_as_scripts is true, don't run anything, but put it all in a bash file
    # Update run date?
    faidx = cfg.reference_fasta_path + ".fai"
    try:
        f = open(faidx, 'r')
        f.close()
    except FileNotFoundError:
        faidx_cmd = "{} faidx {}".format(cfg.samtools_path, cfg.reference_fasta_path)
        run_cmd(faidx_cmd)


    samples['bwa_cmd'] = samples.apply(lambda x: make_bwa_cmd(x, cfg), axis=1)
    samples['view_cmd'] = samples['name'].map(lambda x: make_view_cmd(x, cfg))
    samples['sort_cmd'] = samples['name'].map(lambda x: make_sort_cmd(x, cfg))
    samples['index_cmd'] = samples['name'].map(lambda x: make_index_cmd(x, cfg))
    samples['htseq_cmd'] = samples['name'].map(lambda x: make_htseq_cmd(x, cfg))

    cmd_list = []
    for step in ['bwa_cmd', 'view_cmd', 'sort_cmd', 'index_cmd', 'htseq_cmd']:
        cmd_list.extend(samples[step].tolist())
    print(cmd_list)

    if args.save_as_scripts:
        pass

    print(args)
    print('\n\n')
    samples.to_csv("tmp.txt", sep="\t")
    return
