import argparse
import multiprocessing as mp
import subprocess
import os

from barrelseq import config
from barrelseq import sample


def make_bwa_cmd(sample_row, cfg):
    sam_name = "{0}/{0}.sam".format(sample_row['name'])
    if cfg.opts_bwa_mem is not None:
        optional_args = " " + cfg.opts_bwa_mem + " "
    else:
        optional_args = " "
    return "{} mem -M -t 1{}{} {} -o {}".format(cfg.bwa_path, optional_args, cfg.reference_fasta_path,
                                               " ".join(sample_row['fastq_files'].split(",")), sam_name)


def make_view_cmd(name, cfg):
    sam_name = "{0}/{0}.sam".format(name)
    bam_name = "{0}/{0}.bam".format(name)
    if cfg.opts_samtools_sam2bam is not None:
        optional_args = " " + cfg.opts_samtools_sam2bam + " "
    else:
        optional_args = " "
    return "{} view{}-bt {} -o {} {}".format(cfg.samtools_path, optional_args,
                                             cfg.reference_fasta_path, bam_name, sam_name)


def make_sort_cmd(name, cfg):
    bam_name = "{0}/{0}.bam".format(name)
    sorted_name = "{0}/{0}.sorted.bam".format(name)
    if cfg.opts_samtools_sort is not None:
        optional_args = " " + cfg.opts_samtools_sort + " "
    else:
        optional_args = " "
    return "{} sort {} {} -o {}".format(cfg.samtools_path, optional_args, bam_name, sorted_name)


def make_index_cmd(name, cfg):
    sorted_name = "{0}/{0}.sorted.bam".format(name)
    if cfg.opts_index_mem is not None:
        optional_args = " " + cfg.opts_index_mem + " "
    else:
        optional_args = " "
    return "{} index {} {}".format(cfg.samtools_path, optional_args, sorted_name)


def make_htseq_cmd(name, cfg):
    sorted_name = "{0}/{0}.sorted.bam".format(name)
    summary_name = "{0}/{0}.summary.dat".format(name)
    if cfg.opts_htseq_count is not None:
        optional_args = " " + cfg.opts_htseq_count + " "
    else:
        optional_args = " "
    return "python {}{}-f bam -m intersection-nonempty -s no -t gene -i ID {} {} > " \
           "{}".format(cfg.htseq_count_path, optional_args, sorted_name, cfg.reference_gff_path, summary_name)


def remove_intermediates(name):
    # intermediate_files = ["{0}/{0}.sam".format(name), "{0}/{0}.bam".format(name), "{0}/{0}.sorted.bam".format(name),
    #                     "{0}/{0}.sorted.bam.index".format(name)]
    intermediate_files = ["{0}/{0}.sam".format(name), "{0}/{0}.bam".format(name)]
    for int_file in intermediate_files:
        try:
            os.remove(int_file)
            print("Deleted: {}".format(int_file))
        except FileNotFoundError:
            pass

def run_cmd(cmd):
    err = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
    print(err)
    return cmd


def run(args):
    cfg = config.validate(args)
    samples = sample.validate(cfg)

    os.chdir(os.path.join(cfg.project_dir, "workspace"))

    # now the system is ready to go with a populated dataframe with sample info
    # and all of the system configuration data on the cfg object
    # the prints below demonstrate what attributes are on cfg and what
    # the schema of the sample_info table is

    if args.samples is not None:
        samples = samples.loc[samples['name'].isin(args.samples)]
    if args.samples is not None:
        for smp in args.samples:
            if smp not in samples['name'].values:
                raise KeyError(f"Invalid sample name", f"The sample {smp} is not contained in the list of valid samples.")
    samples = samples.loc[samples['name'].isin(args.samples)]


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
        cmd_list.append(samples[step].tolist())

    if args.save_as_scripts:
        with open("barrelseq.sh", 'w') as b:
            b.write("#!/bin/bash\n\n")
            for step in cmd_list:
                for specific_command in step:
                    b.write(specific_command)
                    b.write("\n")
                b.write("\n")
    else:
        for step in cmd_list:
            if args.processes is None:
                args.processes = 1
            if len(step) < args.processes:
                args.processes = len(step)

            if args.processes > 1:
                with mp.Pool(processes=args.processes) as pool:
                    output = "\n".join(pool.map(run_cmd, step))
            else:
                output = ""
                for specific_command in step:
                    output += run_cmd(specific_command)
            print(output)

    if not args.save_intermediate_files:
        samples['name'].map(remove_intermediates)

    return
