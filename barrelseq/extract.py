import argparse
import os

import pandas as pd

from barrelseq import config
from barrelseq import sample


def extract(args):
    '''
    After running htseq-count on a set of specified RNA-seq samples, 
    combine all of their summary count files into one master dataframe
    of gene rows and experiment sample columns.
    '''

    # validate config and retrieve specified samples (use all by default)
    cfg = config.validate(args)
    samples = sample.validate(cfg)

    os.chdir(os.path.join(cfg.project_dir, "workspace"))

    # If sample names were provided from the command line, check if
    # it exists in the config
    if args.samples is not None:
        for smp in args.samples:
            # raise error if sample name provided was not in the config
            if smp not in samples['name'].values:
                raise KeyError(f"Invalid sample name", f"The sample {smp} is not contained in the list of valid samples.")
        
        # filter to only the samples passed in the args
        samples = samples.loc[samples['name'].isin(args.samples)]

    # initialize empty df to collect experiment sample columns
    master_df = None
    # initialize a row counter for catching changes in upstream software
    rows = 0

    # for every sample
    for smp in samples['name'].values:
        # load up its summary count file
        summary_file = os.path.join(smp, smp + '.summary.dat')
        
        # If this is the first sample being added, initiialize df
        if master_df is None:
            master_df = pd.read_csv(summary_file, sep='\t', names=['locus_tag', smp])
            rows = master_df.shape[0]
        
        # otherwise, merge this new sample in with the master df
        else:
            # make a temp df with the new experiment file
            tmp_df = pd.read_csv(summary_file, sep='\t', names=['locus_tag', smp])
            
            # check if this experiment has a different number of rows in the summary
            # (if so, something probably changed in upstream counting software)
            if rows != tmp_df.shape[0]:
                raise ValueError("Data frame size mismatch", f"The number of rows in the master data frame is {rows} and the number of rows in the data frame for smp {smp} contains {tmp_df.shape[0]}")
            
            # merge in the new experiment, joining on the gene names
            master_df = master_df.merge(tmp_df, how='outer', on='locus_tag')
            
            # check again for a row count mismatch, just to be sure
            if rows != master_df.shape[0]:
                raise(ValueError, "Data frame size mismatch", f"The number of rows in the old master data frame is {rows} and the number of rows in the merged data frame including smp {smp} contains {master_df.shape[0]}")

    # save final merged df a tsv
    print(f'saving output to {args.output}')
    master_df.to_csv(args.output, sep='\t', header=True, index=False)

    return
