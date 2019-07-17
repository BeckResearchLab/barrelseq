import argparse
import os

import pandas as pd

from barrelseq import config
from barrelseq import sample


def extract(args):
    cfg = config.validate(args)
    samples = sample.validate(cfg)

    os.chdir(os.path.join(cfg.project_dir, "workspace"))

    if args.samples is not None:
        for smp in args.samples:
            if smp not in samples['name'].values:
                raise KeyError(f"Invalid sample name", f"The sample {smp} is not contained in the list of valid samples.")
        samples = samples.loc[samples['name'].isin(args.samples)]

    master_df = None
    rows = 0
    for smp in samples['name'].values:
        summary_file = os.path.join(smp, smp + '.summary.dat')
        if master_df is None:
            master_df = pd.read_csv(summary_file, sep='\t', names=['locus_tag', smp])
            rows = master_df.shape[0]
        else:
            tmp_df = pd.read_csv(summary_file, sep='\t', names=['locus_tag', smp])
            if rows != tmp_df.shape[0]:
                raise ValueError("Data frame size mismatch", f"The number of rows in the master data frame is {rows} and the number of rows in the data frame for smp {smp} contains {tmp_df.shape[0]}")
            master_df = master_df.merge(tmp_df, how='outer', on='locus_tag')
            if rows != master_df.shape[0]:
                raise(ValueError, "Data frame size mismatch", f"The number of rows in the old master data frame is {rows} and the number of rows in the merged data frame including smp {smp} contains {master_df.shape[0]}")


    print(f'saving output to {args.output}')
    master_df.to_csv(args.output, sep='\t', header=True, index=False)

    return
