import argparse
import os

from barrelseq import config
from barrelseq import sample


def extract(args):
    cfg = config.validate(args)
    samples = sample.validate(cfg)

    os.chdir(os.path.join(cfg.project_dir, "workspace"))

    if args.samples is not None:
        samples = samples.loc[samples['name'].isin(args.samples)]

    print(f'saving output to {args.output}')

    return
