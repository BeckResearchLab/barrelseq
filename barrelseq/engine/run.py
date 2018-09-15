import argparse

from barrelseq import config
from barrelseq import sample


def run(args):
    cfg = config.validate(args)
    samples = sample.load(cfg)
    # now the system is ready to go with a populated datafram with sample info
    # and all of the system configuration data on the cfg object
    # the prints below demonstrate what attributes are on cfg and what
    # the schema of the sample_info table is
    print(cfg)
    print('\n\n')
    print(samples)
    return
