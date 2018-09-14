import argparse

import yaml


def create(args):
    print(yaml.dump(args))
    return


def validate(args):
    print(args)
    return
