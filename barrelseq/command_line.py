import argparse
import sys

import barrelseq


def main():
    """Main entry point into command line tool for barrelseq.

    This is the main function for the barrelseq package.  It parses the
    arguments from ``sys.argv`` with a call to ``parse_args`` and then
    dispatches the appropriate sub-command function for further processing
    and actual heavy lifting.
    """
    try:
        parser = parse_args(sys.argv[1:])
    except TypeError, e:
        print(repr(e))
        sys.exit(1)
    except ValueError, e:
        print(repr(e))
        sys.exit(1)
    except Exception, e:
        print('{0}: An unhandled exception occured while parsing command '
              'line arguments.'.format(sys.argv[0]))
        print(repr(e))
        sys.exit(1)
    # do something interesting here
    sys.exit(0)


def parse_args(args):
    """Command line argument parser function.

    This function takes a list of command line arguments, usually from
    ``sys.arg[1:]``, but also from constructed lists such as for unit testing
    and does parsing with ``argparse` and handles type checking.  Sanity
    checking of argument values that have side-effects, such as checking
    if a file exists and is readable, is handled downstream.

    Args:
        *args: Variable length argument list.

    Returns:
        ``argparse.Namespace`` object with parsed arguments.

    Raises:
        TypeError: If argument types fail type checking.
        ValueError: If argument values fail pre-checks.

    """
    parser = parser_create()
    return parser.parse_args(args)

def parser_create():
    """Factory function for creating a parser object.

    This factory function creates an argparse parser object for use by
    ``parse_args`` and the unit testing harness.

    Returns:
        ``argparse.ArgumentParser`` object suitable for parsing arguments.

    """
    # create the top-level parser for sub-commands
    parser = argparse.ArgumentParser(prog=SCRIPT_NAME)
    subparsers = parser.add_subparsers(help='command help')
    # config sub-command parser
    parser_config = subparsers.add_parser('config', help='config file utility help')
    # sample sub-command parser
    parser_sample = subparsers.add_parser('sample', help='sample management help')
    # engine sub-command parser
    parser_engine = subparsers.add_parser('engine', help='executation engine help')
    # deseq2 sub-command parser
    parser_deseq2 = subparsers.add_parser('deseq2', help='DESeq2 analysis help')
    # extract sub-command parser
    parser_extract = subparsers.add_parser('extract', help='data extraction help')
    return parser
