import argparse
import sys

import barrelseq


def main():
    """Main entry point into command line tool for 

    This is the main function for the barrelseq package.  It parses the
    arguments from ``sys.argv`` with a call to ``parse_args`` and then
    dispatches the appropriate sub-command function for further processing
    and actual heavy lifting.
    """
    try:
        parser = parse_args(sys.argv[1:])
    #except TypeError as e:
    #    print(repr(e))
    #    sys.exit(1)
    #except ValueError as e:
    #    print(repr(e))
    #    sys.exit(1)
    except Exception as e:
        print('{0}: An unhandled exception occured while parsing command '
              'line arguments.'.format(sys.argv[0]))
        print(repr(e))
        raise e
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
    parser = argparse.ArgumentParser(prog=barrelseq.SCRIPT_NAME)
    subparsers = parser.add_subparsers(
            title='{0} commands'.format(barrelseq.SCRIPT_NAME),
            description='list of commands supported',
            help='command help'
            )
    parser.add_argument('--version', action='version', 
            version='{0} {1}'.format(barrelseq.SCRIPT_NAME,
                barrelseq.__version__
                )
            )
    # config sub-command parser
    parser_config = subparsers.add_parser('config',
            help='config file utility help'
            )

    # sample sub-command parser
    parser_sample = subparsers.add_parser('sample',
            help='sample management help'
            )

    # engine sub-command parser
    parser_engine = subparsers.add_parser('engine',
            help='executation engine help'
            )

    # analysis sub-command parser
    parser_analysis = subparsers.add_parser('analysis',
            help='data analysis help'
            )
    subparser_analysis = parser_analysis.add_subparsers(
            title='analysis commands'.format(barrelseq.SCRIPT_NAME),
            description='list of commands supported',
            help='analysis module help'
            )
    # deseq2 sub-command parser
    parser_deseq2 = subparser_analysis.add_parser('deseq2',
            help='DESeq2 analysis help'
            )
    parser_deseq2.add_argument('--config-file', type=argparse.FileType('r'),
            required=True, help='input configuration file'
            )
    parser_deseq2.add_argument('--output', type=str, required=True,
            help='output file prefix')
    parser_deseq2.add_argument('--group-A', type=str,
            help='name of sample condition of first group')
    parser_deseq2.add_argument('--group-B', type=str,
            help='name of sample condition of second group')
    parser_deseq2.add_argument('--regroup-A', type=str,
            help='ad-hoc first group, e.g. EPS=uMax1,uMax2')
    parser_deseq2.add_argument('--regroup-B', type=str,
            help='ad-hoc second group, e.g. EPS=uMax1,uMax2')
    parser_deseq2.add_argument('--shrinkage', 
            choices=['normal', 'apeglm', 'ashr'],
            default='apeglm',
            help='method for shrinkage estimates, see DESeq2 docs'
            )
    parser_deseq2.add_argument('--transformation',
            choices=['vst', 'rlog'],
            default='vst',
            help='variable stabilizing transformation or regularized log'
            )
    parser_deseq2.add_argument('--generate-figures',
            help='should figures be generated'
            )
    # example sub-command parser as template for future analyses
    parser_example = subparser_analysis.add_parser('example', 
            help='example analysis help')
    parser_example.add_argument('--config-file', type=argparse.FileType('r'),
            required=True, help='input configuration file'
            )
    parser_example.add_argument('--integer-parameter', type=int,
            help='example integer parameter for analysis'
            )

    # extract sub-command parser
    parser_extract = subparsers.add_parser('extract',
            help='data extraction help'
            )
    parser_extract.add_argument('--config-file', type=argparse.FileType('r'),
            required=True, help='input configuration file'
            )
    parser_extract.add_argument('--output', type=argparse.FileType('wb'),
            required=True, help='output filename'
            )
    parser_extract.add_argument('--values',
            choices=['TPM', 'raw', 'RPKM'],
            default='TPM',
            help='what values should be exported'
            )
    parser_extract.add_argument('--format',
            choices=['tsv', 'csv', 'xls'],
            default='tsv',
            help='what file format should be used for the output'
            )
    parser_extract.add_argument('--samples',
            nargs='*',
            help='what samples should be included in the output'
            )
    return parser
