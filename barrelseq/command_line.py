import argparse
import sys

import barrelseq


def main():
    """Main entry point into command line tool for this package.

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


# Modified from argparse-subparser-monolithic-help-output for recursion
class _HelpAction(argparse._HelpAction):
    """This object is used to provide a custom help action.

    This object is designed to produce a complete help listing at the top
    level of the command.  Purely a convienence, but a nice one.
    """

    def _recurse_parser(self, parser, prefix):
        """Method to recurse an argparse parse tree while printing.

        This recursive function enumerates the argparse parse tree
        while printing the help text for each element in the tree.
        Between elements, a line made of ``-`` delineates pages.

        Recursive

        Side effect: prints to stdout
        """
        # retrieve subparsers from parser
        subparsers_actions = [
                action for action in parser._actions
                    if isinstance(action, argparse._SubParsersAction)
                ]
        # there will probably only be one subparser_action,
        # but better safe than sorry
        for subparsers_action in subparsers_actions:
            # get all subparsers and print help
            for choice, subparser in subparsers_action.choices.items():
                if not prefix:
                    new_prefix = choice
                else:
                    new_prefix = prefix + ' ' + choice
                print('\n{0} {1}'.format(new_prefix,
                                    (78 - len(new_prefix)) * '-'))
                print(subparser.format_help())
                # recurse through subparsers if they exist
                self._recurse_parser(subparser, new_prefix)
    
    def __call__(self, parser, namespace, values, option_string=None):
        """Entry point to help action for argparse trees

        This is the top level help printer for an argparse tree of 
        sub-commands (subparsers).  It prints the top level help,
        recurses the tree while printing each elements help text
        and then exits the parser.
        """
        parser.print_help()
        self._recurse_parser(parser, '')
        parser.exit()


def parser_create():
    """Factory function for creating a parser object.

    This factory function creates an argparse parser object for use by
    ``parse_args`` and the unit testing harness.

    Returns:
        ``argparse.ArgumentParser`` object suitable for parsing arguments.

    """
    # create the top-level parser for sub-commands
    parser = argparse.ArgumentParser(prog=barrelseq.SCRIPT_NAME, add_help=False)
    parser.add_argument('--help', action=_HelpAction,
            help='full help listing'
            )
    parser.add_argument('--version', action='version', 
            version='{0} {1}'.format(barrelseq.SCRIPT_NAME,
                barrelseq.__version__
                )
            )
    subparsers = parser.add_subparsers(
            title='{0} commands'.format(barrelseq.SCRIPT_NAME)
            )

    # config sub-command parser
    parser_config = subparsers.add_parser('config',
            help='config file utility help',
            add_help=False
            )
    parser_config.add_argument('--help', action=_HelpAction,
            help='full help listing'
            )
    subparser_config = parser_config.add_subparsers(
            title='config file utility commands',
            help='config file module help'
            )
    # config create sub-command parser
    parser_config_create = subparser_config.add_parser('create',
            help='help for creating a config file',
            add_help=False
            )
    parser_config_create.add_argument('--help', action=_HelpAction,
            help='full help listing'
            )
    parser_config_create.add_argument('--config-file', type=argparse.FileType('w'),
            required=True, help='output configuration file'
            )
    parser_config_create.add_argument('--project-name', required=True,
            help='human readable name of project'
            )
    parser_config_create.add_argument('--bwa-path',
            type=argparse.FileType('r'),
            help='path to bwa executable',
            default='bwa'
            )
    parser_config_create.add_argument('--samtools-path',
            type=argparse.FileType('r'),
            help='path to samtools executable',
            default='samtools'
            )
    parser_config_create.add_argument('--htseq-count-path',
            type=argparse.FileType('r'),
            help='path to htseq-count executable',
            default='htseq-count'
            )
    parser_config_create.add_argument('--R', type=argparse.FileType('r'),
            help='path to R executable',
            default='R'
            )
    parser_config_create.add_argument('--project-dir',
            type=argparse.FileType('r'),
            required=True, help='path to project directory; '
                'this should be a full path, thus it starts with /'
            )
    parser_config_create.add_argument('--reference-name',
            type=argparse.FileType('r'),
            required=True, help='human readable name of reference species'
            )
    parser_config_create.add_argument('--reference-gff-path', 
            type=argparse.FileType('r'),
            required=True, help='path to reference GFF file'
            )
    parser_config_create.add_argument('--reference-fasta-path',
            type=argparse.FileType('r'),
            required=True, help='path to reference nucleotide FASTA file'
            )
    parser_config_create.add_argument('--pair-ended',
            help='this study\'s data is from pair-ended sequencing runs'
            )
    parser_config_create.add_argument('--opts-bwa-mem', type=str,
            help='passthrough options for bwa in mem mode step'
            )
    parser_config_create.add_argument('--opts-htseq-count', type=str,
            help='passthrough options for htseq-count step'
            )
    parser_config_create.add_argument('--opts-samtools-sam2bam', type=str,
            help='passthrough options for samtools SAM to BAM conversion step'
            )
    parser_config_create.add_argument('--opts-samtools-sort', type=str,
            help='passthrough options for samtools BAM sort step'
            )
    parser_config_create.add_argument('--opts-samtools-index', type=str,
            help='passthrough options for samtools BAM index step'
            )
    # config validate sub-command parser
    parser_config_validate = subparser_config.add_parser('validate',
            help='help for validating a config file',
            add_help=False
            )
    parser_config_validate.add_argument('--help', action=_HelpAction,
            help='full help listing'
            )
    parser_config_validate.add_argument('--config-file', type=argparse.FileType('r'),
            required=True, help='input configuration file for validation'
            )

    # sample sub-command parser
    parser_sample = subparsers.add_parser('sample',
            help='sample management help',
            add_help=False
            )
    parser_sample.add_argument('--help', action=_HelpAction,
            help='full help listing'
            )
    subparser_sample = parser_sample.add_subparsers(
            title='sample management commands',
            help='sample management help'
            )
    # sample add sub-command parser
    parser_sample_add = subparser_sample.add_parser('add',
            help='help for adding a sample',
            add_help=False
            )
    parser_sample_add.add_argument('--help', action=_HelpAction,
            help='full help listing'
            )
    parser_sample_add.add_argument('--name', required=True,
            help='name of sample to be added; must start with a letter '
            'and contain only letters, numbers and underscores(_)'
            )
    parser_sample_add.add_argument('--group', required=True,
            help='name of group this sample belong to; will be added if '
            'it doesn\'t exist; must start with a letter and contain only '
            'letters, numbers and underscores'
            )
    parser_sample_add.add_argument('--fastq-files', required=True,
            help='list of full pathes to fastq file; separated by spaces',
            nargs='+'
            )
    parser_sample_add.add_argument('--description', required=True,
            help='complete description of sample; in single quotes, '
            'not containing any tabs or new lines'
            )
    # sample remove sub-command parser
    parser_sample_remove = subparser_sample.add_parser('remove',
            help='help for removing a sample - DANGER',
            add_help=False
            )
    parser_sample_remove.add_argument('--help', action=_HelpAction,
            help='full help listing'
            )
    parser_sample_remove.add_argument('--name', required=True,
            help='name of sample to be removed - DANGER! - this will '
            'remove a sample from the sample metadata file and may '
            'remove alignment and summary files (workspace files)'
            )
    parser_sample_remove.add_argument('--dont-remove-files',
            help='workspace files for sample will be preserved'
            )

    # engine sub-command parser
    parser_engine = subparsers.add_parser('engine',
            help='execution engine help',
            add_help=False
            )
    parser_engine.add_argument('--help', action=_HelpAction,
            help='full help listing'
            )
    subparser_engine = parser_engine.add_subparsers(
            title='execution engine commands',
            help='execution engine module help'
            )

    # analysis sub-command parser
    parser_analysis = subparsers.add_parser('analysis',
            help='data analysis help',
            add_help=False
            )
    parser_analysis.add_argument('--help', action=_HelpAction,
            help='full help listing'
            )
    subparser_analysis = parser_analysis.add_subparsers(
            title='analysis commands',
            help='analysis module help'
            )
    # deseq2 sub-command parser
    parser_deseq2 = subparser_analysis.add_parser('deseq2',
            help='DESeq2 analysis help',
            add_help=False
            )
    parser_deseq2.add_argument('--help', action=_HelpAction,
            help='full help listing'
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
            help='example analysis help',
            add_help=False
            )
    parser_example.add_argument('--help', action=_HelpAction,
            help='full help listing'
            )
    parser_example.add_argument('--config-file', type=argparse.FileType('r'),
            required=True, help='input configuration file'
            )
    parser_example.add_argument('--integer-parameter', type=int,
            help='example integer parameter for analysis'
            )

    # extract sub-command parser
    parser_extract = subparsers.add_parser('extract',
            help='data extraction help',
            add_help=False
            )
    parser_extract.add_argument('--help', action=_HelpAction,
            help='full help listing'
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
