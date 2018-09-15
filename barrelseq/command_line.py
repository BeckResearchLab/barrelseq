import argparse
import os
import sys

import barrelseq
from barrelseq.version import __version__
from barrelseq import config
from barrelseq import sample
from barrelseq import engine
from barrelseq import analysis
from barrelseq import extract


def _value_is_positive(value):
    """Type checker for argparse parser that requires an int >=1.
    """
    ivalue = int(value)
    if ivalue < 1:
        raise argparse.ArgumentTypeError(
                '{0} must be >= 1'.format(ivalue))


PARSE_TREE = [
    ['--version', {
        'action': 'version',
        'version': '{} {}'.format(__name__.split('.')[0],
            __version__)
    }],
    {
        'title': '{} commands'.format(__name__.split('.')[0]),
        'dest': 'command',
        'required': True,
        'help': 'top level commands',
        'subcommands': [{
            'command': 'sample',
            'help': 'sample management',
            'sub_tree': [{
                'dest': 'sample_command',
                'title': 'sample management',
                'help': 'sample management commands',
                'required': True,
                'subcommands': [{
                    'command': 'add',
                    'help': 'adding samples',
                    'title': 'sample add',
                    'func': sample.add,
                    'sub_tree': [
                        ['--config-file', {
                            'type': argparse.FileType('r'),
                            'required': True,
                            'help': 'input configuration file'
                        }],
                        ['--name', {
                            'type': str,
                            'required': True,
                            'help': 'name of sample to be added; '
                            'must start with a letter and contain '
                            'only letters, numbers and underscores(_)'
                        }],
                        ['--group', {
                            'type': str,
                            'required': True,
                            'help': 'name of group this sample belong '
                            'to; will be added if it doesn\'t exist; '
                            'must start with a letter and contain '
                            'only letters, numbers and underscores'
                        }],
                        ['--fastq-files', {
                            'required': True,
                            'nargs': '+',
                            'help': 'list of full pathes to fastq files'
                        }],
                        ['--description', {
                            'required': True,
                            'help': 'complete description of sample; in'
                            'single quotes not containing any tabs or'
                            'new lines'
                        }]
                    ]
                }, {
                    'command': 'remove',
                    'help': 'removing samples - DANGER!',
                    'title': 'sample remove',
                    'func': sample.remove,
                    'sub_tree': [
                        ['--config-file', {
                            'type': argparse.FileType('r'),
                            'required': True,
                            'help': 'input configuration file'
                        }],
                        ['--name', {
                            'type': str,
                            'required': True,
                            'help': 'name of sample to be removed - '
                            'DANGER! - this will remove a sample '
                            'from the sample metadata file and '
                            'may remove alignment and summary '
                            'files (workspace files)'
                        }],
                        ['--dont-remove-files', {
                            'action': 'store_true',
                            'help': 'workspace files for sample will be preserved'
                        }]
                    ]
                }]
            }]
        }, {
            'command': 'engine',
            'help': 'execution engine control',
            'sub_tree': [{
                'dest': 'engine_command',
                'title': 'engine management',
                'help': 'execution engine commands',
                'required': True,
                'subcommands': [{
                    'command': 'run',
                    'help': 'execution engine run mode',
                    'title': 'engine run',
                    'func': engine.run,
                    'sub_tree': [
                        ['--config-file', {
                            'type': argparse.FileType('r'),
                            'required': True,
                            'help': 'input configuration file'
                        }],
                        ['--processes', {
                            'type': _value_is_positive,
                            'help': 'number of parallel processes to run'
                        }],
                        ['--samples', {
                            'nargs': '*',
                            'help': 'list of samples to be included in '
                            'analysis, otherwise all eligible will '
                            'be included'
                        }],
                        ['--save-as-scripts', {
                            'action': 'store_true',
                            'help': 'instead of running the analysis, '
                            'generate shell scripts to run the '
                            'commands; useful for manual inspection '
                            'of the commands or for submission to a '
                            'queueing system'
                        }],
                        ['--save-intermediate-files', {
                            'action': 'store_true',
                            'help': 'intermediate processing files such as '
                            'uncompressed SAM files will be saved for '
                            'later inspection'
                        }]
                    ]
                }]
            }]
        }, {
            'command': 'analysis',
            'help': 'stats and other analyses',
            'sub_tree': [

                {
                    'dest': 'analysis_command',
                    'title': 'analysis modules',
                    'help': 'analysis commands',
                    'required': True,
                    'subcommands': [{
                        'command': 'deseq2',
                        'help': 'DESeq2 fold change analyses in R',
                        'title': 'DESeq2 based analyses',
                        'func': analysis.deseq2,
                        'sub_tree': [
                            ['--config-file', {
                                'type': argparse.FileType('r'),
                                'required': True,
                                'help': 'input configuration file'
                            }],
                            ['--output', {
                                'type': str,
                                'required': True,
                                'help': 'output file prefix, defaults to'
                                '{groupA name}_vs_{groupB name}'
                            }],
                            ['--group-A', {
                                'type': str,
                                'help': 'name of sample condition of first group'
                            }],
                            ['--group-B', {
                                'type': str,
                                'help': 'name of sample condition of second group'
                            }],
                            ['--regroup-A', {
                                'type': str,
                                'help': 'ad-hoc first group, e.g. EPS=uMax1,uMax2'
                            }],
                            ['--regroup-B', {
                                'type': str,
                                'help': 'ad-hoc second group, e.g. EPS=uMax1,uMax2'
                            }],
                            ['--shrinkage', {
                                'type': str,
                                'choices': ['normal', 'apeglm', 'ashr'],
                                'default': 'apeglm',
                                'help': 'method for shrinkage estimates, see DESeq2 docs'
                            }],
                            ['--transformation', {
                                'type': str,
                                'choices': ['vst', 'rlog'],
                                'default': 'vst',
                                'help': 'variable stabilizing transformation or regularized log'
                            }],
                            ['--generate-figures', {
                                'action': 'store_true',
                                'default': False,
                                'help': 'should figures be generated, default is False'
                            }]

                        ]
                    }, {
                        'command': 'example',
                        'help': 'Example template of an analysis module',
                        'func': analysis.example,
                        'title': 'example template',
                        'sub_tree': [
                            ['--config-file', {
                                'type': argparse.FileType('r'),
                                'required': True,
                                'help': 'input configuration file'
                            }],
                            ['--integer-parameter', {
                                'type': int,
                                'help': 'example typed argument to analysis'
                            }]
                        ]
                    }]
                }

            ]
        }, {
            'command': 'extract',
            'help': 'data extraction utility',
            'func': extract.extract,
            'sub_tree': [
                ['--config-file', {
                    'type': argparse.FileType('r'),
                    'required': True,
                    'help': 'input configuration file'
                }],
                ['--output', {
                    'type': argparse.FileType('wb'),
                    'required': True,
                    'help': 'output filename prefix'
                }],
                ['--values', {
                    'choices': ['TPM', 'raw', 'RPKM'],
                    'default': 'TPM',
                    'help': 'what values should be exported'
                }],
                ['--format', {
                    'choices': ['tsv', 'csv', 'xls'],
                    'default': 'tsv',
                    'help': 'what file format should be used for the output'
                }],
                ['--samples', {
                    'nargs': '*',
                    'help': 'what samples should be included in the output'
                }]

            ]
        }, {
            'command': 'config',
            'help': 'config file utility',
            'sub_tree': [{
                'title': 'config file utility',
                'dest': 'config_command',
                'required': True,
                'help': 'config management commands',
                'subcommands': [{
                        'command': 'create',
                        'title': 'config file creation utility',
                        'help': 'creating a config file from scratch',
                        'func': config.create,
                        'sub_tree': [
                            ['--config-file', {
                                'type': argparse.FileType('w+'),
                                'required': True,
                                'help': 'output configuration file'
                            }],
                            ['--project-dir', {
                                'type': str,
                                'required': True,
                                'help': 'path to project directory; '
                                'this should be a full path'
                            }],
                            ['--project-name', {
                                'type': str,
                                'required': True,
                                'help': 'human readable name of project'
                            }],
                            ['--bwa-path', {
                                'type': str,
                                'default': 'bwa',
                                'help': 'path to bwa executable'
                            }],
                            ['--samtools-path', {
                                'type': str,
                                'default': 'samtools',
                                'help': 'path to samtools executable'
                            }],
                            ['--htseq-count-path', {
                                'type': str,
                                'default': 'htseq-count',
                                'help': 'path to htseq-count executable'
                            }],
                            ['--R-path', {
                                'type': str,
                                'default': 'R',
                                'help': 'path to R executable'
                            }],
                            ['--reference-name', {
                                'type': str,
                                'required': True,
                                'help': 'human readable name of '
                                'reference species'
                            }],
                            ['--reference-gff-path', {
                                'type': str,
                                'required': True,
                                'help': 'path to reference GFF file'
                            }],
                            ['--reference-fasta-path', {
                                'type': str,
                                'required': True,
                                'help': 'path to reference nucleotide '
                                'FASTA file'
                            }],
                            ['--pair-ended', {
                                'action': 'store_true',
                                'default': False,
                                'help': 'this study\'s data '
                                'is from pair-ended sequencing runs'
                            }],
                            ['--opts-bwa-mem', {
                                'type': str,
                                'help': 'passthrough options for bwa in '
                                'mem mode step'
                            }],
                            ['--opts-samtools-sam2bam', {
                                'type': str,
                                'help': 'passthrough options for samtools in '
                                'SAM to BAM conversion step'
                            }],
                            ['--opts-samtools-sort', {
                                'type': str,
                                'help': 'passthrough options for samtools in '
                                'BAM sort step'
                            }],
                            ['--opts-index-mem', {
                                'type': str,
                                'help': 'passthrough options for samtools in '
                                'BAM index step'
                            }],
                            ['--opts-htseq-count', {
                                'type': str,
                                'help': 'passthrough options for '
                                'htseq-count step'
                            }]
                        ]
                    },
                    {
                        'command': 'edit',
                        'help': 'help for editing a config file',
                        'func': config.edit,
                        'sub_tree': [
                            ['--config-file', {
                                'type': argparse.FileType('r+'),
                                'required': True,
                                'help': 'output configuration file'
                            }],
                            ['--project-dir', {
                                'type': str,
                                'required': False,
                                'help': 'path to project directory; '
                                'this should be a full path'
                            }],
                            ['--project-name', {
                                'type': str,
                                'required': False,
                                'help': 'human readable name of project'
                            }],
                            ['--bwa-path', {
                                'type': str,
                                'default': None,
                                'required': False,
                                'help': 'path to bwa executable'
                            }],
                            ['--samtools-path', {
                                'type': str,
                                'default': None,
                                'required': False,
                                'help': 'path to samtools executable'
                            }],
                            ['--htseq-count-path', {
                                'type': str,
                                'default': None,
                                'required': False,
                                'help': 'path to htseq-count executable'
                            }],
                            ['--R-path', {
                                'type': str,
                                'default': 'R',
                                'required': False,
                                'help': 'path to R executable'
                            }],
                            ['--reference-name', {
                                'type': str,
                                'required': False,
                                'help': 'human readable name of '
                                'reference species'
                            }],
                            ['--reference-gff-path', {
                                'type': str,
                                'required': False,
                                'help': 'path to reference GFF file'
                            }],
                            ['--reference-fasta-path', {
                                'type': str,
                                'required': False,
                                'help': 'path to reference nucleotide '
                                'FASTA file'
                            }],
                            ['--pair-ended', {
                                'action': 'store_true',
                                'default': None,
                                'help': 'this study\'s data '
                                'is from pair-ended sequencing runs'
                            }],
                            ['--opts-bwa-mem', {
                                'type': str,
                                'help': 'passthrough options for bwa in '
                                'mem mode step'
                            }],
                            ['--opts-samtools-sam2bam', {
                                'type': str,
                                'help': 'passthrough options for samtools in '
                                'SAM to BAM conversion step'
                            }],
                            ['--opts-samtools-sort', {
                                'type': str,
                                'help': 'passthrough options for samtools in '
                                'BAM sort step'
                            }],
                            ['--opts-index-mem', {
                                'type': str,
                                'help': 'passthrough options for samtools in '
                                'BAM index step'
                            }],
                            ['--opts-htseq-count', {
                                'type': str,
                                'help': 'passthrough options for '
                                'htseq-count step'
                            }]
                        ]
                    }
                ]
            }]
        }]
    }
]

def main():
    """Main entry point into command line tool for this package.

    This is the main function for the barrelseq package.  It parses the
    arguments from ``sys.argv`` with a call to ``parse_args`` and then
    dispatches the appropriate sub-command function for further processing
    and actual heavy lifting.
    """
    args = parse_args(sys.argv[1:])
    args.func(args)
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


def parser_create_from_tree(parse_tree, root):
    """Recursive function for building parse tree from data structure.

    The PARSE_TREE object at the top of this file contains the layout
    of the command line user interface (i.e. its parse tree).  That object
    is used by this function to build the parse tree with ``argparse``
    by recursing that data structure's tree.  It is easier to manage the
    user interface by modifying the data structure than several hundred
    function calls in the same pattern.
    """
    for arg in parse_tree:
        if isinstance(arg, list):
            #print(f'{arg[0]}')
            root.add_argument(arg[0], **arg[1])
        elif isinstance(arg, dict):
            subparser = root.add_subparsers(title=arg['title'],
                    help=arg['help'], dest=arg['dest'])
            if 'required' in arg.keys():
                if arg['required']:
                    subparser.required = True
            for subcmd in arg['subcommands']:
                #print(f'{subcmd["command"]}')
                child_parser = subparser.add_parser(subcmd['command'],
                            help=subcmd['help'])
                if 'func' in subcmd.keys():
                    child_parser.set_defaults(func=subcmd['func'])
                if 'sub_tree' in subcmd.keys():
                    parser_create_from_tree(subcmd['sub_tree'], child_parser)
                
        else:
            raise TypeError(f'unexpected object of type {type(toplvl)} '
                'found in parse tree')


def parser_create():
    """Factory function for creating a parser object.

    This factory function creates an argparse parser object for use by
    ``parse_args`` and the unit testing harness.

    Returns:
        ``argparse.ArgumentParser`` object suitable for parsing arguments.

    """
    # create the top-level parser for sub-commands
    parser = argparse.ArgumentParser(prog=barrelseq.SCRIPT_NAME)
    parser_create_from_tree(PARSE_TREE, parser)
    return parser
