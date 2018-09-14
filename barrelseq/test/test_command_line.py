import os
import shutil
import tempfile
import unittest

from barrelseq import command_line


PROJECT_NAME = 'test project'
REFERENCE_NAME = 'test organism'
GFF_PATH = os.path.join(os.path.dirname(__file__), 'test.gff')
FASTA_PATH = os.path.join(os.path.dirname(__file__), 'test.fasta')
BWA_PATH = os.path.join(os.path.dirname(__file__), 'bwa')
BWA_MEM_OPTS = os.path.join(os.path.dirname(__file__), '-v')

class TestCommandLine(unittest.TestCase):
    def setUp(self):
        self.parser = command_line.parser_create()
        fi, self.config_file = tempfile.mkstemp()
        os.fdopen(fi).close()
        self.project_dir = tempfile.mkdtemp()
        args = self.parser.parse_args([
                'config', 'create',
                '--config-file', self.config_file,
                '--project-dir', self.project_dir,
                '--project-name', PROJECT_NAME,
                '--reference-name', REFERENCE_NAME,
                '--reference-gff-path', GFF_PATH,
                '--reference-fasta-path', FASTA_PATH,
                '--bwa-path', BWA_PATH,
                '--opts-bwa-mem', BWA_MEM_OPTS
            ])
        args.config_file.close()

    def tearDown(self):
        shutil.rmtree(self.project_dir)
        os.remove(self.config_file)

    def test_parser_subcommand(self):
        args = self.parser.parse_args(['config', 'edit', '--config-file', '/tmp/b', '--project-name', 'TestZ'])
        self.assertEqual(args.command, 'config')
        args.config_file.close()


if __name__ == '__main__':
    unittest.main()
