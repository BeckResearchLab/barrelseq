import argparse
import os
import tempfile
import types
import unittest

from barrelseq import config


class TestConfig(unittest.TestCase):
    def setUp(self):
        fi, self.config_file = tempfile.mkstemp()
        os.fdopen(fi).close()

    def tearDown(self):
        os.remove(self.config_file)

    def _config_create(self):
        args = argparse.Namespace(config_file=open(self.config_file, 'r+'),
                    test_option='test value',
                    test_option2='test value 2'
                )
        config.create(args)
        args.config_file.close()

    def test_config_create(self):
        self._config_create()
        with open(self.config_file, 'r') as f:
            lines = f.readlines()
        self.assertEqual(len(lines), 2)
        self.assertEqual(lines[0], 'test_option: test value\n')

    def test_config_load(self):
        self._config_create()
        args = argparse.Namespace(config_file=open(self.config_file, 'r'))
        config_data = config.load(args)
        args.config_file.close()
        self.assertEqual(type(config_data), types.SimpleNamespace)
        self.assertEqual(config_data.test_option, 'test value')
        self.assertEqual(config_data.test_option2, 'test value 2')

    def test_config_edit(self):
        self._config_create()
        args = argparse.Namespace(config_file=open(self.config_file, 'r+'),
                test_option='new test value',
                test_option2=None,
                test_option3='new test option')
        config_data = config.edit(args)
        args.config_file.close()
        with open(self.config_file, 'r') as f:
            lines = f.readlines()
        self.assertEqual(len(lines), 3)
        self.assertEqual(lines[0], 'test_option: new test value\n')
        self.assertEqual(lines[1], 'test_option2: test value 2\n')
        self.assertEqual(lines[2], 'test_option3: new test option\n')
        

if __name__ == '__main__':
    unittest.main()
