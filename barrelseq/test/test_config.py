import argparse
import os
import tempfile
import unittest

from barrelseq import config


class TestConfig(unittest.TestCase):
    def setUp(self):
        fi, self.config_file = tempfile.mkstemp()
        os.fdopen(fi).close()

    def tearDown(self):
        os.remove(self.config_file)

    def test_config_create(self):
        args = argparse.Namespace(config_file=open(self.config_file, 'r+'),
                    test_option='test value'
                )
        config.create(args)
        args.config_file.close()
        with open(self.config_file, 'r') as f:
            lines = f.readlines()
        self.assertEqual(len(lines), 1)
        self.assertEqual(lines[0], 'test_option: test value\n')


if __name__ == '__main__':
    unittest.main()
