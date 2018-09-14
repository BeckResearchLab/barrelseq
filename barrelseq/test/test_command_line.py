import unittest

import command_line


class TestCommandLine(unittest.TestCase):
    def setUp(self):
        self.parser = command_line.create_parser()

    def test_parser_subcommand(self):
        """place holder from https://stackoverflow.com/questions/18160078/how-do-you-write-tests-for-the-argparse-portion-of-a-python-module"""
        parsed = self.parser.parse_args(['config'])
        self.assertEqual(parsed.subcommand, 'config')


if __name__ == '__main__':
    unittest.main()
