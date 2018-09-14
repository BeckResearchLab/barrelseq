# this is not yet a real setup.py but a placeholder that indicates how
# the command line portion of this package should be shimmed in
# for details, see this page:
# https://python-packaging.readthedocs.io/en/latest/command-line-scripts.html

setup(
    ...
    entry_points = {
        'console_scripts': ['barrelseq=barrelseq.command_line:main'],
    }
    ...
)
