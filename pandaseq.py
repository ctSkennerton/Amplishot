#!/usr/bin/env python

from cogent.app.parameters import FlagParameter, ValuedParameter
from cogent.app.util import CommandLineApplication


"""Application controller for pandaseq"""


__author__ = "Connor Skennerton"
__copyright__ = "Copyright 2012-2013, Connor Skennerton"
__credits__ = ["Connor Skennerton"]
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Connor Skennerton"
__email__ = "c.skennerton@gmail.com"
__status__ = "Development"


class Pandaseq(CommandLineApplication):
    """Simple pandaseq application controller.
    """
    _parameters = {
        '-6': FlagParameter('-', '6'),
        '-B': FlagParameter('-', 'B'),
        '-C': ValuedParameter('-', 'C', Delimiter=' '),
        '-d': ValuedParameter('-', 'd', Delimiter=' '),
        '-f': ValuedParameter('-', 'f', Delimiter=' ', IsPath=True),
        '-F': FlagParameter('-', 'F'),
        '-j': FlagParameter('-', 'j'),
        '-l': ValuedParameter('-', 'l', Delimiter=' '),
        '-L': ValuedParameter('-', 'L', Delimiter=' '),        
        '-N': FlagParameter('-', 'N'),
        '-o': ValuedParameter('-', 'o', Delimiter=' '),
        '-p': ValuedParameter('-', 'p', Delimiter=' '),
        '-q': ValuedParameter('-', 'q', Delimiter=' '),
        '-r': ValuedParameter('-', 'r', Delimiter=' ', IsPath=True),
        '-t': ValuedParameter('-', 't', Delimiter=' '),
        '-T': ValuedParameter('-', 'T', Delimiter=' '),
        }
    _input_handler = '_input_as_parameter'
    _command = 'pandaseq'

    def _accept_exit_status(self, exit_status):
        """Accept an exit status of 0 for the pandaseq program.
        """
        return exit_status == 0
    
    def _input_as_parameter(self, data):
        """ data is a hash that gets passed through to __call__
        """
        self.Parameters['-f'].on(data['forward'])
        self.Parameters['-r'].on(data['reverse'])
        return ''
