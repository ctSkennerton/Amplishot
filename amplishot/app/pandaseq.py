#!/usr/bin/env python
###############################################################################
#
# pandaseq.py - application controller for pandaseq
#
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__author__ = "Connor Skennerton"
__copyright__ = "Copyright 2013"
__credits__ = ["Connor Skennerton"]
__license__ = "GPL3"
__version__ = "0.3.1"
__maintainer__ = "Connor Skennerton"
__email__ = "c.skennerton@gmail.com"
__status__ = "Development"

###############################################################################
from cogent.app.parameters import FlagParameter, ValuedParameter
from amplishot.app.util import ExtendedCommandLineApplication
#from cogent.app.util import CommandLineApplication

class Pandaseq(ExtendedCommandLineApplication):
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
    _command = 'pandaseq'

    def _accept_exit_status(self, exit_status):
        """Accept an exit status of 0 for the pandaseq program.
        """
        return exit_status == 0
