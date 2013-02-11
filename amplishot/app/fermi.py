#!/usr/bin/env python
###############################################################################
#
# simple application controller for fermi
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
__version__ = "0.1.2"
__maintainer__ = "Connor Skennerton"
__email__ = "c.skennerton@gmail.com"
__status__ = "Development"

###############################################################################
import os
from amplishot.app.util import ExtendedCommandLineApplication
from cogent.app.util import ValuedParameter, FlagParameter, ResultPath

class Fermi(ExtendedCommandLineApplication):
    _parameters = {
            '-c': FlagParameter('-', 'c'),
            '-P': FlagParameter('-', 'P'),
            '-B': FlagParameter('-', 'B'),
            '-e': ValuedParameter('-', 'e', Delimeter=' ', IsPath=True),
            '-p': ValuedParameter('-', 'p', Delimeter=' ', Value='fmdef'),
            '-t': ValuedParameter('-', 't', Delimeter=' '),
            '-l': ValuedParameter('-', 'l', Delimeter=' '),
            '-k': ValuedParameter('-', 'k', Delimeter=' ')
            }
    _command = 'run-fermi.pl'

    _input_handler = '_input_as_paths'

    def _get_result_paths(self):
        results = {}
        prefix = self.Parameters['-p'].Value
        suffixes = ['raw.fmd', 'ec.fq.gz', 'ec.fmd', 'p0.mag.gz', 'p1.mag.gz',
                'p2.mag.gz']
        if self.Parameters['-c'].isOn() or self.Parameters['-P'].isOn():
            suffixes.extend(['p4.fa.gz', 'p5.fq.gz'])
            
        for suffix in suffixes:
            results[suffix] = ResultPath(Path=os.path.join(self.WorkingDir,
                preifx + '.' + suffix), IsWritten=True)
        
        return results
