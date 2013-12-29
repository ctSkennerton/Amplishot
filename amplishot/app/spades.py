#!/usr/bin/env python
###############################################################################
#
# spades.py - application controller for spades 2.5.1
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
__version__ = "0.4.0"
__maintainer__ = "Connor Skennerton"
__email__ = "c.skennerton@gmail.com"
__status__ = "Development"

###############################################################################
import os
from cogent.app.parameters import FlagParameter, ValuedParameter
from cogent.app.util import CommandLineApplication, ResultPath, ApplicationError

class Spades(CommandLineApplication):
    """Simple spades application controller.
    """
    _parameters = {
        '-o': ValuedParameter('-', 'o', Delimiter=' ', IsPath=True),
        '--sc': FlagParameter('--', 'sc'),
        '--12': ValuedParameter('--', '12', Delimiter=' ', IsPath=True),
        '-s': ValuedParameter('-', 's', Delimiter=' ', IsPath=True),
        '-1': ValuedParameter('-', '1', Delimiter=' ', IsPath=True),
        '-2': ValuedParameter('-', '2', Delimiter=' ', IsPath=True),
        '--rectangles': FlagParameter('--', 'rectangles'),
        '--careful': FlagParameter('--', 'careful')
        }
    _synonyms = {
            'rectangles': '--rectangles',
            'careful': '--careful'
            }
    _command = 'spades.py'

    def _accept_exit_status(self, exit_status):
        """Accept an exit status of 0 for the pandaseq program.
        """
        return exit_status == 0
    
    def _get_result_paths(self, data):
        ret = {}
        outdir = self.Parameters['-o'].Value
        ret['contigs'] = ResultPath(Path=os.path.join(self.WorkingDir,outdir,'contigs.fasta'),
                    IsWritten=True)
        if self.Parameters['--12'].isOn() or self.Parameters['-1'].isOn():
            ret['scaffolds'] = ResultPath(Path=os.path.join(self.WorkingDir,outdir,'scaffolds.fasta'),
                        IsWritten=True)
        return ret

    def _handle_app_result_build_failure(self,out,err,exit_status,result_paths):
        return_message = "Problem running spades. exit_status = %s\n%s\n%s\n%s\n" % (str(exit_status),str(out), str(err),
            str(result_paths) )
        for k,v in result_paths.items():
            return_message += "%s\t%s" % (str(k), str(v.Path))
        raise ApplicationError(return_message)
