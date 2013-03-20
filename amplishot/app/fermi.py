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
__version__ = "0.3.3"
__maintainer__ = "Connor Skennerton"
__email__ = "c.skennerton@gmail.com"
__status__ = "Development"

###############################################################################
import os
from amplishot.app.util import ExtendedCommandLineApplication,\
        CompressedResultPath
from cogent.app.util import ValuedParameter, FlagParameter, ResultPath,\
    ApplicationError
import gzip

class Fermi(ExtendedCommandLineApplication):
    _parameters = {
            '-c': FlagParameter('-', 'c'),
            '-P': FlagParameter('-', 'P'),
            '-B': FlagParameter('-', 'B'),
            '-D': FlagParameter('-', 'D'),
            '-e': ValuedParameter('-', 'e', Delimiter=' ', IsPath=True),
            '-p': ValuedParameter('-', 'p', Delimiter=' ', Value='fmdef'),
            '-t': ValuedParameter('-', 't', Delimiter=' '),
            '-l': ValuedParameter('-', 'l', Delimiter=' '),
            '-k': ValuedParameter('-', 'k', Delimiter=' ')
            }
    _synonyms = {
            'paired': '-P',
            'interleaved': '-c',
            'split_build_indexing': '-D',
            'original_fmd_algorithm': '-B',
            'exec_name': '-e',
            'threads': '-t',
            'prefix': '-p',
            'trim_length': '-l',
            'kmer_length': '-k'
            }
    _command = 'run-fermi.pl'

    _input_handler = '_input_as_paths'

    def _decompress_file(self, infile):
        f_out = open(os.path.splitext(infile)[0], 'wb')
        f_in = gzip.open(infile, 'rb')
        f_out.writelines(f_in)
        f_out.close()
        f_in.close()
        return f_out.name
    
    def _get_result_paths(self, data):
        results = {}
        prefix = self.Parameters['-p'].Value
        if self.Parameters['-c'].isOn() or self.Parameters['-P'].isOn():
            results['scaffolds'] =\
            ResultPath(Path=self._decompress_file(os.path.join(self.WorkingDir,prefix +
                '.p5.fq.gz')), IsWritten=True)

        results['contigs'] = ResultPath(Path=self._decompress_file(os.path.join(self.WorkingDir,
            prefix + '.p2.mag.gz')), IsWritten=True)
        return results

    def _handle_app_result_build_failure(self, out, err, exit_status, result_paths):
        raise ApplicationError('exit status: %s\n result_paths:\n%s\n' %
                (str(exit_status), str(result_paths)))
