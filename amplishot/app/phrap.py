#!/usr/bin/env python
###############################################################################
#
# phrap.py - application controller for phrap
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
__version__ = "0.2.2"
__maintainer__ = "Connor Skennerton"
__email__ = "c.skennerton@gmail.com"
__status__ = "Development"

###############################################################################
import os
from amplishot.app.util import ExtendedCommandLineApplication
from cogent.app.util import ApplicationError, ValuedParameter, FlagParameter,\
        FilePath, ResultPath

class Phrap(ExtendedCommandLineApplication):
    """Simple phrap application controller.
    """
    _parameters = {
            #'infile': ValuedParameter('','', Delimiter='', IsPath=True),
            
            '-penalty': ValuedParameter('-', 'penalty', Delimiter=' '), #Mismatch (substitution) penalty for SWAT comparisons.  

            '-gap_init': ValuedParameter('-', 'gap_init', Delimiter=' '),# Gap initiation penalty for SWAT comparisons.           

            '-gap_ext': ValuedParameter('-', 'gap_ext', Delimiter=' '), # Gap extension penalty for SWAT comparisons.            

            '-ins_gap_ext': ValuedParameter('-', 'ins_gap_ext', Delimiter=' '), 

            '-del_gap_ext': ValuedParameter('-', 'del_gap_ext', Delimiter=' '), 

            '-matrix': ValuedParameter('-', 'matrix', Delimiter=' '), 

            '-raw': FlagParameter('-', 'raw'), 

            '-minmatch': ValuedParameter('-', 'minmatch', Delimiter=' '), 

            '-maxmatch': ValuedParameter('-', 'maxmatch', Delimiter=' '), 

            '-max_group_size': ValuedParameter('-', 'max_group_size', Delimiter=' '), 

            '-word_raw': FlagParameter('-', 'word_raw'), 

            '-bandwidth': ValuedParameter('-', 'bandwidth', Delimiter=' '), 

            '-minscore': ValuedParameter('-', 'minscore', Delimiter=' '), 

            '-vector_bound': ValuedParameter('-', 'vector_bound', Delimiter=' '), 

            '-default_qual': ValuedParameter('-', 'default_qual', Delimiter=' '), 

            '-subclone_delim': ValuedParameter('-', 'subclone_delim', Delimiter=' '), 

            '-n_delim': ValuedParameter('-', 'n_delim', Delimiter=' '), 

            '-group_delim': ValuedParameter('-', 'group_delim', Delimiter=' '), 

            '-trim_start': ValuedParameter('-', 'trim_start', Delimiter=' '), 
            
            '-forcelevel': ValuedParameter('-', 'forcelevel', Delimiter=' '), 

            '-bypasslevel': ValuedParameter('-', 'bypasslevel', Delimiter=' '), 
            
            '-maxgap': ValuedParameter('-', 'maxgap', Delimiter=' '), 
            
            '-repeat_stringency': ValuedParameter('-', 'repeat_stringency', Delimiter=' '), 

            '-revise_greedy': FlagParameter('-', 'revise_greedy'), 

            '-shatter_greedy': FlagParameter('-', 'shatter_greedy'), 

            '-preassemble': FlagParameter('-', 'preassemble'), 

            '-force_high': FlagParameter('-', 'force_high'), 

            '-node_seg': ValuedParameter('-', 'node_seg', Delimiter=' '), 

            '-node_space': ValuedParameter('-', 'node_space', Delimiter=' '), 

            '-tags': FlagParameter('-', 'tags'), 

            '-screen': FlagParameter('-', 'screen'), 

            '-old_ace': FlagParameter('-', 'old_ace'), 

            '-new_ace': FlagParameter('-', 'new_ace'), 

            '-ace': FlagParameter('-', 'ace'), 

            '-view': FlagParameter('-', 'view'), 

            '-qual_show': ValuedParameter('-', 'qual_show', Delimiter=' '), 

            '-print_extraneous_matches': FlagParameter('-', 'print_extraneous_matches'), 

            '-retain_duplicates': FlagParameter('-', 'retain_duplicates'), 

            '-max_subclone_size': ValuedParameter('-', 'max_subclone_size', Delimiter=' '), 

            '-trim_penalty': ValuedParameter('-', 'trim_penalty', Delimiter=' '), 

            '-trim_score': ValuedParameter('-', 'trim_score', Delimiter=' '), 

            '-trim_qual': ValuedParameter('-', 'trim_qual', Delimiter=' '), 

            '-confirm_length': ValuedParameter('-', 'confirm_length', Delimiter=' '), 

            '-confirm_trim': ValuedParameter('-', 'confirm_trim', Delimiter=' '), 

            '-confirm_penalty': ValuedParameter('-', 'confirm_penalty', Delimiter=' '), 

            '-confirm_score': ValuedParameter('-', 'confirm_score', Delimiter=' '), 

            '-indexwordsize': ValuedParameter('-', 'indexwordsize', Delimiter=' ')

               }
    _command = 'phrap'

    def _input_as_string(self, data):
        raise ApplicationError("Input data cannot be specified in this way"\
        "please use the prepend function")

    def _input_as_multiline_string(self, data):
        raise ApplicationError("Input data cannot be specified in this way"\
        "please use the prepend function")

    def _input_as_lines(self, data):
        raise ApplicationError("Input data cannot be specified in this way"\
        "please use the prepend function")

    def _input_as_path(self, data):
        raise ApplicationError("Input data cannot be specified in this way"\
        "please use the prepend function")

    def _input_as_paths(self, data):
        raise ApplicationError("Input data cannot be specified in this way"\
        "please use the prepend function")

    def _accept_exit_status(self, exit_status):
        """Accept an exit status of 0 for the phrap program.
        """
        return exit_status == 0
    
    def _get_result_paths(self, data):
        oprefix = self._positionals[0]
        results = {
                'contigs': ResultPath(Path=os.path.join(self.WorkingDir,
                    oprefix + '.contigs'), IsWritten=True),
                'contigs_qual': ResultPath(Path=os.path.join(self.WorkingDir, oprefix +\
                    '.contigs.qual'), IsWritten=True),
                'log': ResultPath(Path=os.path.join(self.WorkingDir, oprefix +\
                    '.log'), IsWritten=True),
                'problems': ResultPath(Path=os.path.join(self.WorkingDir,
                    oprefix + '.problems'), IsWritten=True),
                'problems_qual': ResultPath(Path=os.path.join(self.WorkingDir,
                    oprefix + '.problems.qual'),
                    IsWritten=True),
                'singlets': ResultPath(Path=os.path.join(self.WorkingDir,
                    oprefix + '.singlets'), IsWritten=True)
                }
        return results
