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


from Bio.Application import _Argument, _Switch, _Option, AbstractCommandline
import os

"""Application controller for phrap"""


__author__ = "Connor Skennerton"
__copyright__ = "Copyright 2012-2013, Connor Skennerton"
__credits__ = ["Connor Skennerton"]
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Connor Skennerton"
__email__ = "c.skennerton@gmail.com"
__status__ = "Prototype"


class Phrap(AbstractCommandline):
    """Simple phrap application controller.
    """
    def __init__(self, cmd='phrap', **kwargs):
        self.parameters = [

            _Argument(['infile', 'input', 'file'],
                'Input file containing reads in fasta format',
                filename=True, is_required=True),

            _Option(['-penalty', 'penalty'],'', equate=False),

            _Option(['-gap_init','gap_init'],'', equate=False),

            _Option(['-gap_ext','gap_ext'],'', equate=False),

            _Option(['-ins_gap_ext','ins_gap_ext'],'', equate=False),

            _Option(['-del_gap_ext','del_gap_ext'],'', equate=False),

            _Option(['-matrix','matrix'],'', equate=False),

            _Switch(['-raw','raw'],'',),

            _Option(['-minmatch','minmatch'],'', equate=False),

            _Option(['-maxmatch','maxmatch'],'', equate=False),

            _Option(['-max_group_size','max_group_size'],'', equate=False),

            _Switch(['-word_raw','word_raw'],'',),

            _Option(['-bandwidth','bandwidth'],'', equate=False),

            _Option(['-minscore','minscore'],'', equate=False),

            _Option(['-vector_bound','vector_bound'],'', equate=False),

            _Option(['-default_qual','default_qual'],'', equate=False),

            _Option(['-subclone_delim','subclone_delim'],'', equate=False),

            _Option(['-n_delim','n_delim'],'', equate=False),

            _Option(['-group_delim','group_delim'],'', equate=False),

            _Option(['-trim_start','trim_start'],'', equate=False),
        
            _Option(['-forcelevel','forcelevel'],'', equate=False),

            _Option(['-bypasslevel','bypasslevel'],'', equate=False),
        
            _Option(['-maxgap','maxgap'],'', equate=False),
        
            _Option(['-repeat_stringency','repeat_stringency'],'', equate=False),

            _Switch(['-revise_greedy','revise_greedy'],'',),

            _Switch(['-shatter_greedy','shatter_greedy'],'',),

            _Switch(['-preassemble','preassemble'],'',),

            _Switch(['-force_high','force_high'],'',),

            _Option(['-node_seg','node_seg'],'', equate=False),

            _Option(['-node_space','node_space'],'', equate=False),

            _Switch(['-tags','tags'],'',),

            _Switch(['-screen','screen'],'',),

            _Switch(['-old_ace','old_ace'],'',),

            _Switch(['-new_ace','new_ace'],'',),

            _Switch(['-ace','ace'],'',),

            _Switch(['-view','view'],'',),

            _Option(['-qual_show','qual_show'],'', equate=False),

            _Switch(['-print_extraneous_matches','print_extraneous_matches'],'',),

            _Switch(['-retain_duplicates','retain_duplicates'],'',),

            _Option(['-max_subclone_size','max_subclone_size'],'', equate=False),

            _Option(['-trim_penalty','trim_penalty'],'', equate=False),

            _Option(['-trim_score','trim_score'],'', equate=False),

            _Option(['-trim_qual','trim_qual'],'', equate=False),

            _Option(['-confirm_length','confirm_length'],'', equate=False),

            _Option(['-confirm_trim','confirm_trim'],'', equate=False),

            _Option(['-confirm_penalty','confirm_penalty'],'', equate=False),

            _Option(['-confirm_score','confirm_score'],'', equate=False),

            _Option(['-indexwordsize','indexwordsize'],'', equate=False),

       ]
        AbstractCommandline.__init__(self, cmd, **kwargs)
    def get_result_paths(self, cwd=''):
        results = dict()
        results['contigs'] = os.path.join(cwd,self.file +
               '.contigs')
        results['contigs_qual'] = os.path.join(cwd,self.file +
               '.contigs.qual')
        results['log'] = os.path.join(cwd,self.file +
               '.log')
        results['problems'] = os.path.join(cwd,self.file +
               '.problems')
        results['problems_qual'] = os.path.join(cwd,self.file +
               '.problems.qual')
        results['singlets'] = os.path.join(cwd,self.file +
               '.siglets')
        return results
