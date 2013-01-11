#!/usr/bin/env python

from cogent.app.parameters import FlagParameter, ValuedParameter
from cogent.app.util import CommandLineApplication


"""Application controller for phrap"""


__author__ = "Connor Skennerton"
__copyright__ = "Copyright 2012-2013, Connor Skennerton"
__credits__ = ["Connor Skennerton"]
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Connor Skennerton"
__email__ = "c.skennerton@gmail.com"
__status__ = "Prototype"


class Phrap(CommandLineApplication):
    """Simple phrap application controller.
    """
    _parameters = {

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
    _input_handler = '_input_as_path'
    _command = 'phrap'


    def _accept_exit_status(self, exit_status):
        """Accept an exit status of 0 for the phrap program.
        """
        return exit_status == 0
    
    def _get_result_paths(self, data):
        oprefix = getattr(self,input_handler)(data)
        results = {
                'contigs': ResultPath(Path=oprefix+'.contigs', IsWritten=True),
                'contigsQual': ResultPath(Path=oprefix+'.contigs.qual',
                    IsWritten=True),
                'log': ResultPath(Path=oprefix+'.log', IsWritten=True),
                'problems': ResultPath(Path=oprefix+'.problems',
                    IsWritten=True),
                'problemsQual': ResultPath(Path=oprefix+'.problems.qual',
                    IsWritten=True),
                'singlets': ResultPath(Path=oprefix+'.singlets', IsWritten=True)
                }
        return results


