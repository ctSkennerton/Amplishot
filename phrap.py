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

            '-penalty': ValuedParameter('-', 'penalty', Delimeter=' '), #Mismatch (substitution) penalty for SWAT comparisons.  

            '-gap_init': ValuedParameter('-', 'gap_init', Delimeter=' '),# Gap initiation penalty for SWAT comparisons.           

            '-gap_ext': ValuedParameter('-', 'gap_ext', Delimeter=' '), # Gap extension penalty for SWAT comparisons.            

            '-ins_gap_ext': ValuedParameter('-', 'ins_gap_ext', Delimeter=' '), 

            '-del_gap_ext': ValuedParameter('-', 'del_gap_ext', Delimeter=' '), 

            '-matrix': ValuedParameter('-', 'matrix', Delimeter=' '), 

            '-raw': FlagParameter('-', 'raw'), 

            '-minmatch': ValuedParameter('-', 'minmatch', Delimeter=' '), 

            '-maxmatch': ValuedParameter('-', 'maxmatch', Delimeter=' '), 

            '-max_group_size': ValuedParameter('-', 'max_group_size', Delimeter=' '), 

            '-word_raw': FlagParameter('-', 'word_raw'), 

            '-bandwidth': ValuedParameter('-', 'bandwidth', Delimeter=' '), 

            '-minscore': ValuedParameter('-', 'minscore', Delimeter=' '), 

            '-vector_bound': ValuedParameter('-', 'vector_bound', Delimeter=' '), 

            '-default_qual': ValuedParameter('-', 'default_qual', Delimeter=' '), 

            '-subclone_delim': ValuedParameter('-', 'subclone_delim', Delimeter=' '), 

            '-n_delim': ValuedParameter('-', 'n_delim', Delimeter=' '), 

            '-group_delim': ValuedParameter('-', 'group_delim', Delimeter=' '), 

            '-trim_start': ValuedParameter('-', 'trim_start', Delimeter=' '), 
            
            '-forcelevel': ValuedParameter('-', 'forcelevel', Delimeter=' '), 

            '-bypasslevel': ValuedParameter('-', 'bypasslevel', Delimeter=' '), 
            
            '-maxgap': ValuedParameter('-', 'maxgap', Delimeter=' '), 
            
            '-repeat_stringency': ValuedParameter('-', 'repeat_stringency', Delimeter=' '), 

            '-revise_greedy': FlagParameter('-', 'revise_greedy'), 

            '-shatter_greedy': FlagParameter('-', 'shatter_greedy'), 

            '-preassemble': FlagParameter('-', 'preassemble'), 

            '-force_high': FlagParameter('-', 'force_high'), 

            '-node_seg': ValuedParameter('-', 'node_seg', Delimeter=' '), 

            '-node_space': ValuedParameter('-', 'node_space', Delimeter=' '), 

            '-tags': FlagParameter('-', 'tags'), 

            '-screen': FlagParameter('-', 'screen'), 

            '-old_ace': FlagParameter('-', 'old_ace'), 

            '-new_ace': FlagParameter('-', 'new_ace'), 

            '-ace': FlagParameter('-', 'ace'), 

            '-view': FlagParameter('-', 'view'), 

            '-qual_show': ValuedParameter('-', 'qual_show', Delimeter=' '), 

            '-print_extraneous_matches': FlagParameter('-', 'print_extraneous_matches'), 

            '-retain_duplicates': FlagParameter('-', 'retain_duplicates'), 

            '-max_subclone_size': ValuedParameter('-', 'max_subclone_size', Delimeter=' '), 

            '-trim_penalty': ValuedParameter('-', 'trim_penalty', Delimeter=' '), 

            '-trim_score': ValuedParameter('-', 'trim_score', Delimeter=' '), 

            '-trim_qual': ValuedParameter('-', 'trim_qual', Delimeter=' '), 

            '-confirm_length': ValuedParameter('-', 'confirm_length', Delimeter=' '), 

            '-confirm_trim': ValuedParameter('-', 'confirm_trim', Delimeter=' '), 

            '-confirm_penalty': ValuedParameter('-', 'confirm_penalty', Delimeter=' '), 

            '-confirm_score': ValuedParameter('-', 'confirm_score', Delimeter=' '), 

            '-indexwordsize': ValuedParameter('-', 'indexwordsize', Delimeter=' '), 

               }
    _input_handler = '_input_as_parameter'
    _command = 'phrap'

    def _accept_exit_status(self, exit_status):
        """Accept an exit status of 0 for the phrap program.
        """
        return exit_status == 0
    
    def _input_as_parameter(self,forward,reverse):
        """ Set -f & -r to the forward and reverse reads files
        """
        self.Parameters['-f'].on(forward)
        self.Parameters['-r'].on(reverse)
        return ''

