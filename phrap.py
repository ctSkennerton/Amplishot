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

            '-indexwordsize': ValuedParameter('-', 'indexwordsize', Delimiter=' '), 

               }
    _input_handler = '_input_as_path'
    _command = 'phrap'

    def __call__(self,data=None, remove_tmp=True):
        """Run the application with the specified kwargs on data
        
            data: anything that can be cast into a string or written out to
                a file. Usually either a list of things or a single string or 
                number. input_handler will be called on this data before it 
                is passed as part of the command-line argument, so by creating
                your own input handlers you can customize what kind of data
                you want your application to accept

            remove_tmp: if True, removes tmp files
        """
        input_handler = self.InputHandler
        suppress_stdout = self.SuppressStdout
        suppress_stderr = self.SuppressStderr
        if suppress_stdout:
            outfile = FilePath('/dev/null')
        else:
            outfile = self.getTmpFilename(self.TmpDir)
        if suppress_stderr:
            errfile = FilePath('/dev/null')
        else:
            errfile = FilePath(self.getTmpFilename(self.TmpDir))
        if data is None:
            input_arg = ''
        else:
            input_arg = getattr(self,input_handler)(data)

        # Build up the command, consisting of a BaseCommand followed by
        # input and output (file) specifications
        command = self._command_delimiter.join(filter(None,\
            [str(input_arg), self.BaseCommand,'>',str(outfile),'2>',\
                str(errfile)]))
        if self.HaltExec: 
            raise AssertionError, "Halted exec with command:\n" + command
        # The return value of system is a 16-bit number containing the signal 
        # number that killed the process, and then the exit status. 
        # We only want to keep the exit status so do a right bitwise shift to 
        # get rid of the signal number byte
        exit_status = system(command) >> 8
        
        # Determine if error should be raised due to exit status of 
        # appliciation
        if not self._accept_exit_status(exit_status):
            raise ApplicationError, \
             'Unacceptable application exit status: %s\n' % str(exit_status) +\
             'Command:\n%s\nStdOut:\n%s\nStdErr:\n%s\n' % (command, 
                                                           open(outfile).read(), 
                                                           open(errfile).read())
        
        # open the stdout and stderr if not being suppressed
        out = None
        if not suppress_stdout:
            out = open(outfile,"r")
        err = None        
        if not suppress_stderr:
            err = open(errfile,"r")
            
        try:
            result = CommandLineAppResult(\
             out,err,exit_status,result_paths=self._get_result_paths(data))
        except ApplicationError:
            result = self._handle_app_result_build_failure(\
             out,err,exit_status,self._get_result_paths(data))

        # Clean up the input file if one was created
        if remove_tmp:
            if self._input_filename:
                remove(self._input_filename)
                self._input_filename = None

        return result

    def _accept_exit_status(self, exit_status):
        """Accept an exit status of 0 for the phrap program.
        """
        return exit_status == 0
    


