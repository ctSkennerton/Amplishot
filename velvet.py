#!/usr/bin/env python

from cogent.app.parameters import FlagParameter, ValuedParameter, ParameterError
from cogent.app.util import CommandLineApplication


"""Application controller for velvet"""


__author__ = "Connor Skennerton"
__copyright__ = "Copyright 2013, Connor Skennerton"
__credits__ = ["Connor Skennerton"]
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Connor Skennerton"
__email__ = "c.skennerton@gmail.com"
__status__ = "Development"

class VelvethCategoryError(ParameterError):
    pass


class Velveth(CommandLineApplication):
    """Simple velveth application controller
       Velveth uses an interesting scheme for its command line parameters
       whereby there are 'categories'.  each category has the same options
       available to it, however the number of categories can be infinate, but
       determined at compile time.  Furthermore runtime flags carryover to any
       category, until they are overriden.  

       Example:
       vleveth outdir 63 -fastq -short in.fq -shortPaired in_paired.fq
       | cmd  |      |  |   category 1      |     category 2          |

       In the example above both category 1 and category 2 will share the
       -fastq option

       this kind of option structure doesn't lend itself to a singular mapping
       of options; instead each category will be held in a list that contains a
       mapping of each of the category's options
    """

    _categories = []
    _category_counts = {'short': 1, 'shortPaired': 1}
    _kmers = ""
    _parameters = {
        '-strand_specific': FlagParameter('-', 'strand_specific'),
        '-reuse_Sequences': FlagParameter('-', 'reuse_Sequences'),
        '-noHash': FlagParameter('-', 'noHash'),
        '-create_binary': FlagParameter('-', 'create_binary'),
        }
    _command = 'velveth'
    _accepted_formats = ["fasta", "fastq", "fasta.gz", "fastq.gz", "sam",
            "bam", "eland", "gerald"]
    _accepted_readtypes = ["short", "shortPaired", "long", "longPaired",
            "reference"]
    _accepted_layouts = ["separate", "interleaved"]
    
    def __init__(self,kmer, params=None,InputHandler=None,SuppressStderr=None,\
        SuppressStdout=None, WorkingDir=None, TmpDir='/tmp', \
        TmpNameLen=20, HALT_EXEC=False):
        self.Kmer = kmer
        super(Velveth,self).__init__(params=params,
                InputHandler=InputHandler, SuppressStderr=SuppressStderr,
                SuppressStdout=SuppressStdout,  WorkingDir=WorkingDir,
                TmpDir=TmpDir, TmpNameLen=TmpNameLen, HALT_EXEC=HALT_EXEC)

    def _accept_exit_status(self, exit_status):
        """Accept an exit status of 0 for the velvet program.
        """
        return exit_status == 0
    
    def _get_base_command(self):
        """ Returns the full command string 
            Goes through all of the categories and prints them out
        """
        command_parts = []
        # Append a change directory to the beginning of the command to change 
        # to self.WorkingDir before running the command
        # WorkingDir should be in quotes -- filenames might contain spaces
        cd_command = ''.join(['cd ',str(self.WorkingDir),';'])
        if self._command is None:
            raise ApplicationError, '_command has not been set.'
        command = self._command
        parameters = self.Parameters
        
        command_parts.append(cd_command)
        command_parts.append(command)
        command_parts.append(str(self.WorkingDir))
        command_parts.append(str(self.Kmer))
        for category in self._categories:
            print category
            command_parts.append(self._command_delimiter.join(filter(\
                    None,(map(str, category)))))

        command_parts.append(self._command_delimiter.join(filter(\
            None,(map(str,parameters.values())))))
      
        return self._command_delimiter.join(command_parts).strip()

    BaseCommand = property(_get_base_command)

    def add_category(self, infiles, informat=None, readtype=None,
            layout=None):
        category_opts = []
        # check informat 
        if informat is None:
            category_opts.append(FlagParameter('-', informat))
        elif informat not in self._accepted_formats:
            raise VelvethCategoryError
        else:
            category_opts.append(FlagParameter('-', informat,Value=True))

        # check readtype
        if readtype is None:
            category_opts.append(FlagParameter('-', readtype))
        elif readtype not in self._accepted_readtypes:
            raise VelvethCategoryError
        else:
            try:
                if self._category_counts[readtype] > 1:
                    count = self._category_counts[readtype]
                    self._category_counts[readtype] += 1
                    readtype += str(count)
                else:
                    self._category_counts[readtype] += 1
            finally:
                category_opts.append(FlagParameter('-', readtype, Value=True))


        #check layout
        if layout is None:
            category_opts.append(FlagParameter('-', layout))
        elif layout not in self._accepted_layouts:
            raise VelvethCategoryError
        else:
            category_opts.append(FlagParameter('-', layout, Value=True))

        # add in the files
        if isinstance(infiles,list):
            category_opts.extend(infiles)
        else:
            category_opts.append(infiles)

        self._categories.append(category_opts)
