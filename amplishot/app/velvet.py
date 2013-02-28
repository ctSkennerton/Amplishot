#!/usr/bin/env python

from cogent.app.parameters import FlagParameter, ValuedParameter, ParameterError
from cogent.app.util import CommandLineApplication, ResultPath
from cogent.util.misc import app_path
from os.path import exists, join
import subprocess
import re


"""Application controller for velvet"""


__author__ = "Connor Skennerton"
__copyright__ = "Copyright 2013, Connor Skennerton"
__credits__ = ["Connor Skennerton"]
__license__ = "GPL"
__version__ = "0.3.0"
__maintainer__ = "Connor Skennerton"
__email__ = "c.skennerton@gmail.com"
__status__ = "Development"

class VelvetError(ParameterError):
    pass

class VelvethCategoryError(VelvetError):
    pass

class Velvet(CommandLineApplication):
    """Abstract base class for Velveth and Velvetg
       The function of this class is to overide the
       _error_on_missing_application function such that it also sets variables
       that are set during compile-time with velvet.  Two variables hold
       information about the max kmer length and the number of categories,
       which is subsequently used in determining parameters
    """
    _max_kmer_length = 31
    _max_categories = 2
    _regex = re.compile('CATEGORIES = (\d+)\nMAXKMERLENGTH = (\d+)')
    def __init__(self,params=None,InputHandler=None,SuppressStderr=None,\
        SuppressStdout=None, WorkingDir=None, TmpDir='/tmp', \
        TmpNameLen=20, HALT_EXEC=False):
        
        super(Velvet,self).__init__(params=params,
                    InputHandler=InputHandler, SuppressStderr=SuppressStderr,
                    SuppressStdout=SuppressStdout,  WorkingDir=WorkingDir,
                    TmpDir=TmpDir, TmpNameLen=TmpNameLen, HALT_EXEC=HALT_EXEC)

    def _error_on_missing_application(self, params):
        command = self._command
        if not (exists(command) or app_path(command)):
            raise ApplicationNotFoundError,\
             "Cannot find %s. Is it installed? Is it in your path?"\
             % command
        else:
            velvet_out = subprocess.Popen([command], stdout=subprocess.PIPE)
            text = velvet_out.communicate()[0]
            m = self._regex.search(text)
            if m:
                self._max_kmer_length = int(m.group(2))
                self._max_categories = int(m.group(1))
            else:
                raise VelvetError, "Cannot determine kmer length or categories\
                        from %s help message" % command
                

class Velveth(Velvet):
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

    _kmers = ""
    _parameters = {
        '-strand_specific': FlagParameter('-', 'strand_specific'),
        '-reuse_Sequences': FlagParameter('-', 'reuse_Sequences'),
        '-noHash': FlagParameter('-', 'noHash'),
        '-create_binary': FlagParameter('-', 'create_binary'),
        }
    _command = 'velveth'
    _accepted_formats = ['fasta', 'fastq', 'fasta.gz', 'fastq.gz', 'sam',
            'bam', 'eland', 'gerald', 'fmtAuto']
    _accepted_readtypes = ['short', 'shortPaired', 'long', 'longPaired',
            'reference']
    _accepted_layouts = ['separate', 'interleaved']
    
    def __init__(self,kmer, params=None,InputHandler=None,SuppressStderr=None,\
        SuppressStdout=None, WorkingDir=None, TmpDir='/tmp', \
        TmpNameLen=20, HALT_EXEC=False):
        super(Velveth,self).__init__(params=params,
                InputHandler=InputHandler, SuppressStderr=SuppressStderr,
                SuppressStdout=SuppressStdout,  WorkingDir=WorkingDir,
                TmpDir=TmpDir, TmpNameLen=TmpNameLen, HALT_EXEC=HALT_EXEC)
        if kmer > self._max_kmer_length:
            raise VelvetError, "the kmer length (%i) in greater than allowed by %s (%i)"\
                % (kmer, self._command, self._max_kmer_length)
        self.Kmer = kmer
        self._categories = []
        self._category_counts = {'short': 1, 'shortPaired': 1}

    def _accept_exit_status(self, exit_status):
        """Accept an exit status of 0 for the velveth program.
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
                    if self._category_counts[readtype] > self._max_categories:
                        raise VelvethCategoryError, "Adding more categories than allowed by %s %i > %i"\
                            % (self._command, self._category_counts[readtype], self._max_categories)
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


class Velvetg(Velvet):
    _parameters = {
        '-cov_cutoff': ValuedParameter('-', 'cov_cutoff', Delimiter=' '),
        '-ins_length': ValuedParameter('-', 'ins_length', Delimiter=' '),
        '-read_trkg': ValuedParameter('-', 'read_trkg', Delimiter=' '),
        '-min_contig_lgth': ValuedParameter('-', 'min_contig_lgth', Delimiter=' '),
        '-amos_file': ValuedParameter('-', 'amos_file', Delimiter=' '),
        '-exp_cov': ValuedParameter('-', 'exp_cov', Delimiter=' '),
        '-long_cov_cutoff': ValuedParameter('-', 'long_cov_cutoff', Delimiter=' '),
        '-ins_length_long': ValuedParameter('-', 'ins_length_long', Delimiter=' '),
        '-ins_length_long_sd': ValuedParameter('-', 'ins_length_long_sd', Delimiter=' '),
        '-ins_length_sd': ValuedParameter('-', 'ins_length_sd', Delimiter=' '),
        '-scaffolding': ValuedParameter('-', 'scaffolding', Delimiter=' '),        
        '-max_branch_length': ValuedParameter('-', 'max_branch_length', Delimiter=' '),
        '-max_divergence': ValuedParameter('-', 'max_divergence', Delimiter=' '),
        '-max_gap_count': ValuedParameter('-', 'max_gap_count', Delimiter=' '),
        '-min_pair_count': ValuedParameter('-', 'min_pair_count', Delimiter=' '),
        '-max_coverage': ValuedParameter('-', 'max_coverage', Delimiter=' '),
        '-coverage_mask': ValuedParameter('-', 'coverage_mask', Delimiter=' '),
        '-long_mult_cutoff': ValuedParameter('-', 'long_mult_cutoff', Delimiter=' '),
        '-unused_reads': ValuedParameter('-', 'unused_reads', Delimiter=' '),
        '-alignments': ValuedParameter('-', 'alignments', Delimiter=' '),
        '-exportFiltered': ValuedParameter('-', 'exportFiltered', Delimiter=' '),
        '-clean': ValuedParameter('-', 'clean', Delimiter=' '),
        '-very_clean': ValuedParameter('-', 'very_clean', Delimiter=' '),
        '-paired_exp_fraction': ValuedParameter('-', 'paired_exp_fraction', Delimiter=' '),
        '-conserveLong': ValuedParameter('-', 'conserveLong', Delimiter=' '),
        '-shortMatePaired': ValuedParameter('-', 'shortMatePaired', Delimiter=' '),
        }
    _synonyms = {
        'cov_cutoff': '-cov_cutoff',
        'ins_length': '-ins_length',
        'read_trkg': '-read_trkg',
        'min_contig_lgth': '-min_contig_lgth',
        'amos_file': '-amos_file',
        'exp_cov': '-exp_cov',
        'long_cov_cutoff': '-long_cov_cutoff',
        'ins_length_long': '-ins_length_long',
        'ins_length_long_sd': '-ins_length_long_sd',
        'ins_length_sd': '-ins_length_sd',
        'scaffolding': '-scaffolding',
        'max_branch_length': '-max_branch_length',
        'max_divergence': '-max_divergence',
        'max_gap_count': '-max_gap_count',
        'min_pair_count': '-min_pair_count',
        'max_coverage': '-max_coverage',
        'coverage_mask': '-coverage_mask',
        'long_mult_cutoff': '-long_mult_cutoff',
        'unused_reads': '-unused_reads',
        'alignments': '-alignments',
        'exportFiltered': '-exportFiltered',
        'clean': '-clean',
        'very_clean': '-very_clean',
        'paired_exp_fraction': '-paired_exp_fraction',
        'conserveLong': '-conserveLong',
        'shortMatePaired': '-shortMatePaired',
        }
    _command = 'velvetg'
    def __init__(self,params=None,InputHandler=None,SuppressStderr=None,\
        SuppressStdout=None, WorkingDir=None, TmpDir='/tmp', \
        TmpNameLen=20, HALT_EXEC=False):

        super(Velvetg,self).__init__(params=params,
                InputHandler=InputHandler, SuppressStderr=SuppressStderr,
                SuppressStdout=SuppressStdout,  WorkingDir=WorkingDir,
                TmpDir=TmpDir, TmpNameLen=TmpNameLen, HALT_EXEC=HALT_EXEC)
        for i in range(2,self._max_categories):
            il = '-ins_length'+str(i)
            sd = '-ins_length%i_sd' % i
            mp = '-shortMatePaired'+str(i)
            self._parameters[il] = ValuedParameter('-', il, Delimiter=' ')
            self._parameters[sd] = ValuedParameter('-', sd, Delimiter=' ')
            self._parameters[mp] = ValuedParameter('-', mp, Delimiter=' ')

    def _get_base_command(self):
        """ Returns the full command string 
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

        command_parts.append(self._command_delimiter.join(filter(\
            None,(map(str,parameters.values())))))
      
        return self._command_delimiter.join(command_parts).strip()

    BaseCommand = property(_get_base_command)


    def _accept_exit_status(self, exit_status):
        """Accept an exit status of 0 for the pandaseq program.
        """
        return exit_status == 0
    
    def _get_result_paths(self, data):
        result = {}
        out_files = {
                     'log': 'Log',
                     'contigs': 'contigs.fa', 
                     'stats': 'stats.txt'
                    }
        if self.Parameters['-amos_file'].isOn():
            out_files['amos'] = 'velvet_asm.afg'
        
        if self.Parameters['-unused_reads'].isOn():
            out_files['unused_reads'] = 'UnusedReads.fa'

        for name, outfile in out_files.items():
            result[name] = ResultPath(Path=join(self.WorkingDir,outfile),
                    IsWritten=True)

        return result


