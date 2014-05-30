#!/usr/bin/env python
###############################################################################
#
# taxon_segregator.py - take a bam file and a taxonomy and segregate the reads
#                       into their relative taxanomic divisions
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
import sys
import os
import yaml
import time
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper
###############################################################################
###############################################################################
###############################################################################
###############################################################################
_DEFAULT_CONFIG = '''
---
pipeline:
    - demultiplex
    - overlap
    - map
    - filter
    - reduce
    - assemble
    - otu
    - taxonomy
threads: 1
log_level: INFO
log_file: null
output_directory: "."
input_raw_reads: # path to files from the Illumina machine
    - [null, null] # Each line contains the files for each sample. There should be two files per sample
    - [null, null]
pairtig_read_files: []
    #- null
    #- null
aliases: null # names to give your files in the output.  Must be the same length as 'input_raw_files', comment out the null on this line and uncomment the following lines
    #- sample1
    #- sample2
aligner_template: null
minimum_pairtig_length: 350 # minimum length of the overlapped pairs
pair_overlap_length: 30 # mimimum length of the overlap
mapper: bowtie # program used for read mapping
mapper_database: file.index # path to file of mapper
sam_file: initial_mapping.sam
taxonomy_file: file.tab # path to the file containing taxonomy for the reference sequences
mapping_similarity_cutoffs: [0.88, 0.92, 0.94, 0.98]
initial_mapping_similarity: 0.98 # the sequence similarity required between the reference database and the reads
taxon_coverage: [2, 1000] # list of two numbers. The first is the minimum coverage, the second is the number of bases that need to be covered
minimum_reconstruction_length: 1000 # minimum length of sequences that we define as 'full length'
otu_clustering_method: cdhit
otu_clustering_similarity: 0.99 # the similarity used for clustering full-length sequences from different samples into OTUs
read_mapping_percent: 0.90 # the percent identity that individual reads have to map with to be considered part of the reference
repset_output_file: full_length_sequences.repset.fa
assign_taxonomy_method: blast
minimum_taxon_similarity: 0.90 # sequences that fall below this cutoff will be listed as no taxonomy
'''
class AmplishotConfigError(Exception):
    def __init__(self, msg=''):
        self.msg = msg
    def __str__(self):
        return repr(self.msg)

class AmplishotProgramNotFoundError(AmplishotConfigError):
    def __init__(self, msg=''):
        super(AmplishotProgramNotFoundError, self).__init__(msg)

class AmplishotConfig(object):
    """ Class for reading the config file
        The Amplishot config file in written in YAML and contains all of the
        options for running the pipeline
    """

    def __init__(self, config=_DEFAULT_CONFIG):
        """ Load the YAML config file into an object
        """
        self.data = self._load(config)

    def __str__(self):
        """ dump the YAML to a string
        """
        return yaml.dump(self.data, Dumper=Dumper, default_flow_style=False,
                explicit_start=True)

    def _load(self, data):
        return yaml.load(data, Loader=Loader)

    def _check_program(self, prog):
        def is_exe(fpath):
            return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

        fpath, fname = os.path.split(program)
        if fpath:
            if is_exe(program):
                return True
        else:
            for path in os.environ["PATH"].split(os.pathsep):
                path = path.strip('"')
                exe_file = os.path.join(path, program)
                if is_exe(exe_file):
                    return True

        raise AmplishotProgramNotFoundError('Failed to find %s in your PATH. Is it installed?' % (prog,))

    def _check_required_programs(self):
        ''' check to make sure that programs called by subprocess are present
            in the users path before trying to run them
        '''
        if self.data['otu_clustering_method'] == 'cdhit':
            self._check_program('cd-hit')
        else:
            self._check_program(self.data['otu_clustering_method'])

        if self.data['mapper'] == 'bowtie':
            self._check_program('bowtie2')
        else:
            self._check_program(self.data['mapper'])

        if self.data['assembly_method'] == 'velvet':
            self._check_program('velveth')
        elif self.data['assembly_method'] == 'ray':
            self._check_program('Ray')
        else:
            self._check_program(self.data['assembly_method'])

        self._check_program('samtools')
        self._check_program('blastn')
        self._check_program('pear')


    def populate_from_config_file(self, fp):
        conf_file_data = self._load(fp)
        for key, value in conf_file_data.items():
            self.data[key] = value

    def populate_from_commandline(self, args):
        """ take the args that were input on the command line and overwrite the
            defaults
        """
        for key, value in vars(args).items():
            if value is not None and key in self.data:
                self.data[key] = value

    def write_config(self):
        # check for absolute paths; change them if they are not
        current_time = time.strftime('%Y%m%d%H%M%S')
        outfp = os.path.join(self.data['output_directory'], 'Amplishot_%s_config.yml' %
                current_time)
        with open(outfp, 'w') as fp:
            fp.write(str(self))

    def _set_path_to_absolute(self, fp, create=False):
        fp = os.path.expanduser(fp)
        fp = os.path.abspath(fp)
        if not os.path.exists(fp):
            if create:
                try:
                    os.makedirs(fp)
                except OSError:
                    raise AmplishotConfigError('cannot make the directory %s' % fp)
            else:
                raise AmplishotConfigError('cannot find file: %s' % fp)
        return fp

    def check_config_and_set_output(self, args):

        self.populate_from_commandline(args)
        if args.config is not None:
            self.populate_from_config_file(open(args.config))

        try:
            root_dir = self.data['output_directory'] =\
                self._set_path_to_absolute(self.data['output_directory'],
                        create=True)
        except KeyError:
            root_dir = self.data['output_directory'] = os.getcwd()

        if not os.path.exists(root_dir):
            os.makedirs(root_dir)

        try:
            self.data['nieghbours_file'] =\
                self._set_path_to_absolute(self.data['nieghbours_file'])
            self.data['taxonomy_file'] =\
                self._set_path_to_absolute(self.data['taxonomy_file'])
            self.data['mapper_database'] =\
                self._set_path_to_absolute(self.data['mapper_database'])
            self.data['blast_db'] =\
                self._set_path_to_absolute(self.data['blast_db'])
        except KeyError:
            raise AmplishotConfigError('one or more compulsory files have not \
                    been given on the command line or in the config file. \
                    Please check that you have defined the nieghbours_file, \
                    taxonomy_file, mapper_database, aligner_template and \
                    blast_db and that they point to valid files')

        # not compulsory files but if present fix the paths
        try:
            self.data['aligner_template'] =\
                self._set_path_to_absolute(self.data['aligner_template'])
        except KeyError:
            pass

        for i in range(len(self.data['input_raw_reads'])):
            f = self.data['input_raw_reads'][i]
            if isinstance(f, list):
                if len(f) != 2:
                    raise AmplishotConfigError('The value for each\
                            input_raw_read line must be either a single file\
                            path or a list of two files, one for each end of\
                            the fragment')
                f[0] = self._set_path_to_absolute(f[0])
                f[1] = self._set_path_to_absolute(f[1])
            else:
                f = self._set_path_to_absolute(f)

        self._check_required_programs()

        return True
###############################################################################
###############################################################################
###############################################################################
###############################################################################
