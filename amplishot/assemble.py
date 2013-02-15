###############################################################################
#
# assemble.py - assembly functions used in amplishot
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
__version__ = "0.1.4"
__maintainer__ = "Connor Skennerton"
__email__ = "c.skennerton@gmail.com"
__status__ = "Development"

import logging
import multiprocessing
import os
from cogent.app.util import FilePath
from cogent.app.cd_hit import CD_HIT_EST
from amplishot.app.phrap_pycogent import Phrap
import amplishot.parse.fastx
###############################################################################

def phrap_constructor( workingdir, params=None, infile='cdhitout.fa',
        suppressStdout=True, suppressStderr=True):
    p = amplishot.app.phrap_pycogent.Phrap(params=params, SuppressStdout=suppressStdout,
            SuppressStderr=suppressStderr, WorkingDir=workingdir,
            PrependPositionals=True)
   # p.Parameters['-minscore'].on(300)
   # p.Parameters['-ace'].on()
   # p.Parameters['-penalty'].on(-9)
   # p.Parameters['-gap_init'].on(-11)
   # p.Parameters['-gap_ext'].on(-10)
    p.add_positional_argument(FilePath(infile))
    return p


def cd_hit_reduce(workingdir, infile='in.fa', outfile='cdhitout.fa', similarity=0.98, maxMemory=1000):
    cdhit = CD_HIT_EST(WorkingDir=workingdir)
    cdhit.Parameters['-i'].on(infile)
    cdhit.Parameters['-o'].on(outfile)
    cdhit.Parameters['-c'].on(similarity)
    cdhit.Parameters['-M'].on(maxMemory)
    cdhit()


def process_cd_hit_results(directory, cdhitout='cdhitout.fa', cdhitin='in.fa'):
    clustered_qual = open(os.path.join(directory,cdhitout+'.qual'), 'w')
    clustered_reads = set()
    cd_hit_fasta = open(os.path.join(directory,cdhitout))
    fxparser = amplishot.parse.fastx.FastxReader(cd_hit_fasta)
    for name, seq, qual in fxparser.parse():
        clustered_reads.add(name)

    cd_hit_fasta.close()
    inqual = open(os.path.join(directory,cdhitin+'.qual'))
    qualparser = amplishot.parse.fastx.QualityReader(inqual)
    for name, qual in qualparser.parse():
        if name in clustered_reads:
            clustered_qual.write('>%s\n%s\n' % (name, qual))
    inqual.close()

class AssemblyWrapper(object):
    def __init__(self, assembler, config, preAssembleReduction=True,
            stdout=None, stderr=None):
        ''' Entry point for all denovo assembly methods
            This method will pass off onto the correct assembly method and perform
            the required cleanup and format conversions required for each assembler

            assembler: A string representing the assembler of choice
            taxonDirs: A list of paths to directories containing taxons for
            assembly
            config: amplishot config object
            preAssembleReduction: set to True if the reads in each taxon should be
            clustered before assembly
            stdout: specify a python file like object to output the stdout stream
            from the assembler.  Set to False to suppress stdout, when set to none
            a temporary file will be made
            stderr: specify a python file like object to output the stderr stream
            from the assembler.  Set to False to suppress stderr, when set to none
            a temporary file will be made
        '''
        self.reduce = preAssembleReduction
        self.config = config
        self.fullLengthSeqs = dict()
        if assembler != 'phrap':
            raise RuntimeError('only phrap is supported as an assember at the'\
                    'moment. your choice was: %s' % assembler)
        else:
            self.assembler = assembler
            self.constructor = phrap_constructor
            if stdout is False:
                self.suppressStdout = True
            else:
                self.suppressStdout = False

            if stderr is False:
                self.suppressStderr = True
            else:
                self.suppressStderr = False

        
    def __call__(self, taxons, sampleName):
        infile_name = 'reads.fa'
        if self.reduce:
            # dataset reduction - cd-hit - each partition
            pool = multiprocessing.Pool(processes=self.config.data['threads'])
            logging.info('clustering...')
            for t in taxons:
                pool.apply_async(cd_hit_reduce, (t,), dict(infile=infile_name, 
                    similarity=self.config.data['read_clustering_similarity'],
                    maxMemory=self.config.data['cdhit_max_memory']))
            pool.close()
            pool.join()
            logging.info('assembling...')
            for t in taxons:
                process_cd_hit_results(t, cdhitin=infile_name)

            infile_name = 'cdhitout.fa'

        ## generate overlaps - each partition
        full_length_counter = 1
        for t in taxons:
            try:
                params = self.config.data[self.assembler]
            except KeyError:
                params = None
            inst = self.constructor(t, params, infile=infile_name,
                    suppressStderr=self.suppressStderr,
                    suppressStdout=self.suppressStdout)
            logging.debug('%s (wd: %s)', str(inst), t)
            results = inst()
            # collect all full-length sequences
            fxparser = amplishot.parse.fastx.FastxReader(results['contigs'])
            for name, seq, qual in fxparser.parse(callback=amplishot.parse.fastx.greater_than, 
            length=self.config.data['minimum_reconstruction_length']):
                seq_name = '%s_%i %s' % (sampleName, full_length_counter, name)
                self.fullLengthSeqs[seq_name] = seq
                full_length_counter += 1

    def write(self, fp):
        for name, seq in self.fullLengthSeqs.items():
            fp.write('>%s\n%s\n' %(name, seq))
