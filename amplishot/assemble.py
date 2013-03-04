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
__version__ = "0.3.2"
__maintainer__ = "Connor Skennerton"
__email__ = "c.skennerton@gmail.com"
__status__ = "Development"

import logging
import multiprocessing
import os
from cogent.app.util import FilePath
from cogent.app.cd_hit import CD_HIT_EST
from amplishot.app.phrap import Phrap
from amplishot.app.fermi import Fermi
from amplishot.app.velvet import Velvetg, Velveth
import amplishot.parse.fastx
###############################################################################

def phrap_constructor( workingdir, params=None, infile='cdhitout.fa',
        suppressStdout=True, suppressStderr=True):
    #logger = multiprocessing.log_to_stderr()
    p = amplishot.app.phrap.Phrap(params=params, SuppressStdout=suppressStdout,
            SuppressStderr=suppressStderr, WorkingDir=workingdir,
            PrependPositionals=True, HALT_EXEC=False)
    p.add_positional_argument(FilePath(infile))
    logging.debug('%s (wd: %s)' % (str(p), workingdir))
    return p()


def fermi_constructor(workingdir, params=None, infile='reads.fa',
        suppressStdout=True, suppressStderr=True):
    name = os.path.splitext(infile)[0]
    params['prefix'] = name
    params['split_build_indexing'] = True
    params['original_fmd_algorithm'] = True
    f = Fermi(params=params, SuppressStdout=suppressStdout,
            SuppressStderr=suppressStderr, WorkingDir=workingdir,
            HALT_EXEC=False)
    # FIXME:
    # this is pretty doggy since fermi has to run make after so here I am just
    # injecting shell commands
    f.add_positional_argument(FilePath(infile))
    f.add_positional_argument('| make -f -')
    logging.debug('%s (wd: %s)' % (str(f), workingdir))
    return f()

def velvet_constructor(workingdir, params=None, infile='reads.fa',
        suppressStdout=True, suppressStderr=True):
    kmer_length = params['kmer_length']
    del params['kmer_length']
    name = os.path.splitext(infile)[0]
    workingdir = os.path.join(workingdir,name)
    vh = Velveth(kmer_length,  SuppressStdout=suppressStdout,
            SuppressStderr=suppressStderr, WorkingDir=workingdir,
            HALT_EXEC=False)
    vh.add_category(os.path.join('..',infile))
    vh()
    vg = Velvetg(params=params, SuppressStdout=suppressStdout,
            SuppressStderr=suppressStderr, WorkingDir=workingdir,
            HALT_EXEC=False)
    return vg()


def cd_hit_reduce(workingdir, infile='in.fa', outfile='cdhitout.fa', similarity=0.98, maxMemory=1000):
    cdhit = CD_HIT_EST(WorkingDir=workingdir)
    cdhit.Parameters['-i'].on(infile)
    cdhit.Parameters['-o'].on(outfile)
    cdhit.Parameters['-c'].on(similarity)
    cdhit.Parameters['-M'].on(maxMemory)
    logging.debug(cdhit.BaseCommand)
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


def reduce_and_assemble(taxon, assembler_params, cd_hit_params,
        assembler_constructor,
        suppressStderr=True, suppressStdout=True):
    cd_hit_reduce(taxon, **cd_hit_params)
    process_cd_hit_results(taxon, cdhitin=cd_hit_params['infile'],
            cdhitout=cd_hit_params['outfile'])

    return assemble(taxon, assembler_params, assembler_constructor,
            infile=cd_hit_params['outfile'],
            suppressStderr=suppressStderr,
                suppressStdout=suppressStdout)

def assemble(taxon, assembler_params,
        assembler_constructor, infile='reads.fa',
        suppressStderr=True, suppressStdout=True):
    results = assembler_constructor(taxon, assembler_params,
            infile=infile,
                suppressStderr=suppressStderr,
                suppressStdout=suppressStdout)
    return results['contigs'].name

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
        self.assembler_extension = 'fa'
        if assembler == 'fermi':
            self.constructor = fermi_constructor
        elif assembler == 'phrap':
            self.constructor = phrap_constructor
        elif assembler == 'velvet':
            self.constructor = velvet_constructor
        else:
            raise RuntimeError('your choice of assembler is not supported')

        self.assembler = assembler
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
        pool = multiprocessing.Pool(processes=self.config.data['threads'])
        pool_results = []
        try:
            assembly_params = self.config.data[self.assembler]
        except KeyError:
            assembly_params = None

        if self.reduce:
            # dataset reduction - cd-hit - each partition
            logging.info('clustering...')
            for taxon, cutoffs in taxons.items():
                for c in cutoffs:
                    infile_name = 'reads_%.2f.%s' % (c,
                            self.assembler_extension)
                    outfile_name = 'cdhitout_%.2f.fa' % c
                    
                    cd_hit_params = dict(infile=infile_name, 
                            outfile=outfile_name,
                            similarity=self.config.data['read_clustering_similarity'],
                            maxMemory=self.config.data['cdhit_max_memory'])

                    pool_results.append(pool.apply_async(reduce_and_assemble,
                        (taxon, assembly_params, cd_hit_params, self.constructor),
                        dict(suppressStderr=self.suppressStderr,
                        suppressStdout=self.suppressStdout)))

            pool.close()
            pool.join()

            full_length_counter = 1
            for result in pool_results:
                r = result.get()
                fxparser = amplishot.parse.fastx.FastxReader(open(r))
                for name, seq, qual in fxparser.parse(callback=amplishot.parse.fastx.greater_than, 
                length=self.config.data['minimum_reconstruction_length']):
                    seq_name = '%s_%i %s' % (sampleName, full_length_counter, name)
                    self.fullLengthSeqs[seq_name] = seq
                    full_length_counter += 1
        else:

            # generate overlaps - each partition
            for taxon, cutoffs in taxons.items():
                for c in cutoffs:
                    infile_name = 'reads_%.2f.%s' % (c,
                            self.assembler_extension)
                    try:
                        params = self.config.data[self.assembler]
                    except KeyError:
                        params = None
                    pool_results.append(pool.apply_async(assemble,
                        (taxon, assembly_params, self.constructor),
                        dict(suppressStderr=self.suppressStderr,
                        suppressStdout=self.suppressStdout, infile=infile_name)))
            pool.close()
            pool.join()

            full_length_counter = 1
            for result in pool_results:
                r = result.get()
                fxparser = amplishot.parse.fastx.FastxReader(open(r))
                for name, seq, qual in fxparser.parse(callback=amplishot.parse.fastx.greater_than, 
                length=self.config.data['minimum_reconstruction_length']):
                    seq_name = '%s_%i %s' % (sampleName, full_length_counter, name)
                    self.fullLengthSeqs[seq_name] = seq
                    full_length_counter += 1

    def write(self, fp):
        for name, seq in self.fullLengthSeqs.items():
            fp.write('>%s\n%s\n' %(name, seq))
