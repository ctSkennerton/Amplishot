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
__version__ = "0.4.0"
__maintainer__ = "Connor Skennerton"
__email__ = "c.skennerton@gmail.com"
__status__ = "Development"

import logging
import multiprocessing
import os
import glob
import re
import subprocess
from cogent.app.util import FilePath
from amplishot.app.cd_hit import CD_HIT_EST
from amplishot.app.phrap import Phrap
from amplishot.app.fermi import Fermi
from amplishot.app.velvet import Velvetg, Velveth
from amplishot.app.spades import Spades
from amplishot.app.ray import Ray
import amplishot.parse.fastx
###############################################################################

def phrap_constructor( workingdir, params=None,
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

def spades_constructor(workingdir, params=None,
        suppressStdout=True, suppressStderr=True):
    if '-o' not in params:
        params['-o'] = infile[:-3]

    if 'r1' in params and 'r2' in params:
        params['-1'] = params['r1']
        params['-2'] = params['r2']
        del params['r1']
        del params['r2']
    elif 'r12' in params:
        params['--12'] = params['r12']
        del params['r12']
    
    if 's' in params:
        params['-s'] = params['s']
        del params['s']

    s = Spades(params=params, SuppressStdout=suppressStdout,
            SuppressStderr=suppressStderr, WorkingDir=workingdir,
            HALT_EXEC=False)
    return s()

def ray_constructor(workingdir, params=None,
        suppressStdout=True, suppressStderr=True):
    if '-o' not in params:
        params['-o'] = infile[:-3]

    if 'kmer_size' not in params:
        params['-k'] = 31
    else:
        params['-k'] = params['kmer_size']
        del params['kmer_size']

    if 'r1' in params and 'r2' in params:
        params['-p'] = os.path.join(workingdir,params['r1']) + " " + os.path.join(workingdir, params['r2'])
        del params['r1']
        del params['r2']
    elif 'r12' in params:
        params['-i'] = os.path.join(workingdir,params['r12'])
        del params['r12']
    
    if 's' in params:
        params['-s'] = os.path.join(workingdir,params['s'])
        del params['s']
    cmd = ['Ray']
    cmd.extend([str(k) + " " + str(v) for k,v in params.items()])
    logging.debug(' '.join(cmd))
    with open(os.devnull, 'w') as dn:
        retcode = subprocess.call(' '.join(cmd), stdout=dn, stderr=dn, shell=True)
    if retcode == 0:
        ret = {}
        ret['contigs'] = os.path.join(params['-o'],"Contigs.fasta")
        return ret
    else:
        raise RuntimeError("Ray failed to run. exit code = %d\ncmd: %s" % (retcode, ' '.join(cmd)))
    #r = Ray(params=params, SuppressStdout=suppressStdout,
    #        SuppressStderr=suppressStderr, WorkingDir=workingdir,
    #        HALT_EXEC=False)
    #logging.debug(str(r))
    #return r()

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

def assemble(filePrefix, assembler_params,
        assembler_constructor,
        suppressStderr=True, suppressStdout=True):
    working_dir, file_prefix = os.path.split(filePrefix)
    if 'r1' in assembler_params and 'r2' in assembler_params:
        assembler_params['r1'] = file_prefix + assembler_params['r1']
        assembler_params['r2'] = file_prefix + assembler_params['r2']
    elif 'r12' in assembler_params:
        assembler_params['r12'] = file_prefix + assembler_params['r12']
    
    if 's' in assembler_params:
        assembler_params['s'] = file_prefix + assembler_params['s']

    results = assembler_constructor(working_dir, assembler_params,
                suppressStderr=suppressStderr,
                suppressStdout=suppressStdout)
    try:
        return results['contigs'].name
    except:
        return results['contigs']

class AssemblyWrapper(object):
    singles_re = re.compile('(\d+_0\.\d+)_s\.fast[aq]$')
    pairs_re = re.compile('(\d+_0\.\d+)_R12\.fast[aq]$')
    first_re = re.compile('(\d+_0\.\d+)_R1\.fast[aq]$')

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
        self.fullLengthQuals = dict()
        self.assembler_extension = 'fa'
        if assembler == 'fermi':
            self.constructor = fermi_constructor
        elif assembler == 'phrap':
            self.constructor = phrap_constructor
        elif assembler == 'velvet':
            self.constructor = velvet_constructor
        elif assembler == 'spades':
            self.constructor = spades_constructor
        elif assembler == 'ray':
            self.constructor = ray_constructor
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
            assembly_params = dict()

        if self.reduce:
            # dataset reduction - cd-hit - each partition
            logging.info('clustering...')
            for filepath, cutoffs in taxons.items():
                for read_file in os.listdir(filepath):  #glob.iglob("*fa"):
                    match = self.read_re.search(read_file)
                    if match:
                        infile_name = read_file
                        outfile_name = 'cdhitout_%s.fa' % match.group(1)
                        logging.debug("Found read file: %s\nIn directory: %s",
                                infile_name, filepath)
                    #for c in cutoffs:
                        #infile_name = 'reads_%.2f.%s' % (c,
                        #        self.assembler_extension)
                        #outfile_name = 'cdhitout_%.2f.fa' % c
                        
                        cd_hit_params = dict(infile=infile_name, 
                                outfile=outfile_name,
                                similarity=self.config.data['read_clustering_similarity'],
                                maxMemory=self.config.data['cdhit_max_memory'])

                        pool_results.append(pool.apply_async(reduce_and_assemble,
                            (filepath, assembly_params, cd_hit_params, self.constructor),
                            dict(suppressStderr=self.suppressStderr,
                            suppressStdout=self.suppressStdout)))

            pool.close()
            pool.join()

        else:
            for filepath, read_types in taxons.items():
                local_params = assembly_params.copy()
                local_params.update(read_types)
                local_params['-o'] = filepath

                #single thread mode
                #res = assemble(filepath, assembly_params, self.constructor, suppressStderr=self.suppressStderr, suppressStdout=self.suppressStdout)
                #pool_results.append(res)
                
                #FIXME: There is something super weird with the multithreading
                #mode that causes an error when reading an output file even when
                #the file exists.  This problem does not happen in the single
                #threaded mode.

                pool_results.append(pool.apply_async(assemble,
                    (filepath, local_params, self.constructor),
                    dict(suppressStderr=self.suppressStderr,
                    suppressStdout=self.suppressStdout)))
            pool.close()
            pool.join()

        full_length_counter = 1
        #for r in pool_results:
        for result in pool_results:
            r = result.get()
            have_qual = False
            if os.path.exists(r + '.qual'):
                have_qual = True
                qual_names = dict()
            fxparser = amplishot.parse.fastx.FastxReader(open(r))
            for name, seq, qual in fxparser.parse(callback=amplishot.parse.fastx.greater_than, 
            length=self.config.data['minimum_reconstruction_length']):
                seq_name = '%s_%i %s' % (sampleName, full_length_counter, name)
                self.fullLengthSeqs[seq_name] = seq
                if have_qual:
                    qual_names[name] = seq_name
                logging.debug("Assigning %s from %s", seq_name, r)
                full_length_counter += 1
            
            if have_qual:
                qualparser = amplishot.parse.fastx.QualityReader(open(r + '.qual'))
                for name, qual in qualparser.parse():
                    if name in qual_names:
                        self.fullLengthQuals[qual_names[name]] = qual


    def write(self, fp, qfp = None):
        for name, seq in self.fullLengthSeqs.items():
            fp.write('>%s\n%s\n' %(name, seq))
            if qfp is not None:
                qfp.write('>%s\n%s\n' % (name, self.fullLengthQuals[name]))
