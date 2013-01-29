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
__version__ = "0.0.1"
__maintainer__ = "Connor Skennerton"
__email__ = "c.skennerton@gmail.com"
__status__ = "Development"

###############################################################################

import tempfile
import os
import amplishot.parse.sam
import amplishot.app.bowtie
from qiime.pick_otus  import otu_picking_method_constructors,\
 otu_picking_method_choices, MothurOtuPicker
from qiime.pick_rep_set import (rep_set_picking_methods,
 reference_rep_set_picking_methods)
import biom.table
import numpy as np
###############################################################################
###############################################################################
###############################################################################
###############################################################################


class OTUTableGenerator(object):
    """Get the abundance of each of the representative sequences from our OTUs
    """
    def __init__(self, fullLengthSeqs, workingDir,
            outFilePrefix='amplishot_full_length'):
        """ Save the names and lengths of the full-length sequences for the
        mapping with the sam files.  While we are here we should also index the
        sequences for mapping the reads
        fullLengthSeqs: Path to file containing full length sequences
        """
        self.seqs = dict()
        self.aliases = list()
        self.obsrevation_ids = list()
        self.outdir = workingDir
        self.outprefix = outFilePrefix
        self.biom_table = None
        counter = 0
        for name, seq in fullLengthSeqs.items():
            self.seqs[name] = (len(seq), counter)
            self.observation_ids.append(name)
            counter += 1

        # initialize the first column of the numpy array
        self.fpkc = np.zeros(len(self.seqs))

        self._make_bowtie_index(fullLengthSeqs)

    def __str__(self):
        version_str = 'amplishot %s' % __version__
        return self.biom_table.getBiomFormatJsonString(generated_by=version_str)
    
    def _make_bowtie_index(self, fullLengthSeqs):
        """ write the full-length sequences to file then index them with
        bowtie.  After which the tmpfile can be deleted.
        """
        tmp = tempfile.NamedTemporaryFile(delete=False)
        for name, seq in fullLengthSeqs.items():
            tmp.write('>%s\n%s\n' % name, seq)

        tmp.close()
        bi = amplishot.app.bowtie.BowtieIndex(infile=tmp.name,
                prefix=self.outprefix)
        bi(stdout=False, stderr=False, cwd=self.outdir)

        os.remove(tmp.name)

    def generate_abundance(self, reads, alias=None):
        """ Take in a set of reads and map them with bowtie to the rep set
        reads: should be a list of paths to files containing reads
        Would have been one of the same ones that was originally
        used for taxon segregation
        alias: name to give this file in the OTU table,  default is to take the
        basename of the reads file minus the extension
        """
        stdout, stderr = self._make_sam(reads)
        tmp_fpkc = self._parse_sam(stdout)
        if alias is None:
            alias = os.path.splitext(os.path.basename(reads))[0]
        self.aliases.append(alias)

        # check to see if we need to append on a new column to the 2D-array
        if len(self.aliases) > self.fpkc.shape[1]:
            self.fpkc = np.append(self.fpkc, tmp_fpkc, 1)

    def normalize(self):
        pass

    def generate_biom_table(self, sample_metadata=None, observation_metadata=None):
        self.biom_table = biom.table.table_factory(self.fpkc, self.aliases,
                self.observation_ids, sample_metadata=sample_metadata,
                observation_metadata=observation_metadata)

    def add_metadata(self, sample_metadata=None, observation_metadata=None):
        if sample_metadata is not None:
            self.biom_table.addSampleMetadata(sample_metadata)
        if observation_metadata is not None:
            self.biom_table.addObservationMetadata(observation_metadata)
    
    def _make_sam(self, reads):
        """ Call bowtie on a combined set of full length sequences
        """
        b = amplishot.app.bowtie.Bowtie(index=os.path.join(self.outdir,
            self.outprefix), unpaired_reads=reads)
        return b(seterr=False)

    def _parse_sam(self, sam, percentId=0.97):
        """ count the number of reads that map to a full length sequence
        taking into account for dud reads that need to be filtered
        """
        tmp_counts = dict()
        samfile = amplishot.parse.sam.SamFilereader(sam, parseHeader=False)
        for read in samfile.parse():
            tmp_counts[read.rname] += 1

        return self._get_fpkc(tmp_counts)


    def _get_fpkc(self, counts):
        """ Calculate the fragments mapped per kilobase contig length (fpkc),
        which will be used as our abundance measurement
        """
        tmp = [0.0] * len(self.seqs)
        for name in counts.keys():
            fpkc = float(counts[name] * 1000) / float(self.seqs[name])
            tmp[self.seqs[name][1]] = fpkc
        return tmp

    def _normalize_fpkc(self):
        """ Scale the fpkc values based on the total number of reads that
        mapped in the sam files
        """
        pass


def pick_otus(inputSeqsFilepath, otuPickingMethod='cdhit', similarity=0.98,
        maxCdhitMemory=1000, outputFileName='full_length_otus.txt'):
    """ Call out to qiime to pick the OTUs
    Qiime has a number of OTU picking methods, which I intend to support at
    some point.  However at the moment I think that I'll be opinionated and say
    that you have to use CD-HIT
    """
    # code copied from pick_otus.py in the qiime scripts dicectory
    otu_picker_constructor =\
     otu_picking_method_constructors[otuPickingMethod]
    params = {'Similarity': similarity,
            '-M': maxCdhitMemory}
    otu_picker = otu_picker_constructor(params)

    otu_picker(inputSeqsFilepath, result_path=outputFileName)
    return outputFileName


def pick_rep_set(inputSeqsFilpath, inputOtuMapFilePath,
        outputFilePath='full_length_rep_set.fa', pickingMethod='longest',
        sortBy='otu'):
    """ Call out to qiime for picking the representative set
    """
    # code copied from pick_rep_set.py in qiime scripts diectory

    #if reference_seqs_filepath:
    #    rep_set_picker =\
    #            reference_rep_set_picking_methods[opts.rep_set_picking_method]
    #    rep_set_picker(input_seqs_filepath,
                #input_otu_filepath,
                #reference_seqs_filepath,
                #result_path=result_path,
                #log_path=log_path,
                #sort_by=opts.sort_by)
    #else:
    rep_set_picker = rep_set_picking_methods[pickingMethod]
    rep_set_picker(inputSeqsFilepath,
        inputOtuMapFilePath,
        result_path=outputFilePath,
        sort_by=sortBy)
    return outputFilePath