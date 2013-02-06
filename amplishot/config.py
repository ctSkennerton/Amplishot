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
import logging
import sys
import os
import yaml
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper
###############################################################################
###############################################################################
###############################################################################
###############################################################################
_DEFAULT_CONFIG = '''
# Example Amplishot pipeline configuration file
# This is a YAML 1.4 file
#
# The config file is separated into blocks based on indentation
# Each level contains a controlled vocabulary of options for 
# each section of the Amplishot pipeline.
# NOTE: not all steps need to be listed in the config file, if they
# are missing it is assumed that they have been done previously/manually
---
overlap: # assemble the pairs of reads into a consensus 'pairtig'
    #input_files: # path to files from the Illumina machine
    #    - [one.1.fq, one.2.fq] # Each line contains the files for each sample. There should be two files per sample 
    #    - [two.1.fq, two.2.fq]
    minimum_read_length: 350 # minimum length of the overlapped pairs
    overlap_length: 30 # mimimum length of the overlap
    #output_files: # explicitly state the output file names.
    #    - one.olap.fq
    #    - two.olap.fq
map: # map pairtigs against a reference 16S database then seqgregate them based on taxonomy
    #input_files: # fastq files of overlapped reads
    #    - one.olap.fq
    #    - two.olap.fq
    #aliases: # names to give your files in the output.  Must be the same length as 'files'
    #    - sample1
    #    - sample2
    mapper: bowtie # program used for read mapping 
    #database: &BOWTIEDB file.index # path to file of mapper index file
    output_directory_prefix: &GLOBALOUTDIR root # root directory for taxonomic partitioning.  directories will be in the form: $output_directory_prefix/Bacteria/Proteobacteria/...
    output_read_names: &READS reads.fa # name for file containing reads of a taxon. eg. $output_directory_prefix/Bacteria/Proteobacteria/.../$output_read_names
    similarity: 0.98 # the sequence similarity required between the reference database and the reads
    taxon_coverage: [2, 1000] # list of two numbers. The first is the minimum coverage, the second is the number of bases that need to be covered
reduce: # after partitioning, remove excessive coverage by read clustering
    method: cdhit
    input_directory_prefix: *GLOBALOUTDIR
    similarity: 0.98 # sequence similarity between two reads to be clustered together
    cdhit_max_memory: 1000 # maximum memory allowed for cdhit - does not apply to other reduction methods
    input_file_name: *READS
    output_file_name: &CDHITREADS cdhitout.fa
assemble:
    method: phrap # choose a genome assembler  
    input_directory_prefix: *GLOBALOUTDIR
    input_file_name: *CDHITREADS
    phrap_minscore: 300 # minimum alignment score
    assemble_unknowns: false # choose whether to assemble reads that had no match during mapping. (VERY SLOW WITH PHRAP)
repset: # pick representative sequences from all samples and assign taxonomy
    input_directory_prefix: *GLOBALOUTDIR
    output_file: full_length_sequences.fa
    minimum_length: 1000 # minimum length of sequences that we define as 'full length'
    similarity: 0.97 # the similarity used for clustering full-length sequences from different samples into OTUs
    assign_taxonomy_method: bowtie
    #database: *BOWTIEDB
    similarity: 0.90 # sequences that fall below this cutoff will be listed as no taxonomy
    normalize: true # output a normalized OTU table as well as non-normalized
'''
class AmplishotConfig(object):
    """ Class for reading the config file
        The Amplishot config file in written in YAML and contains all of the
        options for running the pipeline
    """
    def __init__(self, config=_DEFAULT_CONFIG):
        """ Load the YAML config file into an object
        """
        self.data = yaml.load(config, Loader=Loader)

    def __str__(self):
        """ dump the YAML to a string
        """
        return yaml.dump(self.data, Dumper=Dumper, default_flow_style=False,
                explicit_start=True)

    def populate_from_commandline(self, args):
        """ take the args that were input on the command line and overwrite the
            defaults
        """
        for key, value in vars(args).items():
            (pipeline_stage, key) = key.split('_', 1)
            self.data[pipeline_stage][key] = value
###############################################################################
###############################################################################
###############################################################################
###############################################################################
