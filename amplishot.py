#!/usr/bin/env python
###############################################################################
#
# amplishot.py - Pipeline for generating full-length 16S sequences and performing
#                community abundance measurements on the output#
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
__copyright__ = "Copyright 2012-2013"
__credits__ = ["Connor Skennerton"]
__license__ = "GPL3"
__version__ = "0.0.1"
__maintainer__ = "Connor Skennerton"
__email__ = "c.skennerton@gmail.com"
__status__ = "Development"

###############################################################################

import argparse
import sys
from pandaseq import Pandaseq
from cogent import DNA, LoadSeqs
import cogent.core.moltype
from cogent.parse.fastq import MinimalFastqParser

#import os
#import errno

#import numpy as np
#np.seterr(all='raise')     

#import matplotlib as mpl
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import axes3d, Axes3D
#from pylab import plot,subplot,axis,stem,show,figure


###############################################################################
###############################################################################
###############################################################################
###############################################################################

###############################################################################
###############################################################################
###############################################################################
###############################################################################
def panda(args):
    panda = Pandaseq()
    panda.Parameters['-o'].on(args.olap)
    panda.Parameters['-l'].on(args.minlen)
    panda.Parameters['-f'].on(args.read_1)  
    panda.Parameters['-r'].on(args.read_2)
    return panda()

def disambiguate(sequence):
    index = sequence.firstDegenerate()
    if index is None:
        return sequence

    head = sequence[:index]
    tail = sequence[index+1:]
    ret_seqs = list()

    for unambiguous in cogent.core.moltype.IUPAC_DNA_ambiguities[sequence[index]]:
        ret = disambiguate(head+unambiguous+tail)
        if isinstance(ret, list):
            ret_seqs.extend(ret)
        else:
            ret_seqs.append(ret)
    return ret_seqs

def generate_unambiguous_sequences():
    _CONSERVED_SEQUENCES = {
        '357': 'CTGAGAYACGGHCCARACTCCTACGGGAGGCAGCAG',
        '530': 'GGCTAAYTHYGTGCCAGCAGCCGCGGTAAKAC',
        '787': 'CRAACRGGATTAGATACCCYGGTAGTCCW',
        '926': 'AAACTYAAAKGAATTGRCGG',
        '1114': 'GGTGSTGCATGGYTGTCGTCAGCTCGTGYCGTGA',
        '1392': 'GGTGAATACGTTCCCGGGYCTTGYACNCAC'
        }

    unambiguous_conserved_sequences = list()
    rev_unambiguous_conserved_sequences = list()
    for seq in _CONSERVED_SEQUENCES.values():
        dnaseq = DNA.makeSequence(seq)
        #print dnaseq, dnaseq.__class__
        ret = disambiguate(dnaseq)
        if isinstance(ret, list):
            unambiguous_conserved_sequences.extend(ret)
        else:
            unambiguous_conserved_sequences.append(ret)

    for i in unambiguous_conserved_sequences:
        rev_unambiguous_conserved_sequences.append(i.rc())

    unambiguous_conserved_sequences.extend(rev_unambiguous_conserved_sequences)
    
    return unambiguous_conserved_sequences

def doWork( args ):
    """ Main wrapper"""

    # generate fragment consensus sequences - pandaseq
    #panda_out = panda(args)

    # dataset partitioning based on conserved primer presence
    seqs = generate_unambiguous_sequences()
    for i in seqs:
        tree.add(str(i))
    tree.make()
    count = 0
    with open(args.infile) as fp:
        for name, seq, qual in MinimalFastqParser(fp):
            pass
    print count
    # dataset reduction - cd-hit - each partition - threaded

    # generate overlaps - phrap - each partition - threaded

    #combined assemblies - phrap

    #chimera checking - uchime, decipher, chimeraslayer

    #read alignment - bwa

    #generate abundance based on coverage

    #assign taxonomy

    #generte OTU table

    return 0

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    #parser.add_argument('positional_arg', help="Required")
    #parser.add_argument('positional_arg2', type=int, help="Integer argument")
    #parser.add_argument('positional_arg3', nargs='+', help="Multiple values")
    parser.add_argument('-i', '--infile', dest='infile', 
            help="input file containing overlapped reads from pandaseq")
    
    # parse the arguments
    args = parser.parse_args()        

    # do what we came here to do
    doWork(args)

###############################################################################
###############################################################################
###############################################################################
###############################################################################
