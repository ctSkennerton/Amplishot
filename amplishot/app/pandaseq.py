#!/usr/bin/env python
###############################################################################
#
# pandaseq.py - application controller for pandaseq
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
__version__ = "0.1.0"
__maintainer__ = "Connor Skennerton"
__email__ = "c.skennerton@gmail.com"
__status__ = "Development"

###############################################################################
from Bio.Application import _Switch, _Option, AbstractCommandline

class Pandaseq(AbstractCommandline):
    """Simple pandaseq application controller.
    """
    def __init__(self, cmd='pandaseq', **kwargs):
        self.parameters = [
            _Switch(['-6', '6', 'six'],
                'Use PHRED+64 (CASAVA 1.3-1.7) instead of PHRED+33 (CASAVA 1.8+)'),
            _Switch(['-B', 'B', 'unbarcoded'],
                'Allow unbarcoded sequences (try this for BADID errors)'),
            _Option(['-C', 'C', 'module'],
                'Load a sequence validation module', equate=False),
            _Option(['-d', 'd', 'logging_flags'],
            'Control the logging messages. Capital to enable; small to disable.\
    		(R)econstruction detail.\
    		Sequence (b)uilding information.\
    		(F)ile processing.\
    		(k)-mer table construction.\
    		Show every (m)ismatch.\
    		Optional (s)tatistics.', equate=False),
            _Option(['-f', 'f', 'forward_reads'],'Input FASTQ file containing forward reads.',
                 filename=True, is_required=True, equate=False),
            _Switch(['-F', 'F', 'fastq'], 
                'Output FASTQ instead of FASTA.'),
            _Switch(['-j', 'j', 'bzipped'],
                'Input files are bzipped.'),
            _Option(['-l', 'l', 'minlength'],
                'Minimum length for a sequence', equate=False),
            _Option(['-L', 'L', 'maxlength'],
                'Maximum length for a sequence', equate=False),        
            _Switch(['-N', 'N', 'no_unknown'],
                'Eliminate all sequences with unknown nucleotides in the output.'),
            _Option(['-o', 'o', 'overlap'],
                'Minimum overlap between forward and reverse reads (default = 1)',
                equate=False),
            _Option(['-p', 'p', 'forward_primer'], 
                'Forward primer sequence or number of bases to be removed.',
                equate=False),
            _Option(['-q', 'q', 'reverse_primer'],
                'Reverse primer sequence or number of bases to be removed.',
                equate=False),
            _Option(['-r', 'r', 'reverse_reads'],
                'Input FASTQ file containing reverse reads.', 
                filename=True, is_required=True, equate=False),
            _Option(['-t', 't', 'primer_probability'],
                'The minimum probability that a sequence must have to match a primer. (default = 6.000000e-01)'),
            _Option(['-T', 'T', 'threads'],
                'Run with a number of parallel threads.', equate=False),
            ]
        super(Pandaseq, self).__init__(self, cmd, **kwargs)
        

