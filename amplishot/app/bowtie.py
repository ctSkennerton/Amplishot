#!/usr/bin/env python
###############################################################################
#
# bowtie.py - app controller for bowtie
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

from Bio.Application import _Switch, _Option, AbstractCommandline


class Bowtie(AbstractCommandline):
    """ Simple bowtie2 application controller
    """
    def __init__(self, cmd='bowtie2', **kwargs):

        self.parameters = [
            _Option(['-x','bt2-idx','x'],
            '''Index filename prefix (minus trailing .X.bt2).
         NOTE: Bowtie 1 and Bowtie 2 indexes are not compatible.''',
         filename=True, is_required=True)
            _Option(['-1','1','m1'],
            '''Files with #1 mates, paired with files in <m2>.
         Could be gzip'ed (extension: .gz) or bzip2'ed (extension:
         .bz2).''', filename=True)
            _Option(['-2','2','m2'],
            '''Files with #2 mates, paired with files in <m1>.
         Could be gzip'ed (extension: .gz) or bzip2'ed (extension:
         .bz2).''', filename=True)
            _Option(['-U','U'],
            '''Files with unpaired reads.
         Could be gzip'ed (extension: .gz) or bzip2'ed (extension:
         .bz2).''', filename=True)
            _Option(['-S','sam',
                '''File for SAM output (default: stdout)''',
                filename=True)

            _Switch(['-q','q'],
                '''query input files are FASTQ .fq/.fastq (default)''')
            _Switch(['--qseq','qseq'],
                '''query input files are in Illumina's qseq format''')
            _Switch(['-f','f'],
                '''query input files are (multi-)FASTA .fa/.mfa''')
            _Switch(['-r','r'],
                '''query input files are raw one-sequence-per-line''')
            _Switch(['-c','c'],
                '''<m1>, <m2>, <r> are sequences themselves, not files''')
            _Option(['-s','s','skip'],
                '''skip the first <int> reads/pairs in the input (none)''')
            _Option(['-u','u','upto'],
                '''stop after first <int> reads/pairs (no limit)''')
            _Option(['-5','5','trim5'],
                '''trim <int> bases from 5'/left end of reads (0)''')
            _Option(['-3','3','trim3'],
                '''trim <int> bases from 3'/right end of reads (0)''')
            _Switch(['--phred33', 'phred33'],
                '''qualities are Phred+33 (default)''')
            _Switch(['--phred64', 'phred64'],
                '''qualities are Phred+64''')
            _Switch(['--int-quals', 'int-quals'],
                '''qualities encoded as space-delimited integers''')

            _Switch(['--very-fast', 'very-fast'],
                '''-D 5 -R 1 -N 0 -L 22 -i S,0,2.50''')
            _Switch(['--fast', 'fast'],
                '''-D 10 -R 2 -N 0 -L 22 -i S,0,2.50''')
            _Switch(['--sensitive', 'sensitive'],
                '''-D 15 -R 2 -N 0 -L 22 -i S,1,1.15 (default)''')
            _Switch(['--very-sensitive', 'very-sensitive'],
                '''-D 20 -R 3 -N 0 -L 20 -i S,1,0.50''')

            _Switch(['--very-fast-local', 'very-fast-local',
                '''-D 5 -R 1 -N 0 -L 25 -i S,1,2.00''')
            _Switch(['--fast-local', 'fast-local'],
                '''-D 10 -R 2 -N 0 -L 22 -i S,1,1.75''')
            _Switch(['--sensitive-local', 'sensitive-local'],
                '''-D 15 -R 2 -N 0 -L 20 -i S,1,0.75 (default)''')
            _switch(['--very-sensitive-local', 'very-sensitive-local' ],
                '''-D 20 -R 3 -N 0 -L 20 -i S,1,0.50''')

            _Option(['-N', 'N'],
                '''max # mismatches in seed alignment; can be 0 or 1 (0)''')
            _Option(['-L', 'L'],
                '''length of seed substrings; must be >3, <32 (22)''')
            _Option(['-i', 'i'],
                '''interval between seed substrings w/r/t read len (S,1,1.15)''')
            _Option(['--n-ceil', 'n-ceil'],
                '''func for max # non-A/C/G/Ts permitted in aln (L,0,0.15)''')
            _Option(['--dpad', 'dpad'],
                '''extra ref chars on sides of DP table (15)''')
            _Option(['--gbar', 'gbar'],
                '''nucs of read extremes (4)''')
            _Switch(['--ignore-quals', 'ignore-quals'],
                '''treat all quality values as 30 on Phred scale (off)''')
            _Switch(['--nofw', 'nofw'],
                '''do not align forward (original) version of read (off)''')
            _Switch(['--norc', 'norc'],
                '''do not align reverse-complement version of read (off)''')

            _Switch(['--end-to-end', 'end-to-end'],
                '''entire read must align; no clipping (on)''')
            _Switch(['--local', 'local'],
                '''local alignment; ends might be soft clipped (off)''')

            _Option(['--ma', 'ma'],
                '''match bonus (0 for --end-to-end, 2 for --local) ''')
            _Option(['--mp', 'mp'],
                '''max penalty for mismatch; lower qual = lower penalty (6)''')
            _Option(['--np', 'np'],
                '''penalty for non-A/C/G/Ts in read/ref (1)''')
            _Option(['--rdg', 'rdg'],
                '''read gap open, extend penalties (5,3)''')
            _Option(['--rfg', 'rfg'],
                '''reference gap open, extend penalties (5,3)''')
            _Option(['--score-min', 'score-min'],
                '''min acceptable alignment score w/r/t
                read length (G,20,8 for local, L,-0.6,-0.6 for end-to-end)''')

            _Option(['-k', 'k'],
                '''report up to <int> alns per read; MAPQ not meaningful''')
            _Switch['-a', 'a', 'all'],
            '''report all alignments; very slow, MAPQ not meaningful''')
            _Option(['-D', 'D'],
                    '''give up extending after <int> failed extends in a row (15)''')
            _Option(['-R', 'R'],
                    '''for reads w/ repetitive seeds, try <int> sets of seeds (2)''')

            _Option(['-I','I','minins'],
                    '''minimum fragment length (0)''')
            _Option(['-X','X','maxins'],
                    '''maximum fragment length (500)''')
            _Switch(['--fr', 'fr', 'forward-reverse'],
                    '''mates are in forward/reverse orientation ''')
            _Switch(['--rf', 'rf', 'reverse-forward'],
                    ''' mates are in revese/forward orientation''')
            _Switch(['--ff', 'ff', 'forward-forward'],
                    '''mates are in forward/forward orientation'''
            _Switch(['--no-mixed', 'no-mixed'],
                '''suppress unpaired alignments for paired reads''')
            _Switch(['--no-discordant', 'no-discordant'],
                '''suppress discordant alignments for paired reads''')
            _Switch(['--no-dovetail', 'no-dovetail'],
                '''not concordant when mates extend past each other''')
            _Switch(['--no-contain', 'no-contain'],
                '''not concordant when one mate alignment contains other''')
            _Switch(['--no-overlap', 'no-overlap'],
                '''not concordant when mates overlap at all''')

            _Switch(['-t', 't', 'time'],
                '''print wall-clock time taken by search phases''')
            _Option(['--un', 'un'],
                '''write unpaired reads that didn't align to <path>''')
            _Option(['--al', 'al'],
                '''write unpaired reads that aligned at least once to <path>''')
            _Option(['--un-conc', 'un-conc'],
                '''write pairs that didn't align concordantly to <path>''') 
            _Option(['--al-conc', 'al-conc'],
                '''write pairs that aligned concordantly at least once to <path>''') 
            _Option(['--un-gz', 'un-gz'],
                '''to gzip compress output, or add '-bz2' to bzip2 compress output.)''') 
            _Switch(['--quiet', 'quiet'],
                '''print nothing to stderr except serious errors''')
            _Option(['--met-file', 'met-file'],
                '''send metrics to file at <path> (off)''')
            _Switch(['--met-stderr', 'met-stderr'],
                '''send metrics to stderr (off)''') 
            _Option(['--met', 'met'], '''secs (1)''')
            _Switch(['--no-head', 'no-head'],
                '''supppress header lines, i.e. lines starting with @''')
            _Switch(['--no-sq', 'no-sq'],
                '''supppress @SQ header lines''')
            _Option(['--rg-id', 'rg-id'],
                '''set read group id, reflected in @RG line and RG:Z: opt field''')
            _Option(['--rg', 'rg'],
                '''("lab:value") to @RG line of SAM header.''')
            _Switch(['--omit-sec-seq','omit-sec-seq'],
                '''put '*' in SEQ and QUAL fields for secondary alignments.''')
            _Option(['-o','o','offrate'],
                '''override offrate of index; must be >= index's offrate''')
            _Option(['-p','p','threads'],
                '''number of alignment threads to launch (1)''')
            _Switch(['--reorder', 'reorder'],
                '''force SAM output order to match order of input reads''')
            _Switch(['--mm', 'mn'],
                '''use memory-mapped I/O for index; many 'bowtie's can share''')

            _Switch(['--qc-filter', 'qc-filter'],
                '''filter out reads that are bad according to QSEQ filter''')
            _Option(['--seed', 'seed'],
                '''seed for random number generator (0)''')
            _Switch(['--non-deterministic', 'non-deterministic'], 
                '''seed rand. gen. arbitrarily instead of using read attributes''')
            _Switch(['--version', 'version'], 
                '''print version information and quit''')
            _Switch(['-h', 'h','help'], '''print this usage message''')
        ]
