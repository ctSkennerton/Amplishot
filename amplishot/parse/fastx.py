#!/usr/bin/env python
###############################################################################
#
# fastx.py - take a bam file and a taxonomy and segregate the reads
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
__version__ = "0.3.2"
__maintainer__ = "Connor Skennerton"
__email__ = "c.skennerton@gmail.com"
__status__ = "Development"

###############################################################################
import collections

def decode_quality(qual, offset=33):
    qual = map(lambda x: ord(x) - offset, qual)
    return qual

def calcGc(seq):
    seq2 = seq.upper()
    d = collections.defaultdict(int)
    for c in seq2:
        d[c] += 1
    return float(d['G'] + d['C'])/float(len(seq))

def greater_than(name, seq, qual, length=None, gc=None, inclusive=True):
    """ Returns a list of sequences that are >(=) to the parameters set 
    """
    if inclusive:
        if length and not len(seq) >= length:
            return False
        if gc and not calcGc(seq) >= gc:
            return False
    else:
        if length and not len(seq) > length:
            return False
        if gc and not calcGc(seq) > gc:
            return False
    return True

def less_than(name, seq, qual, length=None, gc=None, inclusive=True):
    """ Returns a list of sequences that are <(=) to the parameters set 
    """
    if inclusive:
        if length and not len(seq) <= length:
            return False
        if gc and not calcGc(seq) <= gc:
            return False
    else:
        if length and not len(seq) < length:
            return False
        if gc and not calcGc(seq) < gc:
            return False
    return True

def include(name, seq, qual, headers=None):
    """Returns sequences that contain the headers
    """
    return name in headers

def not_include(name, seq, qual, headers=None):
    """Returns sequences that contain the headers
    """
    return name not in headers

class FastxReader(object):
    """ Class for extracting sequences from a fastx file 
    """
    def __init__(self, fp):
        super(FastxReader, self).__init__()
        self.fp = fp
        if isinstance(fp, str):
            self.fp = open(self.fp)


    def parse(self, callback=None, **kwargs):
        last = None # this is a buffer keeping the last unprocessed line
        while True: # mimic closure; is it a bad idea?
            if not last: # the first record or a record following a fastq
                for l in self.fp: # search for the start of the next record
                    if l[0] in '>@': # fasta/q header line
                        last = l[:-1] # save this line
                        break
            if not last: break
            name, seqs, last = last[1:].split()[0], [], None
            for l in self.fp: # read the sequence
                if l[0] in '@+>':
                    last = l[:-1]
                    break
                seqs.append(l[:-1])
            if not last or last[0] != '+': # this is a fasta record
                if callback:
                    if callback(name, ''.join(seqs), None, **kwargs):
                        yield name, ''.join(seqs), None # yield a fasta record
                else:
                    yield name, ''.join(seqs), None # yield a fasta record
                if not last: break
            else: # this is a fastq record
                seq, leng, seqs = ''.join(seqs), 0, []
                for l in self.fp: # read the quality
                    seqs.append(l[:-1])
                    leng += len(l) - 1
                    if leng >= len(seq): # have read enough quality
                        last = None
                        if callback:
                            if callback(name, seq, ''.join(seqs), **kwargs):
                                yield name, seq, ''.join(seqs); # yield a fastq record
                        else:
                            yield name, seq, ''.join(seqs); # yield a fastq record
                        break
                if last: # reach EOF before reading enough quality
                    if callback:
                        if callback(name, seq, None, **kwargs):
                            yield name, seq, None # yield a fasta record instead
                    else:
                        yield name, seq, None # yield a fasta record instead
                    break

class QualityReader(object):
    """ Class for extracting sequences from a fastx file 
    """
    def __init__(self, fp):
        super(QualityReader, self).__init__()
        self.fp = fp
        if isinstance(fp, str):
            self.fp = open(self.fp)

    def parse(self, callback = None, **kwargs):
        last = None # this is a buffer keeping the last unprocessed line
        while True: # mimic closure; is it a bad idea?
            if not last: # the first record or a record following a fastq
                for l in self.fp: # search for the start of the next record
                    if l[0] in '>@': # fasta/q header line
                        last = l[:-1] # save this line
                        break
            if not last: break
            name, seqs, last = last[1:].split()[0], [], None
            for l in self.fp: # read the sequence
                if l[0] in '@+>':
                    last = l[:-1]
                    break
                seqs.append(l[:-1])
            if not last or last[0] != '+': # this is a fasta record
                if callback:
                    if callback(name, ' '.join(seqs), **kwargs):
                        yield name, ' '.join(seqs) # yield a fasta record
                else:
                    yield name, ' '.join(seqs) # yield a fasta record
                if not last: break
            else: # this is a fastq record
                raise RuntimeError, 'This appears to be a fastq record, not\
                fasta quality'
