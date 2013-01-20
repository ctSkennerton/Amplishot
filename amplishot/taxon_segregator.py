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

import sys
#import pysam
import os
import re
import subprocess
import tempfile
###############################################################################
###############################################################################
###############################################################################
###############################################################################
DNA_COMPLEMENT_TABLE = {'A':'T',
                        'T':'A',
                        'C':'G',
                        'G':'C',
                        'a':'t',
                        't':'a',
                        'c':'g',
                        'g':'c',
                        'N':'N',
                        'n':'n'
                        }
def reverse_complement(seq):
    seq = seq[::-1]
    ret = str()
    for char in seq:
        ret += DNA_COMPLEMENT_TABLE[char]
    return ret


def decode_quality(qual, offset=33):
    qual = map(lambda x: ord(x) - offset, qual)
    return qual

class SamRead(object):
    
    def __init__(self, fields):
        super(SamRead, self).__init__()
        self.qname = fields[0]
        self.flag = int(fields[1])
        self.rname = fields[2]
        self.pos = int(fields[3])
        self.mapq = int(fields[4])
        self.cigar = fields[5]
        self.rnext = fields[6]
        self.pnext = int(fields[7])
        self.tlen = int(fields[8])
        self.seq = fields[9]
        self.qlen = len(self.seq)
        self.qual = fields[10]
        self.tags = {}
        tags_split = fields[11].split('\t')
        for tag in tags_split:
            tag_name, value_type, value = tag.split(":")
            if value_type == 'i':
                value = int(value)
            self.tags[tag_name] = value
    
    def has_multiple_segments(self):
        return self.flag & 0x1

    def is_properly_aligned(self):
        return self.flag & 0x2

    def is_unmapped(self):
        return self.flag & 0x4

    def is_next_segment_unmapped(self):
        return self.flag & 0x8

    def is_reversed(self):
        return self.flag & 0x10

    def is_next_segment_reversed(self):
        return self.flag & 0x20

    def is_first_segment(self):
        return self.flag & 0x40

    def is_last_segment(self):
        return self.flag & 0x80

    def is_secondary_alignment(self):
        return self.flag & 0x100

    def is_qc_fail(self):
        return self.flag & 0x200

    def is_duplicate(self):
        return self.flag & 0x400

class TaxonSegregator(object):
    """Take a taxonomy and a samfile and segregate the reads into their taxonomic
       divisions
    """

    # greengenes separate file format regex
    ggs_re = re.compile('\w+\t([kpcofg]__.* ){6}s__.*')
    # greengenes combined file format regex
    ggc_re = re.compile('\w+\t([kpcofg]__.*){6}s__.*')
    # silva file format regex
    sil_re = re.compile('.+\t(\w+;)*\w+')
    def __init__(self, taxonfile):
        """
           taxonfile: path to a file containing the mapping between reference
                      identifiers and the taxonomy string
        """
        super(TaxonSegregator, self).__init__()
        self.done_segregation = False

        self.taxon_mapping = dict()
        self.ref_taxon_mapping = dict()
        self.taxon_format = ''

        # this dict will hold the samfile header information for each of the 
        # segregated taxonomies.  The header dictionary should be constructed
        # as the following shows:
        #header = { 'HD': {'VN': '1.0'},
        #    'SQ': [{'LN': 1575, 'SN': 'chr1'},
        #           {'LN': 1584, 'SN': 'chr2'}] }
        self.taxon_header = dict()
        self.taxon_references = set()

        with open(taxonfile) as fp:
            self._parse_taxon_file(fp)

    def _check_taxonomy_file_format(self, taxonfp):
        """check which taxonomy format is being used
        """
        start = taxonfp.tell()
        firstline = taxonfp.readline()
        taxonfp.seek(start)
        if self.ggs_re.match(firstline):
            return 'ggs'
        elif self.ggc_re.match(firstline):
            return 'ggc'
        elif self.sil_re.match(firstline):
            return 'sil'
        else:
            raise RuntimeError, "Taxon file format does not match any known"

    def _parse_greengenes_sparate(self, taxonfp):
        """greengenes taxonomy file where each taxonomic division is space 
           separated
           12345    k__x; p__x; c__x; o__x; f__x; g__x; s__x
        """
        for line in taxonfp:
            line = line.rstrip()
            (refid, taxon_string) = line.split('\t')
            taxon_divisions = taxon_string.split('; ')
            names = []
            for rank in taxon_divisions:
                names.append(rank[3:])

            self._generate_mapping(refid, names)

    def _parse_greengenes_combined(self, taxonfp):
        """greengenes taxonomy file where the taxonomy string has no spaces
           12345    k__x;p__x;c__x;o__x;f__x;g__x;s__x
        """
        for line in taxonfp:
            line = line.rstrip()
            (refid, taxon_string) = line.split('\t')
            taxon_divisions = taxon_string.split(';')
            names = []
            for rank in taxon_divisions:
                names.append(rank[3:])

            self._generate_mapping(refid, names)

    def _parse_silva(self, taxonfp):
        """silva formated taxonomy file - no taxon level designation
           12345    a;b;c;d;e;f;g
        """
        for line in taxonfp:
            line = line.rstrip()
            (refid, taxon_string) = line.split('\t')
            taxon_divisions = taxon_string.split(';')
            self._generate_mapping(refid, taxon_divisions)

    def _generate_mapping(self, refid, taxon_divisions):
        """take the input from one of the parsers and generate a mapping 
           of tuple of taxon divisions to a reference id a the opposite
           reference id to tuple
        """
        taxon_divisions = filter(None, taxon_divisions)
        t = tuple(taxon_divisions)
        self.taxon_mapping[t] = []
        self.ref_taxon_mapping[refid] = t        

    def _parse_taxon_file(self, taxonfp):
        """firse check the file format then call the correct file parser
        """
        format = self._check_taxonomy_file_format(taxonfp)
        if format == 'ggs':
            self._parse_greengenes_sparate(taxonfp)
        elif format == 'ggc':
            self._parse_greengenes_combined(taxonfp)
        elif format == 'sil':
            self._parse_silva(taxonfp)

    def _check_taxon_coverage(self, taxon, minCount=1000, minCoverage=2):
        """calculate the per base coverage for this taxon
           If a taxon does not contain sufficient coverage across the length of
           the reference sequences then it should be merged into a higher taxonomy
           Since the reference sequences should be of approximate length, there 
           should not be too much problem when dealing with multiple references
           from the same taxon.  However it is likely that this coverage information
           will be a little bit 'fuzzy' since there is some variation

           taxon: A tuple containing the taxon string
           minCount: the minimum number of positions that must have coverage
           minCoverage: the minimum coverage allowed for the covered positions
        """
        #assume that the sequences are not greater than 1600bp
        coverage = [0]*1600
        
        for read in self.taxon_mapping[taxon]:
            for i in range(read.pos, read.pos + read.qlen):
                coverage[i] += 1

        count = 0
        for pos in coverage:
            if pos > minCoverage:
                count += 1

        return count >= minCount
    
    def _sam_to_fastx(self, alignedRead, fasta=None, qual=None, fastq=None):
        """Take a pysam AlignedRead and convert it into a fasta, qual or fastq

           alignedRead: A SamRead object
           fasta:       python file object for the sequence in fasta format.  When set
                        to None, no output will be given
           qual:        python file object for the qualaty scores.  Wneh set to None, no  
                        output will be given
           fastq:       python file object for fastq output.  When set to None, no output
                        will be given
        """
        name = alignedRead.qname
        seq = alignedRead.seq
        quality = alignedRead.qual
        if alignedRead.is_reversed():
            seq = reverse_complement(alignedRead.seq)
            quality = quality[::-1]

        if fasta is not None:
            fasta.write('>%s\n%s\n' %(name, seq))

        if fastq is not None:
            fastq.write('@%s\n%s\n+\n%s' % (name,seq,quality))
            
        if qual is not None:
            quality = decode_quality(quality)
            quality = ' '.join(str(x) for x in quality)
            qual.write('>%s\n%s\n' %(name, quality))

    def _open_sam(self, samfile):
        """create an pysam 'samfile' object to parse

           samfile: name of a file in either sam or bam format
        """
        pass
#        return pysam.Samfile(samfile)

    def parse_sam(self, sam, percentId = 0.98):
        """iterate through the records in a samfile and place them into 
           one of the taxonomies based on the mapping

           sam: an opened samfile generated with pysam 
           percentId: The percent identity between the read and the reference
        """
        if self.done_segregation:
            raise RuntimeError, 'Segregation has already taken place.  Please\
            parse all samfiles at one and then call segregate at the end'
        alignments = tempfile.TemporaryFile()
        subprocess.call(['samtools','view', sam], stdout=alignments,
                stderr=open(os.devnull, 'w'))
        alignments.seek(1)

        for line in alignments:
            line = line.rstrip()
            fields = line.split('\t', 11)
            read = SamRead(fields)
            
            # check if the read is mapped
            if not read.is_unmapped():
                if float(read.qlen - read.tags['NM'])/ float(read.qlen) >= percentId:
                    t = self.ref_taxon_mapping[read.rname]
                    self.taxon_mapping[t].append(read)
        #if not isinstance(sam, pysam.Samfile):
        #    sam = self._open_sam(sam)
        #callback = ReadBinner(self.taxon_mapping, self.ref_taxon_mapping,
        #        sam.references)
        #sam.fetch(callback=callback)
        #for read in sam.fetch():
        #    if not count % 1000:
        #        print count
        #    if not read.is_unmapped:
        #        edit_dist = 0
        #        for tag in read.tags:
        #            if tag[0] =='NM':
        #                edit_dist = tag[1]
        #                break
        #        if float(qlen - edit_dist)/ float(qlen) >= percentId:
        #            t = self.ref_taxon_mapping[sam.getrname(read.tid)]
        #            self.taxon_mapping[t].append(read)
        #    count += 1
              #  try:
              #      self.taxon_header[t]['SQ'] = list()
              #  except KeyError:
              #      self.taxon_header[t] = dict()
              #      self.taxon_header[t]['SQ'] = list()

              #  if read.tid not in self.taxon_references:
              #      self.taxon_header[t]['SQ'].append(
              #          { 'LN': sam.header['SQ'][read.tid]['LN'], 
              #          'SN': sam.getrname(read.tid)
              #          })
              #      self.taxon_header[t]['HD'] = {'VN': sam.header['HD']['VN']}
              #      self.taxon_references.add(read.tid)
              #  
        #sam.close()


    def segregate(self, root='root'):
        """Partition all reads into separate files
           Ideally this should be called after all sam files have been parsed
           This function will make a directory structure equal to the taxonomy
           strings for each of the reference sequences.  It will overwrite any
           files already in the directory structure.

           root: The root directory name.  By default creates a directory called 
                 'root' in the current directory
        """
        sorted_taxons = self.taxon_mapping.keys()
        sorted_taxons = sorted(sorted_taxons, reverse=True, cmp=lambda x,y: cmp(len(x), len(y)))
        for taxon_ranks in sorted_taxons:
            if not self._check_taxon_coverage(taxon_ranks) and len(taxon_ranks) > 5:
                print "merge up", taxon_ranks
                new_taxon = []
                for i in range(len(taxon_ranks) - 1):
                    new_taxon.append(taxon_ranks[i])
                
                t = tuple(new_taxon)
                try:
                    self.taxon_mapping[t].extend(self.taxon_mapping[taxon_ranks])
                except KeyError:
                    self.taxon_mapping[t] = self.taxon_mapping[taxon_ranks]

                del self.taxon_mapping[taxon_ranks]
                continue

            reads = self.taxon_mapping[taxon_ranks]


            dir_path = os.path.join(root, *taxon_ranks)
            if not os.path.exists(dir_path):
                os.makedirs(dir_path)

            try:
                fafp = open(os.path.join(dir_path,'reads.fa'), 'w')
                qualfp = open(os.path.join(dir_path,'reads.fa.qual'), 'w')
                for aligned_read in reads:
                    self._sam_to_fastx(aligned_read, fasta=fafp, qual=qualfp)
                #    samfp.write(aligned_read)

            finally:
                fafp.close()
                qualfp.close()
                #samfp.close()

        self.done_segregation = True




###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    #parser.add_argument('positional_arg', help="Required")
    parser.add_argument('taxonfile', 
                        help="Path to file containing the taxonomy of references")
    parser.add_argument('samfile', nargs='+', 
                        help="path to either sam or bam files for segregation")
    #parser.add_argument('-X', '--optional_X', action="store_true", default=False, help="flag")
    
    # parse the arguments
    args = parser.parse_args()        

    t = TaxonSegregator(args.taxonfile)

    # do what we came here to do
    for sam in args.samfile:
        t.parse_sam(sam)
    t.segregate()

###############################################################################
###############################################################################
###############################################################################
###############################################################################
