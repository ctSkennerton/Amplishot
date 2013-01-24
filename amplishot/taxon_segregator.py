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
import os
import re
import Bio.Seq

import amplishot.parse.sam
import amplishot.parse.fastx
###############################################################################
###############################################################################
###############################################################################
###############################################################################
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
            print "taxonomy file is ggs"
            self._parse_greengenes_sparate(taxonfp)
        elif format == 'ggc':
            print "taxonomy file is ggc"
            self._parse_greengenes_combined(taxonfp)
        elif format == 'sil':
            print "taxonomy file is sil"
            self._parse_silva(taxonfp)

    def _check_taxon_coverage(self, taxon, minCount, minCoverage):
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
            seq = Bio.Seq.reverse_complement(alignedRead.seq)
            quality = quality[::-1]

        if fasta is not None:
            fasta.write('>%s\n%s\n' %(name, seq))

        if fastq is not None:
            fastq.write('@%s\n%s\n+\n%s' % (name,seq,quality))
            
        if qual is not None:
            quality = amplishot.parse.fastx.decode_quality(quality)
            quality = ' '.join(str(x) for x in quality)
            qual.write('>%s\n%s\n' %(name, quality))


    def parse_sam(self, sam, percentId = 0.98):
        """iterate through the records in a samfile and place them into 
           one of the taxonomies based on the mapping

           sam: an opened samfile generated with pysam 
           percentId: The percent identity between the read and the reference
        """
        if self.done_segregation:
            raise RuntimeError, 'Segregation has already taken place.  Please\
            parse all samfiles at one and then call segregate at the end'
        samf = amplishot.parse.sam.SamFileReader(sam, parseHeader=False)
        for read in samf():
            if not read.is_unmapped():
                if float(read.qlen - read.tags['NM'])/ float(read.qlen) >= percentId:
                    t = self.ref_taxon_mapping[read.rname]
                    self.taxon_mapping[t].append(read)

    def segregate(self, root='root', mergeUpLevel=5, minCount=1000,
            minCoverage=2, fasta=True, qual=True,
            fastq=False, sam=False, prefix='reads'):
        """Partition all reads into separate files
           Ideally this should be called after all sam files have been parsed
           This function will make a directory structure equal to the taxonomy
           strings for each of the reference sequences.  It will overwrite any
           files already in the directory structure.
           Returns two lists of tuples that contain taxons with suitable
           coverage and those that do not

           root: The root directory name.  By default creates a directory called 
                 'root' in the current directory
           mergeUpLevel: The taxonomic level at which to stop merging up if
                there is not enough coverage.  By default this is level 5, which in
                greengenes is the family level
           minCount: the minimum number of positions that must have coverage
           minCoverage: the minimum coverage allowed for the covered positions
           fasta: output a fasta file containing segregated reads
           qual: output a fasta.qual file containing integer quality values for
                segregated reads
           fastq: output a fastq file containing seqregated reads
           sam: output segregated reads in sam format
           prefix: the name for the output files without extension (don't put
                directories here, they are created automatically)
        """
        complete_taxons = []
        incomplete_taxons = []
        sorted_taxons = self.taxon_mapping.keys()
        sorted_taxons = sorted(sorted_taxons, reverse=True, cmp=lambda x,y: cmp(len(x), len(y)))
        for taxon_ranks in sorted_taxons:
            if not self._check_taxon_coverage(taxon_ranks, monCount,
                    minCoverage):
                if len(taxon_ranks) > mergeUpLevel:
                    new_taxon = []
                    for i in range(len(taxon_ranks) - 1):
                        new_taxon.append(taxon_ranks[i])
                    
                    t = tuple(new_taxon)
                    try:
                        self.taxon_mapping[t].extend(self.taxon_mapping[taxon_ranks])
                    except KeyError:
                        # no taxon above this level therefore no point in
                        # deleting this taxon
                        pass
                    else:
                        del self.taxon_mapping[taxon_ranks]
                else:
                    incomplete_taxons.append(taxon_ranks)
            else:
                complete_taxons.append(taxon_ranks)
                reads = self.taxon_mapping[taxon_ranks]

                dir_path = os.path.join(root, *taxon_ranks)
                if not os.path.exists(dir_path):
                    os.makedirs(dir_path)

                fafp = None
                qualfp = None
                fqfp = None
                samfp = None
                try:
                    if fasta:
                        fafp = open(os.path.join(dir_path,'%s.fa' % prefix), 'w')
                    if qual:
                        qualfp = open(os.path.join(dir_path,'%s.fa.qual' % prefix), 'w')
                    if fastq:
                        fqfp = open(os.path.join(dir_path,'%s.fq' % prefix), 'w')
                    if sam:
                        samfp = open(os.path.join(dir_path, '%s.sam' % prefix),
                                'w')

                    for aligned_read in reads:
                        self._sam_to_fastx(aligned_read, fasta=fafp,
                                qual=qualfp, fastq=fqfp)
                        if sam:
                            samfp.write('%s\n' % str(aligned_read))

                finally:
                    if fafp is not None: 
                        fafp.close()
                    if qualfp is not None:
                        qualfp.close()
                    if fqfp is not None:
                        fqfp.close()
                    if samfp is not None:
                        samfp.close()

        self.done_segregation = True
        return complete_taxons, incomplete_taxons

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
