#!/usr/bin/env python
###############################################################################
#
# taxon_segregator.py - take a sam file and a taxonomy and segregate the reads
#                       into their relative taxonomic divisions
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

###############################################################################
import logging
import sys
import os
import re
import operator
from collections import defaultdict

import amplishot.parse.sam
import amplishot.parse.fastx
from amplishot.util import reverse_complement
###############################################################################
###############################################################################
###############################################################################
###############################################################################

class Cluster (object):
    """Holding class for information about a cluster of reads
    """
    def __init__(self, _ref, _cutoff, _first=None, _second=None, _singles=None):
        """
            _ref: ggID for the reference sequence that this cluster is formed
            around
            _cutoff: The mapping cutoff for reads to the reference sequence in
            the original mapping
            _first: a list of AlignedRead objects that are first in matched
            pairs
            _second: a list of AlignedRead objects that are second in matched
            pairs
            _singles: a list of AlignedRead objects that are singletons
        """
        super(Cluster, self).__init__()
        self.ref = _ref
        self.cutoff = _cutoff
        self.first = _first
        self.second = _second
        self.singles = _singles
        if _first is None:
            self.first = []
        if _second is None:
            self.second = []
        if _singles is None:
            self.singles = []


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

    #unusual chars
    un_re = re.compile('[\s\[\]\{\}]')

    def __init__(self, taxonfile, cutoffs=[0.8,0.87,0.92,0.98],
            neighboursfile=None):
        """
           taxonfile: path to a file containing the mapping between reference
                      identifiers and the taxonomy string
        """
        super(TaxonSegregator, self).__init__()
        self.done_segregation = False

        self.taxon_mapping = dict()
        self.ref_taxon_mapping = dict()
        self.taxon_format = ''

        self.neighbours = dict()
        self.ggRefDistFile = neighboursfile

        cutoffs.sort()
        self.cutoffs = cutoffs

        for i in self.cutoffs:
            self.taxon_mapping[i] = {}
        self.taxon_mapping[-1] = {}

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

    def _remove_unusual_chars(self, tax):
        return self.un_re.sub('', tax)

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
                names.append(self._remove_unusual_chars(rank[3:]))

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
                names.append(self._remove_unusual_chars(rank[3:]))

            self._generate_mapping(refid, names)

    def _parse_silva(self, taxonfp):
        """silva formated taxonomy file - no taxon level designation
           12345    a;b;c;d;e;f;g
        """
        for line in taxonfp:
            line = line.rstrip()
            (refid, taxon_string) = line.split('\t')
            taxon_divisions = taxon_string.split(';')
            taxon_divisions = [self._remove_unusual_chars(x) for x in taxon_divisions]
            self._generate_mapping(refid, taxon_divisions)

    def _generate_mapping(self, refid, taxon_divisions):
        """take the input from one of the parsers and generate a mapping
           of tuple of taxon divisions to a reference id a the opposite
           reference id to tuple
        """
        taxon_divisions = filter(None, taxon_divisions)
        t = tuple(taxon_divisions)
        #self.taxon_mapping[t] = []
        self.ref_taxon_mapping[refid] = t

    def _parse_taxon_file(self, taxonfp):
        """firse check the file format then call the correct file parser
        """
        format = self._check_taxonomy_file_format(taxonfp)
        if format == 'ggs':
    #        print "taxonomy file is ggs"
            self._parse_greengenes_sparate(taxonfp)
        elif format == 'ggc':
    #        print "taxonomy file is ggc"
            self._parse_greengenes_combined(taxonfp)
        elif format == 'sil':
    #        print "taxonomy file is sil"
            self._parse_silva(taxonfp)

    def _check_taxon_coverage(self, taxon, minCount, minCoverage):
        """calculate the per base coverage for this taxon
           If a taxon does not contain sufficient coverage across the length of
           the reference sequences then it should be merged into a higher taxonomy
           Since the reference sequences should be of approximate length, there
           should not be too much problem when dealing with multiple references
           from the same taxon.  However it is likely that this coverage information
           will be a little bit 'fuzzy' since there is some variation

           taxon: A list containing aligned reads
           minCount: the minimum number of positions that must have coverage
           minCoverage: the minimum coverage allowed for the covered positions
        """
        #assume that the sequences are not greater than 1600bp
        coverage = [0]*5000

        for read in taxon:
            for i in range(read.pos, read.pos + read.qlen):
                coverage[i] += 1

        count = 0
        for pos in coverage:
            if pos > minCoverage:
                count += 1

        return count >= minCount

    def _sam_to_fastx2(self, cluster, prefix, fasta=True, qual=True, fastq=False,
            pairs=False, sep12=False):
        """Take a pysam AlignedRead and convert it into a fasta, qual or fastq

           cluster: A Cluster object which contains lists of paired and single
                    reads
           prefix:  file path prefix
        """
        if fasta and fastq:
            raise RuntimeError("cannot output both fasta and fastq at the same time, choose one or the other")
        #logging.debug("requrested the following read"
        #        "outputs:\nfasta=%s\nfastq=%s\npairs=%s\nsep12=%s\n",
        #        str(fasta), str(fastq), str(pairs), str(sep12))
        #logging.debug("cluster reads: first = %d, second = %d, singles = %d",
        #        len(cluster.first), len(cluster.second), len(cluster.singles))
        written_files = {}
        # for seqs - fasta, fastq
        fp1 = None
        fp2 = None
        fps = None
        fpp = None
        # for quals
        qfp1 = None
        qfp2 = None
        qfps = None
        qfpp = None
        if len(cluster.first):
            # we have reads in pairs
            if fasta:
                if pairs:
                    if sep12:
                        suffix1 = "_R1.fasta"
                        suffix2 = "_R2.fasta"
                        fp1 = open(prefix+suffix1, 'w')
                        fp2 = open(prefix+suffix2, 'w')
                        written_files['r1'] = suffix1
                        written_files['r2'] = suffix2
                    else:
                        suffix12 = "_R12.fasta"
                        fpp = open(prefix+suffix12, 'w')
                        written_files['r12'] = suffix12
                else:
                    suffixs = "_s.fasta"
                    fps = open(prefix+suffixs, 'w')
                    written_files['s'] = suffixs

            if qual:
                if pairs:
                    if sep12:
                        suffix1 = "_R1.fasta.qual"
                        suffix2 = "_R2.fasta.qual"
                        qfp1 = open(prefix+suffix1, 'w')
                        qfp2 = open(prefix+suffix2, 'w')
                        written_files['q1'] = suffix1
                        written_files['q2'] = suffix2
                    else:
                        suffix12 = "_R12.fasta.qual"
                        qfpp = open(prefix+suffix12, 'w')
                        written_files['q12'] = suffix12
                else:
                    suffixs = "_s.fasta.qual"
                    qfps = open(prefix+suffixs, 'w')
                    written_files['qs'] = suffixs

            if fastq:
                if pairs:
                    if sep12:
                        suffix1 = "_R1.fastq"
                        suffix2 = "_R2.fastq"
                        fp1 = open(prefix+suffix1, 'w')
                        fp2 = open(prefix+suffix2, 'w')
                        written_files['r1'] = suffix1
                        written_files['r2'] = suffix2
                    else:
                        suffix12 = "_R12.fastq"
                        fpp = open(prefix+suffix12, 'w')
                        written_files['r12'] = suffix12
                else:
                    suffixs = "_s.fastq"
                    fps = open(prefix+suffixs, 'w')
                    written_files['s'] = suffixs

            for ar1, ar2 in zip(cluster.first, cluster.second):
                name1 = ar1.qname
                seq1 = ar1.seq
                quality1 = ar1.qual
                name2 = ar2.qname
                seq2 = ar2.seq
                quality2 = ar2.qual

                if ar1.is_reversed():
                    seq1 = reverse_complement(seq1)
                    quality1 = quality1[::-1]
                if ar2.is_reversed():
                    seq2 = reverse_complement(seq2)
                    quality2 = quality2[::-1]

                if fp1 is not None and fp2 is not None:
                    if fasta:
                        fp1.write(">%s\n%s\n" % (name1, seq1))
                        fp2.write(">%s\n%s\n" % (name2, seq2))
                    if fastq:
                        fp1.write('@%s\n%s\n+\n%s\n' % (name1,seq1,quality1))
                        fp2.write('@%s\n%s\n+\n%s\n' % (name2,seq2,quality2))
                    if qual:
                        quality = amplishot.parse.fastx.decode_quality(quality1)
                        quality = ' '.join(str(x) for x in quality)
                        qfp1.write('>%s\n%s\n' %(name1, quality))

                        quality = amplishot.parse.fastx.decode_quality(quality2)
                        quality = ' '.join(str(x) for x in quality)
                        qfp2.write('>%s\n%s\n' %(name2, quality))
                elif fpp is not None:
                    if fasta:
                        fpp.write(">%s\n%s\n" % (name1, seq1))
                        fpp.write(">%s\n%s\n" % (name2, seq2))
                    if fastq:
                        fpp.write('@%s\n%s\n+\n%s\n' % (name1,seq1,quality1))
                        fpp.write('@%s\n%s\n+\n%s\n' % (name2,seq2,quality2))
                    if qual:
                        quality = amplishot.parse.fastx.decode_quality(quality1)
                        quality = ' '.join(str(x) for x in quality)
                        qfpp.write('>%s\n%s\n' %(name1, quality))

                        quality = amplishot.parse.fastx.decode_quality(quality2)
                        quality = ' '.join(str(x) for x in quality)
                        qfpp.write('>%s\n%s\n' %(name2, quality))
                elif fps is not None:
                    if fasta:
                        fps.write(">%s\n%s\n" % (name1, seq1))
                        fps.write(">%s\n%s\n" % (name2, seq2))
                    if fastq:
                        fps.write('@%s\n%s\n+\n%s\n' % (name1,seq1,quality1))
                        fps.write('@%s\n%s\n+\n%s\n' % (name2,seq2,quality2))
                    if qual:
                        quality = amplishot.parse.fastx.decode_quality(quality1)
                        quality = ' '.join(str(x) for x in quality)
                        qfps.write('>%s\n%s\n' %(name1, quality))

                        quality = amplishot.parse.fastx.decode_quality(quality2)
                        quality = ' '.join(str(x) for x in quality)
                        qfps.write('>%s\n%s\n' %(name2, quality))


        if len(cluster.singles):
            # only reads in singles
            if fasta:
                suffixs = "_s.fasta"
                if fps is None:
                    fps = open(prefix+suffixs, 'w')
                    written_files['s'] = suffixs
            if qual:
                if qfps is None:
                    qfps = open(prefix+"_s.fasta.qual", 'w')
                    written_files['qs'] = "_s.fasta.qual"
            if fastq:
                if fps is None:
                    fps = open(prefix+"_s.fastq",'w')
                    written_files['s'] = "_s.fastq"

            for ar in cluster.singles:
                name = ar.qname
                seq = ar.seq
                quality = ar.qual

                if ar.is_reversed():
                    seq = reverse_complement(seq)
                    quality = quality[::-1]

                if fasta:
                    fps.write(">%s\n%s\n" % (name, seq))
                if fastq:
                    fps.write('@%s\n%s\n+\n%s\n' % (name,seq,quality))
                if qual:
                    quality = amplishot.parse.fastx.decode_quality(quality)
                    quality = ' '.join(str(x) for x in quality)
                    qfps.write('>%s\n%s\n' %(name, quality))

        if fp1:
            fp1.close()
        if fp2:
            fp2.close()
        if fpp:
            fpp.close()
        if fps:
            fps.close()
        if qfp1:
            qfp1.close()
        if qfp2:
            qfp2.close()
        if qfpp:
            qfpp.close()
        if qfps:
            qfps.close()

        return written_files

    def clear(self):
        """ Delete all of the reads from the taxonomy hash
        """
        for key in self.taxon_mapping.keys():
            for i in self.taxon_mapping[key].keys():
                del self.taxon_mapping[key][i][:]
        self.done_segregation = False

    def sortHits(self):
        ''' Returns a matrix.  First dimension is the different cutoff
            values. The second dimension is a sorted list of ggIds
        '''
        result_arrays = [[] * n for n in range(len(self.cutoffs))]
        for cutoff_index, cutoff in enumerate(self.cutoffs):
            result_arrays[cutoff_index].extend(sorted(self.taxon_mapping[cutoff].keys(),
                key=lambda k: len(self.taxon_mapping[cutoff][k]), reverse=True))

        return result_arrays

    def getNeighbours(self, seqIdentityThreshld=0.03):
        # determine GG reference genes within sequence identity threshold
        for line in open(self.ggRefDistFile):
            lineSplit = line.split('\t')
            refId = lineSplit[0].rstrip()
            similarRefSeqs = [x.rstrip() for x in lineSplit[1:]]

            clusteredSeqs = set()
            for item in similarRefSeqs:
                itemSplit = item.split(':')
                if float(itemSplit[1]) < seqIdentityThreshld:
                    clusteredSeqs.add(itemSplit[0])

            self.neighbours[refId] = clusteredSeqs

    def separateReadsIntoConsistentPairs(self, clusteredReads,
            removeInconsistent=True):
        """ Take a single list of AlignedRead objects from a cluster of
        greengenes IDs and determine if any of the mates do not fall into the
        cluster and separate the reads into proper pairs and singleton reads
        """
        pairs = {}
        for ar in clusteredReads:
            if ar.has_multiple_segments():
                if  ar.rnext == '=' or ar.rnext in self.neighbours[ar.rname]:
                    if ar.is_first_segment():
                        if ar.qname in pairs:
                            pairs[ar.qname][0] = ar
                        else:
                            pairs[ar.qname] = [ar, None]
                    else:
                        if ar.qname in pairs:
                            pairs[ar.qname][1] = ar
                        else:
                            pairs[ar.qname] = [None, ar]
                elif not removeInconsistent:
                    pairs[ar.qname] = [ar, None]
            else:
                pairs[ar.qname] = [ar, None]

        first = []
        second = []
        singles = []
        for canonical_name, p in pairs.items():
            if p[1] is None:
                # single
                singles.append(p[0])
            elif p[0] is None:
                singles.append(p[1])
            else:
                first.append(p[0])
                second.append(p[1])
        return first, second, singles


    def clusterHits(self, sortedHits):
        # clusters reference sequences within sequence identity threshold
        # dict of cluster refID to cluster objects
        ggClusters = []
        for n in range(len(self.cutoffs)):
            processedIds = set()
            for i in xrange(0, len(sortedHits[n])):
                ggIdI = sortedHits[n][i]
                if ggIdI in processedIds:
                    continue

                processedIds.add(ggIdI)
                clusteredReads = []
                clusteredReads.extend(self.taxon_mapping[self.cutoffs[n]][ggIdI])

                # -1 means unmapped which means no taxon and no neighbours
                if ggIdI != -1:
                    for j in xrange(i + 1, len(sortedHits[n])):
                        ggIdJ = sortedHits[n][j]
                        if ggIdJ in processedIds:
                            continue

                        if ggIdJ in self.neighbours[ggIdI]:
                            clusteredReads.extend(self.taxon_mapping[self.cutoffs[n]][ggIdJ])
                            processedIds.add(ggIdJ)

                first, second, single = self.separateReadsIntoConsistentPairs(clusteredReads)
                clust = Cluster(ggIdI, self.cutoffs[n], first, second, single)
                ggClusters.append(clust)

        return ggClusters

    def parse_sam3(self, sam):
        ''' pysam version for using bam files
        '''
        track = pysam.Samfile(sam, "rb")
            for aln in track.fetch( until_eof = True ):
                perc = float(self.qlen - self.tags['NM']) / float(self.qlen)
                this_read_cutoff_index = 0
                for i in self.cutoffs:
                    if perc >= i:
                        this_read_cutoff_index = i
                    else:
                        break

                msr = amplishot.parse.sam.MiniSamRead(aln.qname,
                        aln.seq,
                        aln.qual,
                        track.getrname(aln.tid),
                        aln.flag,
                        aln.pos,
                        track.getrname(aln.rnext),
                        aln.pnext)

                try:
                    self.taxon_mapping[this_read_cutoff_index][msr.rname].append(msr)
                except KeyError:
                    self.taxon_mapping[this_read_cutoff_index][msr.rname] = [msr]


    def parse_sam2(self, sam):
        """iterate through the records in a samfile and place them into
           clusters based on the clusters of reference sequences. Requires
           that a mapping file be already parsed

           sam: an opened samfile generated with pysam
        """
        if self.done_segregation:
            raise RuntimeError, 'Segregation has already taken place.  Please\
            parse all samfiles at one and then call segregate at the end'

        samf = amplishot.parse.sam.SamFileReader(sam, parseHeader=False)
        for read in samf.parse():
            if not read.is_unmapped():
                percent_id = read.percent_identity()
                this_read_cutoff_index = 0
                for i in self.cutoffs:
                    if percent_id >= i:
                        this_read_cutoff_index = i
                    else:
                        break

                    try:
                        self.taxon_mapping[this_read_cutoff_index][read.rname].append(read)
                    except KeyError:
                        self.taxon_mapping[this_read_cutoff_index][read.rname] = [read]

    def segregate2(self, root='root', minCount=1000,
            minCoverage=2, fasta=True, qual=True,
            fastq=False, sam=False, prefix='reads',
            pairs=False, sep12=False):
        """Partition all reads into separate files
           Ideally this should be called after all sam files have been parsed
           This function will make a directory containing a number of files
           each containing reads associated with a particular cluster of OTUs.
           It will overwrite any files already in the directory structure.
           Returns a list of tuples that contain taxons with suitable
           coverage.

           root: The root directory name.  By default creates a directory called
                 'root' in the current directory
           minCount: the minimum number of positions that must have coverage
           minCoverage: the minimum coverage allowed for the covered positions
           fasta: output a fasta file containing segregated reads
           qual: output a fasta.qual file containing integer quality values for
                segregated reads
           fastq: output a fastq file containing seqregated reads
           sam: output segregated reads in sam format
           prefix: the name for the output files without extension (don't put
                directories here, they are created automatically)
           pairs: separate reads into paired and singleton reads
           sep12: separate paired reads into read1 and read2 files
        """
        if not self.neighbours:
            self.getNeighbours()

        sorted_hits = self.sortHits()
        otu_clusters = self.clusterHits(sorted_hits)
        complete_taxons = defaultdict(list)

        for cluster in otu_clusters:
            try:
                taxon_ranks = self.ref_taxon_mapping[cluster.ref]
            except:
                taxon_ranks = ()

            if self._check_taxon_coverage(cluster.first + cluster.second + cluster.singles, minCount, minCoverage):

                dir_path = os.path.join(root, *taxon_ranks)
                if not os.path.exists(dir_path):
                    os.makedirs(dir_path)

                prefix = os.path.join(dir_path,'%s_%.2f' % (str(cluster.ref), cluster.cutoff))
                suffixs = self._sam_to_fastx2(cluster, prefix, fasta=fasta, qual=qual,
                        sep12=sep12, fastq=fastq, pairs=pairs)

                complete_taxons[prefix] = suffixs

        self.done_segregation = True
        return complete_taxons

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
