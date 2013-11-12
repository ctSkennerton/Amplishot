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
__version__ = "0.3.3"
__maintainer__ = "Connor Skennerton"
__email__ = "c.skennerton@gmail.com"
__status__ = "Development"

###############################################################################
import logging
import sys
import os
import re
from collections import defaultdict

import amplishot.parse.sam
import amplishot.parse.fastx
from amplishot.util import reverse_complement
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

    #unusual chars
    un_re = re.compile('[\s\[\]\{\}]')
    def __init__(self, taxonfile, cutoffs=[0.8,0.87,0.92,0.98]):
        """
           taxonfile: path to a file containing the mapping between reference
                      identifiers and the taxonomy string
        """
        super(TaxonSegregator, self).__init__()
        self.done_segregation = False

        self.taxon_mapping = dict()
        self.ref_taxon_mapping = dict()
        self.taxon_format = ''

        cutoffs.sort()
        self.cutoffs = cutoffs

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
        self.taxon_mapping[t] = []
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

           taxon: A tuple containing the taxon string
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
            quality = amplishot.parse.fastx.decode_quality(quality)
            quality = ' '.join(str(x) for x in quality)
            qual.write('>%s\n%s\n' %(name, quality))


    def clear(self):
        """ Delete all of the reads from the taxonomy hash
        """
        for key in self.taxon_mapping.keys():
            for i in range(len(self.taxon_mapping[key])):
                del self.taxon_mapping[key][i][:]
        self.done_segregation = False

    def parse_sam(self, sam):
        """iterate through the records in a samfile and place them into 
           one of the taxonomies based on the mapping

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
                for i in range(len(self.cutoffs)):
                    if percent_id >= self.cutoffs[i]:
                        this_read_cutoff_index = i
                    else:
                        break

                    t = self.ref_taxon_mapping[read.rname]
                    try:
                        self.taxon_mapping[t][this_read_cutoff_index].append(read)
                    except IndexError:
                        tmp_array = [[] * n for n in range(len(self.cutoffs))]
                        self.taxon_mapping[t] = tmp_array
                        self.taxon_mapping[t][this_read_cutoff_index].append(read)
            else:
                try:
                    self.taxon_mapping[tuple()][0].append(read)
                except KeyError:
                    self.taxon_mapping[tuple()] = [[read]]


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
        complete_taxons = defaultdict(list)
        incomplete_taxons = defaultdict(list)
        sorted_taxons = self.taxon_mapping.keys()
        sorted_taxons = sorted(sorted_taxons, reverse=True, cmp=lambda x,y: cmp(len(x), len(y)))
        for taxon_ranks in sorted_taxons:
            for cutoff_index in range(len(self.taxon_mapping[taxon_ranks])):
                if not self._check_taxon_coverage(self.taxon_mapping[taxon_ranks][cutoff_index],\
                        minCount, minCoverage):
                    if len(taxon_ranks) > mergeUpLevel:
                        new_taxon = []
                        for i in range(len(taxon_ranks) - 1):
                            new_taxon.append(taxon_ranks[i])
                        
                        t = tuple(new_taxon)
                        try:
                            self.taxon_mapping[t][cutoff_index].extend(self.taxon_mapping[taxon_ranks][cutoff_index])
                        except (IndexError, KeyError):
                            # no taxon above this level therefore no point in
                            # deleting this taxon
                            pass
                        else:
                            del self.taxon_mapping[taxon_ranks][cutoff_index][:]
                    else:
                        incomplete_taxons[taxon_ranks].append(self.cutoffs[cutoff_index])
                else:
                    complete_taxons[taxon_ranks].append(self.cutoffs[cutoff_index])
                    reads = self.taxon_mapping[taxon_ranks][cutoff_index]

                    dir_path = os.path.join(root, *taxon_ranks)
                    if not os.path.exists(dir_path):
                        os.makedirs(dir_path)

                    fafp = None
                    qualfp = None
                    fqfp = None
                    samfp = None
                    try:
                        if fasta:
                            fafp = open(os.path.join(dir_path,'%s_%.2f.fa' %
                                (prefix, self.cutoffs[cutoff_index])), 'w')
                        if qual:
                            qualfp = open(os.path.join(dir_path,'%s_%.2f.fa.qual' %
                                (prefix, self.cutoffs[cutoff_index])), 'w')
                        if fastq:
                            fqfp = open(os.path.join(dir_path,'%s_%.2f.fq' %
                                (prefix, self.cutoffs[cutoff_index])), 'w')
                        if sam:
                            samfp = open(os.path.join(dir_path, '%s_%.2f.sam' %
                                (prefix, self.cutoffs[cutoff_index])), 'w')

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

class ReferenceSeqHit(object):
  def __init__(self, refId):
    self.refId = refId
    self.numHits = 0
    self.pairs = {}
    self.singles = {}

  def addPairHit(self, filename1, filename2, queryId1, queryId2):
    self.numHits += 2

    if filename1 not in self.pairs:
      self.pairs[filename1] = set()

    if filename2 not in self.pairs:
      self.pairs[filename2] = set()

    self.pairs[filename1].add(queryId1)
    self.pairs[filename2].add(queryId2)

  def addSingleHit(self, filename, queryId):
    self.numHits += 1

    if filename not in self.singles:
      self.singles[filename] = set()
    self.singles[filename].add(queryId)

class IdentifyRecoverable16S(object):
  def __init__(self):
    self.ggRefDist = '/srv/db/gg/2013_05/gg_13_5_otus/dist/dist_##_otus.tsv'

    self.bQuiet = False
    pass

  def ggIdFromTaxonomy(self, taxonomy):
    if '(' in taxonomy[7]:
      return taxonomy[7][4:taxonomy[7].rfind('(')]
    else:
      return taxonomy[7][4:].rstrip()

  def readClassifications(self, classificationFile):
    classifications = {}

    for line in open(classificationFile):
      lineSplit = line.split()
      seqId = lineSplit[0]
      taxa = [x.strip() for x in lineSplit[1].split(';') if x.strip() != '']

      classifications[seqId] = taxa

    return classifications

  def identifyConsistentPairs(self, referenceSeqHits, pairFile1, pairFile2, classificationFile1, classificationFile2, neighbours, bPairsAsSingles, bSingleEnded):
    if not self.bQuiet:
      print '  Reading classification file.'

    classifications1 = self.readClassifications(classificationFile1)
    classifications2 = self.readClassifications(classificationFile2)

    # get read ids
    if not self.bQuiet:
      print '  Identifying consistent pairs.'

    # find pairs that agree on classification
    numSingletons = 0
    numPairsInAgreement = 0
    numPairsInDisagreement = 0
    numUnclassified = 0
    for seqId1 in classifications1:
      ggId1 = ggId2 = None

      seqId2 = seqId1[0:-1] + '2'

      ggId1 = self.ggIdFromTaxonomy(classifications1[seqId1])
      if ggId1 == 'unclassified' or ggId1 == 'unmapped':
        numUnclassified += 1
  
        ggId2 = self.ggIdFromTaxonomy(classifications2[seqId2])
        if ggId2 == 'unclassified' or ggId2 == 'unmapped':
          numUnclassified += 1
        else:
          numSingletons += 1
          if bSingleEnded:
            referenceSeqHits[ggId2] = referenceSeqHits.get(ggId2, ReferenceSeqHit(ggId2))
            referenceSeqHits[ggId2].addSingleHit(pairFile2, seqId2)
  
        continue

      ggId2 = self.ggIdFromTaxonomy(classifications2[seqId2])
      if ggId2 == 'unclassified' or ggId2 == 'unmapped':
        numUnclassified += 1

        ggId1 = self.ggIdFromTaxonomy(classifications1[seqId1])
        numSingletons += 1
        if bSingleEnded:
          referenceSeqHits[ggId1] = referenceSeqHits.get(ggId1, ReferenceSeqHit(ggId1))
          referenceSeqHits[ggId1].addSingleHit(pairFile1, seqId1)

        continue

      if ggId1 == ggId2 or ggId1 in neighbours[ggId2]:
        referenceSeqHits[ggId1] = referenceSeqHits.get(ggId1, ReferenceSeqHit(ggId1))
        referenceSeqHits[ggId1].addPairHit(pairFile1, pairFile2, seqId1, seqId2)
        numPairsInAgreement += 1
      else:
        numPairsInDisagreement += 1
        if bPairsAsSingles:
          referenceSeqHits[ggId1] = referenceSeqHits.get(ggId1, ReferenceSeqHit(ggId1))
          referenceSeqHits[ggId1].addSingleHit(pairFile1, seqId1)
          referenceSeqHits[ggId2] = referenceSeqHits.get(ggId2, ReferenceSeqHit(ggId2))
          referenceSeqHits[ggId2].addSingleHit(pairFile2, seqId2)

    if not self.bQuiet:
      print '    Classified reads: ' + str(len(classifications1) + len(classifications2))
      print '      Singletons: ' + str(numSingletons) + ' (%.2f' % (numSingletons*100.0/(len(classifications1) + len(classifications2))) + ')'
      print '      Reads in pairs with similar classifications: ' + str(2*numPairsInAgreement) + ' (%.2f' % (2*numPairsInAgreement*100.0/(len(classifications1) + len(classifications2))) + ')'
      print '      Reads in pairs with different classifications: ' + str(2*numPairsInDisagreement) + ' (%.2f' % (2*numPairsInDisagreement*100.0/(len(classifications1) + len(classifications2))) + ')'
      print '      Unclassified/unmapped reads: ' + str(numUnclassified) + ' (%.2f' % (numUnclassified*100.0/(len(classifications1) + len(classifications2))) + ')'
      print ''

  def addSingletons(self, referenceSeqHits, single, classificationFile):
    if not self.bQuiet:
      print '  Reading classifications.'

    classifications = self.readClassifications(classificationFile)

    if not self.bQuiet:
      print '  Processing single-ended reads.'

    numUnclassified = 0
    for seqId in classifications:
      ggId = self.ggIdFromTaxonomy(classifications[seqId])

      if ggId == 'unclassified' or ggId == 'unmapped':
        numUnclassified += 1
      else:
        referenceSeqHits[ggId] = referenceSeqHits.get(ggId, ReferenceSeqHit(ggId))
        referenceSeqHits[ggId].addSingleHit(single, seqId)

    if not self.bQuiet:
      print '    Classified single-ended reads: ' + str(len(classifications) - numUnclassified)
      print '    Unclassified/unmapped single-ended reads: ' + str(numUnclassified)
      print ''

  def sortHits(self, referenceSeqHits):
    ggHits = {}
    for ggId in referenceSeqHits:
      ggHits[ggId] = referenceSeqHits[ggId].numHits

    sortedHits = sorted(ggHits.iteritems(), key=operator.itemgetter(1), reverse=True)

    return sortedHits

  def getNeighbours(self, ggRefDistFile, seqIdentityThreshld):
    # determine GG reference genes within sequence identity threshold
    neighbours = {}
    for line in open(ggRefDistFile):
      lineSplit = line.split('\t')
      refId = lineSplit[0].rstrip()
      similarRefSeqs = [x.rstrip() for x in lineSplit[1:]]

      clusteredSeqs = set()
      for item in similarRefSeqs:
        itemSplit = item.split(':')
        if float(itemSplit[1]) < seqIdentityThreshld:
          clusteredSeqs.add(itemSplit[0])

      neighbours[refId] = clusteredSeqs

    return neighbours

  def clusterHits(self, sortedHits, neighbours, minSeqCutoff):
    # clusters reference sequences within sequence identity threshold
    ggClusters = []
    processedIds = set()
    for i in xrange(0, len(sortedHits)):
      ggIdI = sortedHits[i][0]
      if ggIdI in processedIds:
        continue

      processedIds.add(ggIdI)

      clusteredIds = [ggIdI]
      totalHits = sortedHits[i][1]

      if not self.bQuiet:
        reportStr = '  Combining ' + str(ggIdI) + ' with ' + str(totalHits) + ' hits to: '

      for j in xrange(i+1, len(sortedHits)):
        ggIdJ = sortedHits[j][0]
        if ggIdJ in processedIds:
          continue

        if ggIdJ in neighbours[ggIdI]:
          clusteredIds.append(ggIdJ)
          totalHits += sortedHits[j][1]
          processedIds.add(ggIdJ)

          if not self.bQuiet:
            reportStr += str(ggIdJ) + ' (' + str(sortedHits[j][1]) + ' hits) '

      if not self.bQuiet and totalHits > minSeqCutoff:
        print reportStr

      ggClusters.append([totalHits, clusteredIds])

    return ggClusters

  def extractClusteredReads(self, clusters, referenceSeqHits, minSeqCutoff):
    readsInFiles = {}

    for cluster in clusters:
      hits = cluster[0]
      ggIds = cluster[1]

      if hits > minSeqCutoff:
        for ggId in ggIds:
          refSeqHit = referenceSeqHits[ggId]

          for pairFile in refSeqHit.pairs:
            readsInFiles[pairFile] = readsInFiles.get(pairFile, set()).union(refSeqHit.pairs[pairFile])

          for singleFile in refSeqHit.singles:
            readsInFiles[singleFile] = readsInFiles.get(singleFile, set()).union(refSeqHit.singles[singleFile])

    seqsInFiles = {}
    for filename in readsInFiles:
      seqIds = set()
      for readId in  readsInFiles[filename]:
        seqIds.add(readId[0:readId.rfind('/')])
      seqsInFiles[filename] = extractSeqs(filename, seqIds)
      
    return seqsInFiles

  def extractRecoverable16S(self, referenceSeqHits, neighbours, minSeqCutoff, outputDir):
    # sort hits to GreenGene reference sequences
    if not self.bQuiet:
      print 'Sorting reference sequences by number of hits.'

    sortedHits = self.sortHits(referenceSeqHits)

    if not self.bQuiet:
      print '  Initial clusters: ' + str(len(sortedHits)) + '\n'

    # greedily combine hits to similar GreenGene reference sequences
    if not self.bQuiet:
      print 'Clustering similar reference sequences.'

    ggClusters = self.clusterHits(sortedHits, neighbours, minSeqCutoff)

    if not self.bQuiet:
      print '\n  Final clusters: ' + str(len(ggClusters)) + '\n'

    # get sequences within clusters
    if not self.bQuiet:
      print 'Extracting clustered reads from file.\n'

    seqsInFiles = self.extractClusteredReads(ggClusters, referenceSeqHits, minSeqCutoff)

    # write out clusters containing a sufficient number of hits
    if not self.bQuiet:
      print 'Writing out reads belonging to putative 16S genes: '

    numRecoverable16S = 0
    for cluster in ggClusters:
      hits = cluster[0]
      ggIds = cluster[1]

      if hits > minSeqCutoff:
        numRecoverable16S += 1

        if not self.bQuiet:
          print '  Cluster ' + str(numRecoverable16S) + ': ' + ggIds[0] + ' (' + str(hits) + ' reads)'

        # write out singletons
        singletonsOut = open(outputDir + ggIds[0] + '.singletons.fasta', 'w')
        pairOut1 = open(outputDir + ggIds[0] + '.1.fasta', 'w')
        pairOut2 = open(outputDir + ggIds[0] + '.2.fasta', 'w')
        
        for ggId in ggIds:
          refSeqHit = referenceSeqHits[ggId]

          for singleFile in refSeqHit.singles:
            for readId in refSeqHit.singles[singleFile]:
              seqId = readId[0:readId.rfind('/')]
              singletonsOut.write('>' + readId + '\n')
              singletonsOut.write(seqsInFiles[singleFile][seqId][1] + '\n')

          for pairFile in sorted(list(refSeqHit.pairs)):
            for readId in sorted(list(refSeqHit.pairs[pairFile])): # ensure reads are in the same order for both files
              seqId = readId[0:readId.rfind('/')]
              if '/1' in readId:
                pairOut1.write('>' + readId + '\n')
                pairOut1.write(seqsInFiles[pairFile][seqId][1] + '\n')
              elif '/2' in readId:
                pairOut2.write('>' + readId + '\n')
                pairOut2.write(seqsInFiles[pairFile][seqId][1] + '\n')
              else:
                print "[Error] Unrecognized file format. Pairs must be denoted by '/1' and '/2'"
                sys.exit()

        singletonsOut.close()
        pairOut1.close()
        pairOut2.close()

    if not self.bQuiet:
      print ''
      print '  Number of recoverable 16S genes identified: ' + str(numRecoverable16S)

  def run(self, configFile, otu, seqIdentityThreshold, minSeqCutoff, bPairsAsSingles, bSingleEnded, bQuiet):
    self.bQuiet = bQuiet

    rc = ReadConfig()
    projectParams, sampleParams = rc.readConfig(configFile, outputDirExists = True)

    ggRefDistFile = self.ggRefDist.replace('##', str(otu))
    neighbours = self.getNeighbours(ggRefDistFile, seqIdentityThreshold)

    # create directory to store putative 16S genes
    dirPutative16S = projectParams['output_dir'] + 'putativeSSU/'
    if not os.path.exists(dirPutative16S):
      os.makedirs(dirPutative16S)
    else:
      rtn = raw_input('Remove previously recovered 16S reads (Y or N)? ')
      if rtn.lower() == 'y' or rtn.lower() == 'yes':
        files = os.listdir(dirPutative16S)
        for f in files:
          if f.endswith('fasta'):
            os.remove(dirPutative16S + '/' + f)
      else:
        sys.exit()

    referenceSeqHits = {}
    for sample in sampleParams:
      if not self.bQuiet:
        print ''
        print sample + ':'

      extractedPrefix = projectParams['output_dir'] + 'extracted/' + sample
      classifiedPrefix = projectParams['output_dir'] + 'classified/' + sample
      pairs = sampleParams[sample]['pairs']
      singles = sampleParams[sample]['singles']

      for i in xrange(0, len(pairs), 2):
        pair1Base = ntpath.basename(pairs[i])
        pair2Base = ntpath.basename(pairs[i+1])
        
        classificationFile1 = classifiedPrefix + '.' + pair1Base[0:pair1Base.rfind('.')] + '.intersect.16S.tsv'
        classificationFile2 = classifiedPrefix + '.' + pair2Base[0:pair2Base.rfind('.')] + '.intersect.16S.tsv'

        if not self.bQuiet:
          print '  Processing files: '
          print '    ' + classificationFile1
          print '    ' + classificationFile2

        pairFile1 = extractedPrefix + '.' + pair1Base[0:pair1Base.rfind('.')] + '.intersect.SSU.fasta'
        pairFile2 = extractedPrefix + '.' + pair2Base[0:pair2Base.rfind('.')] + '.intersect.SSU.fasta'

        self.identifyConsistentPairs(referenceSeqHits, pairFile1, pairFile2, classificationFile1, classificationFile2, neighbours, bPairsAsSingles, bSingleEnded)
        
        if bSingleEnded:
          classificationFile = classifiedPrefix + '.' + pair1Base[0:pair1Base.rfind('.')] + '.difference.16S.tsv'
          
          if not self.bQuiet:
            print '  Processing file: ' + classificationFile
          
          singleFile = extractedPrefix + '.' + pair1Base[0:pair1Base.rfind('.')] + '.difference.SSU.fasta'
          self.addSingletons(referenceSeqHits, singleFile, classificationFile)
          
      if bSingleEnded:
        for single in singles:
          singleBase = ntpath.basename(single)
          classificationFile = classifiedPrefix + '.' + singleBase[0:singleBase.rfind('.')] + '.16S.tsv'

          if not self.bQuiet:
            print '  Processing file: ' + classificationFile

          singleFile = extractedPrefix + '.' + singleBase[0:singleBase.rfind('.')] + '.SSU.fasta'
          self.addSingletons(referenceSeqHits, singleFile, classificationFile)

    self.extractRecoverable16S(referenceSeqHits, neighbours, minSeqCutoff, dirPutative16S)
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
