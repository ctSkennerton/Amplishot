#!/usr/bin/env python
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


# the following is a port of the 16S rRNA orientation checker found in the
# ReadSeq package from RDP (https://github.com/rdpstaff/ReadSeq)
# it is a translation of the classes found in edu.msu.cme.rdp.readseq.utils.orientation
import os

class TrainingInfo(object):

    def __init__(self, datafile):

        self.NUM_OF_WORDS = 0x10000
        self.RNA_BASES = 4
        self.intComplementLookup = (1,0,3,2)
        self.wordPairPriorDiffArr = [None] * self.NUM_OF_WORDS

        self.createLogWordPriorArr(datafile)

    def createLogWordPriorArr(self, datafile):
        logWordPriorArr = [None] * self.NUM_OF_WORDS
        with open(datafile) as fp:
            for line in fp:
                (index, logprior) = line.split()
                logWordPriorArr[int(index)] = float(logprior)

        origWord = [0] * 8
        self.generateWordPairDiffArr(logWordPriorArr, origWord, 0);


    def generateWordPairDiffArr(self,  logWordPriorArr, word, beginIndex):
        if(beginIndex < 0 or beginIndex > len(word)):
            return

        origWordIndex = self.getWordIndex(word)
        revWordIndex = self.getWordIndex(self.getReversedWord(word))

        origWordPrior = logWordPriorArr[origWordIndex]
        revWordPrior = logWordPriorArr[revWordIndex]
        self.wordPairPriorDiffArr[origWordIndex] = origWordPrior - revWordPrior

        for i in range(beginIndex, len(word)):
            origBase = word[i]
            for j in range(4):
                if(word[i] != j):
                    word[i] = j
                    self.generateWordPairDiffArr(logWordPriorArr, word, i + 1)
                    word[i] = origBase

    def getWordIndex(self, word):
        wordIndex = 0;
        for w in range(len(word)):
            wordIndex <<= 2
            wordIndex &= 0xffff
            wordIndex |= word[w]

        return wordIndex


    def getReversedWord(self, word):
        length = len(word)
        reverseWord = [None] * len(word)
        for w in range(length):
            reverseWord[length - 1 - w] = self.intComplementLookup[word[w]]

        return reverseWord;

    def getWordPairPriorDiff(self, wordIndex):
        return self.wordPairPriorDiffArr[wordIndex]


class WordGenerator(object):

    def __init__(self, seq):
        self.position = 0
        self.nextWord = 0
        self.seqString = seq
        self.invalid = False

        self.charLookup = [-1] * 128
        self.charLookup[65] = 0
        self.charLookup[84] = 1
        self.charLookup[85] = 1
        self.charLookup[71] = 2
        self.charLookup[67] = 3
        self.charLookup[97] = 0
        self.charLookup[116] = 1
        self.charLookup[117] = 1
        self.charLookup[103] = 2
        self.charLookup[99] = 3

    def hasNext(self):
        baseCount = 0
        charIndex = 0
        while baseCount < 8 and self.position < len(self.seqString) - 1:
            self.position += 1
            try:
                nextBase = self.seqString[self.position]
            except IndexError, e:
                print "gone past the end: %d\t%d" % (len(self.seqString), self.position)
                raise e

            charIndex = self.charLookup[ord(nextBase)]
            if(charIndex == -1):
                baseCount = -1
                charIndex = 0
                self.invalid = True

            baseCount += 1
            self.nextWord <<= 2
            self.nextWord &= 0xffff
            self.nextWord |= charIndex

        if(baseCount < 8):
            return False;
        else:
            return True;

    def next(self):
        if(self.invalid):
            raise Exception("Attempt to call WordGenerator.next() when no more words.");
        else:
            return self.nextWord


class OrientationChecker(object):
    def __init__(self, data=os.path.join(os.path.dirname(__file__),"data", "logWordPrior.txt")):
        self.traininfo = TrainingInfo(data)

    def isSeqReversed(self, seqString):

        generator = WordGenerator(seqString)
        reverse = False
        priorDiff = 0.0
        wordIndex = 0
        while generator.hasNext():
            wordIndex = generator.next()
            priorDiff += self.traininfo.getWordPairPriorDiff(wordIndex)

        if (priorDiff < 0.0):
            reverse = True

        return reverse

if __name__ == '__main__':
    import argparse
    import sys
    from Bio import SeqIO

    parser = argparse.ArgumentParser()
    parser.add_argument('fasta', help='file containing fasta sequences to orietate')
    parser.add_argument('priors', help='file containing prior probabilities of 8mers')

    args = parser.parse_args()

    checker = OrientationChecker(args.priors)
    for seq_record in SeqIO.parse(args.fasta, "fasta"):
        if checker.isSeqReversed(seq_record.seq):
            seq_record.reverse_complement()
            #print "reversed"
        SeqIO.write(seq_record, sys.stdout, 'fasta')
