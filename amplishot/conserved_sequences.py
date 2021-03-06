import logging
from cogent import DNA
import cogent.core.moltype
from amplishot.util import reverse_complement
from amplishot.search import wumanber
__author__ = "Connor Skennerton"
__copyright__ = "Copyright 2013"
__credits__ = ["Connor Skennerton"]
__license__ = "GPL3"
__version__ = "0.4.0"
__maintainer__ = "Connor Skennerton"
__email__ = "c.skennerton@gmail.com"
__status__ = "Development"

class ConservedSequence(object):
    def __init__(self, dnaseq, pos, rc=False):
        super(ConservedSequence,self).__init__()
        self.dnaseq = dnaseq
        self.pos = int(pos)
        self.rc = rc


class ConservedSequences(object):
    
    def __init__(self):
        super(ConservedSequences,self).__init__()
        self._CONSERVED_SEQUENCES = {
            '357': 'CTGAGAYACGGHCCARACTCCTACGGGAGGCAGCAG',
            '530': 'GGCTAAYTHYGTGCCAGCAGCCGCGGTAAKAC',
            '787': 'CRAACRGGATTAGATACCCYGGTAGTCCW',
            '926': 'AAACTYAAAKGAATTGRCGG',
            '1114': 'GGTGSTGCATGGYTGTCGTCAGCTCGTGYCGTGA',
            '1392': 'GGTGAATACGTTCCCGGGYCTTGYACNCAC'
        }
        self.conserved_sequences = {}
        self.positions = self._CONSERVED_SEQUENCES.keys()
        self._generate_unambiguous_sequences()
        self.wu = wumanber.WuManber(self.conserved_sequences.keys())

    def _disambiguate(self, sequence):
        index = sequence.firstDegenerate()
        if index is None:
            return sequence

        head = sequence[:index]
        tail = sequence[index+1:]
        ret_seqs = list()

        for unambiguous in cogent.core.moltype.IUPAC_DNA_ambiguities[sequence[index]]:
            ret = self._disambiguate(head+unambiguous+tail)
            if isinstance(ret, list):
                ret_seqs.extend(ret)
            else:
                ret_seqs.append(ret)
        return ret_seqs

    def _generate_unambiguous_sequences(self):
        unambiguous_conserved_sequences = dict()
        rev_unambiguous_conserved_sequences = dict()
        for pos,seq in self._CONSERVED_SEQUENCES.items():
            dnaseq = DNA.makeSequence(seq)
            ret = self._disambiguate(dnaseq)
            if isinstance(ret, list):
                for dnaseq_r in ret:
                    self.conserved_sequences[str(dnaseq_r)] = ConservedSequence(dnaseq_r, pos)
            else:
                self.conserved_sequences[str(ret)] = ConservedSequence(ret, pos)

        for seq,con_seq in self.conserved_sequences.items():
            rc_seq = DNA.makeSequence(seq)
            rc_seq.rc()
            self.conserved_sequences[str(rc_seq)] = ConservedSequence(rc_seq, con_seq.pos, 
                    rc=True)

    def orientate(self, seq):
        ret = self.wu.search_text(seq)
        if ret is not None:
            if self.conserved_sequences[ret[1]].rc:
                return reverse_complement(seq)
            else:
                return seq
        else:
            logging.error('cannot identify any 16S conserved sequences'\
                    'in:\n%s\nReturning original sequence',
                    (seq,))
            return seq


