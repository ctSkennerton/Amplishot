#!/usr/bin/env python
###############################################################################
#
# util.py - simple functions needed throughout amplishot
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

# Code copied from Biopython
# Copyright 2000-2002 Brad Chapman.
# Copyright 2004-2005 by M de Hoon.
# Copyright 2007-2010 by Peter Cock.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import string
import sys

ambiguous_dna_complement = {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
        "M": "K",
        "R": "Y",
        "W": "W",
        "S": "S",
        "Y": "R",
        "K": "M",
        "V": "B",
        "H": "D",
        "D": "H",
        "B": "V",
        "X": "X",
        "N": "N",
        }
def _maketrans(complement_mapping):
    """Makes a python string translation table (PRIVATE).

    Arguments:
     - complement_mapping - a dictionary such as ambiguous_dna_complement
       and ambiguous_rna_complement from Data.IUPACData.

    Returns a translation table (a string of length 256) for use with the
    python string's translate method to use in a (reverse) complement.

    Compatible with lower case and upper case sequences.

    For internal use only.
    """
    before = ''.join(complement_mapping.keys())
    after = ''.join(complement_mapping.values())
    before = before + before.lower()
    after = after + after.lower()
    if sys.version_info[0] == 3:
        return str.maketrans(before, after)
    else:
        return string.maketrans(before, after)

_dna_complement_table = _maketrans(ambiguous_dna_complement)

def reverse_complement(sequence):
    return sequence.translate(_dna_complement_table)[::-1]
