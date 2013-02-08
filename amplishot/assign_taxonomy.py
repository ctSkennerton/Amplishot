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
__version__ = "0.1.0"
__maintainer__ = "Connor Skennerton"
__email__ = "c.skennerton@gmail.com"
__status__ = "Development"

###############################################################################
class GreengenesFormatingError(Exception):
    pass

GREENGENES_PREFIXES = ['k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']

def greengenes_format(taxon_ranks):
    if len(taxon_ranks) > 7:
        raise GreengenesFormatingError, "There must be exactly 7 levels of\
    taxonomy for correctly formated greengenes tax string. %i found" %\
        len(taxon_ranks)
    tax = [''] * 7
    for i in range(7):
        if i < len(taxon_ranks):
            tax[i] = GREENGENES_PREFIXES[i] + taxon_ranks[i]
        else:
            tax[i] = GREENGENES_PREFIXES[i]
    return tax


def get_greengenes_tax_string(taxon_ranks):
    return ';'.join(greengenes_format(taxon_ranks))
