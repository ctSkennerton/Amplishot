#!/usr/bin/env python
###############################################################################
#
# sam.py - simple parsing of samfiles - much faster than pysam
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

import subprocess
import tempfile
###############################################################################
###############################################################################
###############################################################################
###############################################################################

class SamRead(object):
    
    def __init__(self, fields):
        super(SamRead, self).__init__()
        if isinstance(fields, str):
            fields = fields.split('\t', 11)
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

class SamFile(object):
    def __init__(self, filepath):
        super(SamFile, self).__init__()
        self.filepath = filepath
        self.fp = None
        self._open()

    def _open(self):
        self.fp = tempfile.TemporaryFile()
        if self.filepath.endswith('sam'):
            subprocess.call(['samtools','view', '-S', self.filepath], stdout=alignments,
                    stderr=open(os.devnull, 'w'))
            self.fp.seek(1)
        else:
            subprocess.call(['samtools','view', self.filepath], stdout=alignments,
                    stderr=open(os.devnull, 'w'))
            self.fp.seek(1)
        
    def parse(self):
        for line in self.fp:
            line = line.rstrip()
            fields = line.split('\t', 11)
            yield SamRead(fields)
