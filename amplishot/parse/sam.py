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
###############################################################################
###############################################################################
###############################################################################
class SamReadError(Exception):
    pass


class AlignmentTagError(SamReadError):
    pass

class SamRead(object):
    
    def __init__(self, fields, parseTags=True):
        super(SamRead, self).__init__()
        if isinstance(fields, str):
            fields = fields.split('\t', 11)
            if len(fields) != 12:
                raise SamFileError
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
        self.tags = None
        self.tags_available = parseTags
        if parseTags:
            self.parse_tags(fields[11])
        else:
            self.tags = fields[11]

    def parse_tags(self, tagString):
        self.tags = dict()
        tags_split = tagString.split('\t')
        for tag in tags_split:
            tag_name, value_type, value = tag.split(":")
            if value_type == 'i':
                value = int(value)
            elif value_type == 'f':
                value = float(value)
            elif value_type == 'B':
                value = value.split(',')
                if value[0] == 'f':
                    value = [float(x) for x in value[1:]]
                else:
                    value = [int(x) for x in value[1:]]
            self.tags[tag_name] = value

    def _generate_tag_string(self):
        string = ''
        for name, value in self.tags.items():
            string += '%s:' % name
            if isinstance(value, int):
                string += 'i:%i' % value
            elif isinstance(value, float):
                string += 'f:%f' % value
            elif isinstance(value, list):
                string+= 'B:'
                if isinstance(value[0], float):
                    string += 'f'
                    string += ','.join(value)
                else:
                    string += 'i' + ','.join(value)
            else:
                string += 'Z:%s' % value
            string += '\t'

        return string.rstrip('\t')



    def __str__(self):
        string = '%s\t%i\t%s\t%i\t%i\t%s\t%s\t%i\t%i\t%s\t%s\t' % (self.qname,
                self.flag, self.rname, self.pos, self.mapq, self.cigar,
                self.rnext, self.pnext, self.tlen, self.seq, self.qual)
        if isinstance(self.tags, dict):
            string += self._generate_tag_string()
        else:
            string += self.tags
        return string
    
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

    def percent_identity(self):
        if not self.tags_available or not self.tags.has_key('NM'):
            raise AlignmentTagError('Cannot calculate percent identity as\
                    either the tags were not parsed or the NM tag is not\
                    present')
        return float(self.qlen - self.tags['NM']) / float(self.qlen)


class SamFileError(Exception): pass

class BamFileReader(object):
    pass


class SamFileReader(object):
    def __init__(self, f, parseHeader=True, parseTags=True, parseString=False):
        super(SamFileReader, self).__init__()
        try:
            if parseString:
                self.fp = f.splitlines()
            elif isinstance(f, file):
                self.fp = f
            else:
                self.fp = open(f)
            self.header = dict()
            self._parse_header(parseHeader)
        except OSError:
            raise SamFileError, 'Cannot open Samfile'

    def _parse_header_line(self, line):
        header_code = line[1:3]
        if header_code == 'HD':
            fields = line.split('\t')
            fields_dict = dict()
            for tag in fields[1:]:
                code, value = tag.split(':')
                fields_dict[code] = value
            self.header[header_code] = fields_dict
        elif header_code != 'CO':
            fields = line.split('\t')
            fields_dict = dict()
            for tag in fields[1:]:
                code, value = tag.split(':')
                fields_dict[code] = value
            try:
                self.header[header_code].append(fields_dict)
            except KeyError:
                self.header[header_code] = [fields_dict]
        else:
            try:
                self.header[header_code].append(line[3:])
            except KeyError:
                self.header[header_code] = [line[3:]]

    def _parse_header(self, saveInfo):
        for line in self.fp:
            line = line.rstrip()
            if line[0] != '@':
                break
            elif saveInfo:
                self._parse_header_line(line)


    def parse(self):
        for line in self.fp:
            fields = line.split('\t', 11)
            if len(fields) != 12:
                print 'malformed alignment.\nnumber of fields: %i\nline is: %s' % (len(fields), line)
                raise SamFileError
            yield SamRead(fields)
