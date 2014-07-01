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
from __future__ import division

__author__ = "Connor Skennerton"
__copyright__ = "Copyright 2013"
__credits__ = ["Connor Skennerton"]
__license__ = "GPL3"
__version__ = "0.4.0"
__maintainer__ = "Connor Skennerton"
__email__ = "c.skennerton@gmail.com"
__status__ = "Development"

from qiime.util import get_rdp_jarpath
from os import system, remove, path, mkdir
from os.path import split, splitext
import tempfile
import logging
from qiime.assign_taxonomy import (
        BlastTaxonAssigner, MothurTaxonAssigner, RdpTaxonAssigner,
        RtaxTaxonAssigner, Tax2TreeTaxonAssigner, validate_rdp_version,
        TaxonAssigner
        )

from amplishot.parse.sam import SamFileReader
from amplishot.app.bowtie import Bowtie2

###############################################################################
class BowtieTaxonAssigner(TaxonAssigner):
    ''' Assign taxon by mapping with bowtie2
    '''
    def __init__(self, params):
        _params = {
                'index': None,
                'percentId': 0.85,
                'id_to_taxonomy_fp': None,
                'threads': 1,
                'fasta': True,
                }
        print self.Params
        _params.update(params)
        print self.Params
        TaxonAssigner.__init__(self, _params)
        print self.Params
    
    def __call__(self, seq_path, result_path=None, log_path=None,
            sam_path=None):
        print self.Params
        if self.Params['index'] is None or self.Params['id_to_taxon_fp'] is None:
            raise ValueError('To use the Bowtie taxon assigner both a bowtie2'\
                    ' index file must be provided and a file containing taxon'\
                    ' sequence identifiers to taxonomy strings')

        try:
            logger = self._get_logger(log_path=log_path)

            if sam_path is not None:
                fp = open(sam_path, 'w+b')
                logger.info('Saving alignments in SAM format to: %s' % sam_path)
            else:
                fp = tempfile.TemporaryFile()
            bowtie_params = dict()
            bowtie_params['-U'] = seq_path
            bowtie_params['-x'] = self.Params['index']
            bowtie_params['-p'] = self.Params['threads']
            if self.Params['fasta']:
                bowtie_params['-f'] = True
            b = Bowtie2(params=bowtie_params)
            logger.info('running bowtie with cmd: %s' % str(b))
            results = b(stdout=fp)
            fp.seek(0)
            results = self._generate_taxon_map(fp)
            if result_path is not None:
                logger.info('writing taxonomy map to: %s' % result_path)
                with open(result_path, 'w') as fp:
                    for seq_id, data in results.items():
                        fp.write('%s\t' % seq_id)
                        fp.write('%s\n' % '\t'.join(map(str, data)))
        finally:
            fp.close()

        return results

    def _generate_taxon_map(self, sam_fp):
        taxon_map =\
        self._parse_id_to_taxonomy_file(open(self.Params['id_to_taxon_fp'],'U'))
        results = dict()
        sam_reader = SamFileReader(sam_fp, parseHeader=False)
        for alignment in sam_reader.parse():
            if not alignment.is_unmapped():
                gg_tax = taxon_map[alignment.rname]
                confidence = alignment.percent_identity()
            else:
                gg_tax = get_greengeses_tax_string(tuple())
                confidence = 0.0
            results[alignment.qname] = (gg_tax, confidence, alignment.rname,
                    alignment.cigar)

        return results

    def _get_logger(self, log_path=None):
            if log_path is not None:
                handler = logging.FileHandler(log_path, mode='w')
            else:
                class NullHandler(logging.Handler):
                    def emit(self, record): pass
                handler = NullHandler()
            logger = logging.getLogger("BowtieTaxonAssigner logger")
            logger.addHandler(handler)
            logger.setLevel(logging.INFO)
            return logger

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
    return '; '.join(greengenes_format(taxon_ranks))


assignment_method_constructors = {
        'blast': BlastTaxonAssigner,
        'mothur': MothurTaxonAssigner,
        'rdp': RdpTaxonAssigner,
        'rtax': RtaxTaxonAssigner,
        'tax2tree': Tax2TreeTaxonAssigner,
        'bowtie': BowtieTaxonAssigner
        }

assignment_method_choices = ['rdp','blast','rtax','mothur', 'tax2tree', 'bowtie']

def assign_taxonomy(assigner, infile, config, result_file=None,
        log_file=None):
    ''' This is the global entry point for taxonomy assignment
        This function is essentially a copy of the Qiime assign taxonomy but
        plugged into the amplishot configuration module.  There is also the
        ability to use bowtie for taxonomy assignment
    '''
    params = dict()
    if assigner not in assignment_method_choices:
        raise ValueError('your choice of taxonomy assigner is invalid.'\
                ' please choose one of the following %s' %\
                str(assignment_method_choices))
    
    try:
        params = config.data[assigner]  #{ 'blast_db': config.data['blast_db'] }
    except KeyError, e:
        if assigner == 'bowtie':
            params['index'] = config.data['mapper_database']
            params['id_to_taxonomy_fp'] = config.data['taxonomy_file']
        else:
            raise e
    finally:
        if 'id_to_taxonomy_fp' not in params:
            params['id_to_taxonomy_fp'] = config.data['taxonomy_file']
        if assigner == 'blast':
            params['id_to_taxonomy_filepath'] = params['id_to_taxonomy_fp']
            try:
                params['Max E value'] = params['evalue']
            except KeyError:
                pass
            else:
                del params['evalue']

    taxon_assigner_constructor =\
         assignment_method_constructors[assigner]
    taxon_assigner = taxon_assigner_constructor(params)
    if assigner == 'bowtie':
        return taxon_assigner(infile,\
                     result_path=result_file, log_path=log_file,
                     sam_path=path.join(config.data['output_directory'],
                         'full_length_alignments.sam')
                     )
    else:
        return taxon_assigner(infile,\
                     result_path=result_file, log_path=log_file)
