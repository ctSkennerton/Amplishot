#!/usr/bin/env python
###############################################################################
#
# bowtie.py - app controller for bowtie
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
__version__ = "0.1.2"
__maintainer__ = "Connor Skennerton"
__email__ = "c.skennerton@gmail.com"
__status__ = "Development"

###############################################################################

from cogent.app.parameters import FlagParameter, ValuedParameter
from amplishot.app.util import ExtendedCommandLineApplication,\
ValuedParameter


class Bowtie2Build(ExtendedCommandLineApplication):
    """ Simple wrapper for bowtie index
    """
    _parameters = {
            '-f': FlagParameter('-', 'f'),
            '-c': FlagParameter('-', 'c'),
            '-a': FlagParameter('-', 'a'),
            '-p': FlagParameter('-', 'p'),
            '--bmax': ValuedParameter('--', 'bmax', Delimiter=' '),
            '--bmaxdivn': ValuedParameter('--', 'bmaxdivn', Delimiter=' '),
            '--dcv': ValuedParameter('--', 'dcv', Delimiter=' '),
            '--nodc': ValuedParameter('--', 'nodc', Delimiter=' '),
            '-r': FlagParameter('-', 'r'),
            '-3': FlagParameter('-', '3'),
            '-o': ValuedParameter('-', 'o', Delimiter=' '),
            '-t': ValuedParameter('-', 't', Delimiter=' '),
            '--seed': ValuedParameter('--', 'seed', Delimiter=' '),
            '--cutoff': ValuedParameter('--', 'cutoff', Delimiter=' '),
            '-q': FlagParameter('-', 'q'),
            '-h': FlagParameter('-', 'h'),
            '--version': FlagParameter('--', 'version')
            }

    _command = 'bowtie2-build'
    _input_handler = '_input_as_paths'


class Bowtie2(ExtendedCommandLineApplication):
    """ Simple bowtie2 application controller
    """
    _parameters = {
            '-x': ValuedParameter('-', 'x', Delimiter=' ', IsPath=True),
            '-1': ValuedParameter('-', '1', Delimiter=' ', IsPath=True),
            '-2': ValuedParameter('-', '2', Delimiter=' ', IsPath=True),
            '-U': ValuedParameter('-', 'U', Delimiter=' ', IsPath=True),
            '-S': ValuedParameter('-', 'S', Delimiter=' ', IsPath=True),
            '-q': FlagParameter('-', 'q'),
            '-f': FlagParameter('-', 'f'),
            '-r': FlagParameter('-', 'r'),
            '-c': FlagParameter('-', 'c'),
            '-s': ValuedParameter('-', 's', Delimiter=' '),
            '-u': ValuedParameter('-', 'u', Delimiter=' '),
            '-5': ValuedParameter('-', '5', Delimiter=' '),
            '-3': ValuedParameter('-', '3', Delimiter=' '),
            '--phred33': FlagParameter('--', 'phred33'),
            '--phred64': FlagParameter('--', 'phred64'),
            '--int-quals': FlagParameter('--', 'int-quals'),
            '--very-fast': FlagParameter('--', 'very-fast'),
            '--fast': FlagParameter('--', 'fast'),
            '--sensitive': FlagParameter('--', 'sensitive'),
            '--very-sensitive': FlagParameter('--', 'very-sensitive'),
            '--very-fast-local': FlagParameter('--', 'very-fast-local'),
            '--fast-local': FlagParameter('--', 'fast-local'),
            '--sensitive-local': FlagParameter('--', 'sensitive-local'),
            '--very-sensitive-local': FlagParameter('--',
                'very-sensitive-local'),
            '-N': ValuedParameter('-', 'N', Delimiter=' '),
            '-L': ValuedParameter('-', 'L', Delimiter=' '),
            '-i': ValuedParameter('-', 'i', Delimiter=' '),
            '--n-ceil': ValuedParameter('--', 'n-ceil', Delimiter=' '),
            '--dpad': ValuedParameter('--', 'dpad', Delimiter=' '),
            '--gbar': ValuedParameter('--', 'gbar', Delimiter=' '),
            '--ignore-quals': ValuedParameter('--', 'ignore-quals', Delimiter=' '),
            '--nofw': ValuedParameter('--', 'nofw', Delimiter=' '),
            '--norc': ValuedParameter('--', 'norc', Delimiter=' '),
            '--end-to-end': FlagParameter('--', 'end-to-end'),
            '--local': FlagParameter('--', 'local'),
            '--ma': ValuedParameter('--', 'ma', Delimiter=' '),
            '--mp': ValuedParameter('--', 'mp', Delimiter=' '),
            '--np': ValuedParameter('--', 'np', Delimiter=' '),
            '--rdg': ValuedParameter('--', 'rdg', Delimiter=' '),
            '--rfg': ValuedParameter('--', 'rfg', Delimiter=' '),
            '--score-min': ValuedParameter('--', 'score-min', Delimiter=' '),
            '-k': ValuedParameter('-', 'k', Delimiter=' '),
            '-a': FlagParameter('-', 'a'),
            '-D': ValuedParameter('-', 'D', Delimiter=' '),
            '-R': ValuedParameter('-', 'R', Delimiter=' '),
            '-I': ValuedParameter('-', 'I', Delimiter=' '),
            '-X': ValuedParameter('-', 'X', Delimiter=' '),
            '--fr': FlagParameter('--', 'fr'),
            '--rf': FlagParameter('--', 'rf'),
            '--ff': FlagParameter('--', 'ff'),
            '--no-mixed': FlagParameter('--', 'no-mixed'),
            '--no-discordant': FlagParameter('--', 'no-discordant'),
            '--no-dovetail': FlagParameter('--', 'no-dovetail'),
            '--no-contain': FlagParameter('--', 'no-contain'),
            '--no-overlap': FlagParameter('--', 'no-overlap'),
            '-t': FlagParameter('-', 't'),
            '--un': ValuedParameter('--', 'un', Delimiter=' ', IsPath=True),
            '--al': ValuedParameter('--', 'al', Delimiter=' ', IsPath=True),
            '--un-conc': ValuedParameter('--', 'un-conc', Delimiter=' ',
                    IsPath=True),
            '--al-conc': ValuedParameter('--', 'al-conc', Delimiter=' ',
                    IsPath=True),
            '--un-gz': FlagParameter('--', 'un-gz'),
            '--un-bz2': FlagParameter('--', 'un-bz2'),
            '--quiet': FlagParameter('--', 'quiet'),
            '--met-file': ValuedParameter('--', 'met-file', Delimiter=' ',
                    IsPath=True),
            '--met-stderr': FlagParameter('--', 'met-stderr'),
            '--met': ValuedParameter('--', 'met', Delimiter=' '),
            '--no-head': FlagParameter('--', 'no-head'),
            '--no-sq': FlagParameter('--', 'no-sq'),
            '--rg-id': ValuedParameter('--', 'rg-id', Delimiter=' '),
            '--rg': ValuedParameter('--', 'rg', Delimiter=' '),
            '--omit-sec-seq': FlagParameter('--', 'omit-sec-seq'),
            '-o': ValuedParameter('-', 'o', Delimiter=' '),
            '-p': ValuedParameter('-', 'p', Delimiter=' '),
            '--reorder': FlagParameter('--', 'reorder'),
            '--mn': FlagParameter('--', 'mn'),
            '--qc-filter': FlagParameter('--', 'qc-filter'),
            '--seed': ValuedParameter('--', 'seed', Delimiter=' '),
            '--non-deterministic': FlagParameter('--', 'non-deterministic'),
            '--version': FlagParameter('--', 'version'),
            '-h': FlagParameter('-', 'h')
            }

    _command = 'bowtie2'

