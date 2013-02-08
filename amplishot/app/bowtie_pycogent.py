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
__version__ = "0.1.0"
__maintainer__ = "Connor Skennerton"
__email__ = "c.skennerton@gmail.com"
__status__ = "Development"

###############################################################################

from cogent.app.util import CommandLineApplication
from cogent.app.parameters import FlagParameter, ValuedParameter
from amplishot.app.util import ExtendedCommandLineApplication,\
RepeatedParameter


class Bowtie2Build(CommandLineApplication):
    """ Simple wrapper for bowtie index
    """
    _parameters = {
            '-f': FlagParameter('-', 'f'),
            '-c': FlagParameter('-', 'c'),
            '-a': FlagParameter('-', 'a'),
            '-p': FlagParameter('-', 'p'),
            '--bmax': ValuedParameter('--', 'bmax', Delimeter=' '),
            '--bmaxdivn': ValuedParameter('--', 'bmaxdivn', Delimeter=' '),
            '--dcv': ValuedParameter('--', 'dcv', Delimeter=' '),
            '--nodc': ValuedParameter('--', 'nodc', Delimeter=' '),
            '-r': FlagParameter('-', 'r'),
            '-3': FlagParameter('-', '3'),
            '-o': ValuedParameter('-', 'o', Delimeter=' '),
            '-t': ValuedParameter('-', 't', Delimeter=' '),
            '--seed': ValuedParameter('--', 'seed', Delimeter=' '),
            '--cutoff': ValuedParameter('--', 'cutoff', Delimeter=' '),
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
            '-x': ValuedParameter('-', 'x', Delimeter=' ', IsPath=True),
            '-1': RepeatedParameter('-', '1', Delimeter=' ', IsPath=True),
            '-2': RepeatedParameter('-', '2', Delimeter=' ', IsPath=True),
            '-U': RepeatedParameter('-', 'U', Delimeter=' ', IsPath=True),
            '-S': ValuedParameter('-', 'S', Delimeter=' ', IsPath=True),
            '-q': FlagParameter('-', 'q'),
            '-f': FlagParameter('-', 'f'),
            '-r': FlagParameter('-', 'r'),
            '-c': FlagParameter('-', 'c'),
            '-s': ValuedParameter('-', 's', Delimeter=' '),
            '-u': ValuedParameter('-', 'u', Delimeter=' '),
            '-5': ValuedParameter('-', '5', Delimeter=' '),
            '-3': ValuedParameter('-', '3', Delimeter=' '),
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
            '-N': ValuedParameter('-', 'N', Delimeter=' '),
            '-L': ValuedParameter('-', 'L', Delimeter=' '),
            '-i': ValuedParameter('-', 'i', Delimeter=' '),
            '--n-ceil': ValuedParameter('--', 'n-ceil', Delimeter=' '),
            '--dpad': ValuedParameter('--', 'dpad', Delimeter=' '),
            '--gbar': ValuedParameter('--', 'gbar', Delimeter=' '),
            '--ignore-quals': ValuedParameter('--', 'ignore-quals', Delimeter=' '),
            '--nofw': ValuedParameter('--', 'nofw', Delimeter=' '),
            '--norc': ValuedParameter('--', 'norc', Delimeter=' '),
            '--end-to-end': FlagParameter('--', 'end-to-end'),
            '--local': FlagParameter('--', 'local'),
            '--ma': ValuedParameter('--', 'ma', Delimeter=' '),
            '--mp': ValuedParameter('--', 'mp', Delimeter=' '),
            '--np': ValuedParameter('--', 'np', Delimeter=' '),
            '--rdg': ValuedParameter('--', 'rdg', Delimeter=' '),
            '--rfg': ValuedParameter('--', 'rfg', Delimeter=' '),
            '--score-min': ValuedParameter('--', 'score-min', Delimeter=' '),
            '-k': ValuedParameter('-', 'k', Delimeter=' '),
            '-a': FlagParameter('-', 'a'),
            '-D': ValuedParameter('-', 'D', Delimeter=' '),
            '-R': ValuedParameter('-', 'R', Delimeter=' '),
            '-I': ValuedParameter('-', 'I', Delimeter=' '),
            '-X': ValuedParameter('-', 'X', Delimeter=' '),
            '--fr': FlagParameter('--', 'fr'),
            '--rf': FlagParameter('--', 'rf'),
            '--ff': FlagParameter('--', 'ff'),
            '--no-mixed': FlagParameter('--', 'no-mixed'),
            '--no-discordant': FlagParameter('--', 'no-discordant'),
            '--no-dovetail': FlagParameter('--', 'no-dovetail'),
            '--no-contain': FlagParameter('--', 'no-contain'),
            '--no-overlap': FlagParameter('--', 'no-overlap'),
            '-t': Flagparameter('-', 't'),
            '--un': ValuedParameter('--', 'un', Delimeter=' ', IsPath=True),
            '--al': ValuedParameter('--', 'al', Delimeter=' ', IsPath=True),
            '--un-conc': ValuedParameter('--', 'un-conc', Delimeter=' ',
                    IsPath=True),
            '--al-conc': ValuedParameter('--', 'al-conc', Delimeter=' ',
                    IsPath=True),
            '--un-gz': FlagParameter('--', 'un-gz'),
            '--un-bz2': FlagParameter('--', 'un-bz2'),
            '--quiet': FlagParameter('--', 'quiet'),
            '--met-file': ValuedParameter('--', 'met-file', Delimeter=' ',
                    IsPath=True),
            '--met-stderr': FlagParameter('--', 'met-stderr'),
            '--met': ValuedParameter('--', 'met', Delimeter=' '),
            '--no-head': FlagParameter('--', 'no-head'),
            '--no-sq': FlagParameter('--', 'no-sq'),
            '--rg-id': ValuedParameter('--', 'rg-id', Delimeter=' '),
            '--rg': ValuedParameter('--', 'rg', Delimeter=' '),
            '--omit-sec-seq': FlagParameter('--', 'omit-sec-seq'),
            '-o': ValuedParameter('-', 'o', Delimeter=' '),
            '-p': ValuedParameter('-', 'p', Delimeter=' '),
            '--reorder': FlagParameter('--', 'reorder'),
            '--mn': FlagParameter('--', 'mn'),
            '--qc-filter': FlagParameter('--', 'qc-filter'),
            '--seed': ValuedParameter('--', 'seed', Delimeter=' '),
            '--non-deterministic': FlagParameter('--', 'non-deterministic'),
            '--version': FlagParameter('--', 'version'),
            '-h': FlagParameter('-', 'h')
            }

    _command = 'bowtie2'
