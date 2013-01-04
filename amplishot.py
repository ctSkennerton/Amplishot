#!/usr/bin/env python
###############################################################################
#
# amplishot.py - Pipeline for generating full-length 16S sequences and performing
#                community abundance measurements on the output#
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
__copyright__ = "Copyright 2012-2013"
__credits__ = ["Connor Skennerton"]
__license__ = "GPL3"
__version__ = "0.0.1"
__maintainer__ = "Connor Skennerton"
__email__ = "c.skennerton@gmail.com"
__status__ = "Development"

###############################################################################

import argparse
import sys
from pandaseq import Pandaseq
#import os
#import errno

#import numpy as np
#np.seterr(all='raise')     

#import matplotlib as mpl
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import axes3d, Axes3D
#from pylab import plot,subplot,axis,stem,show,figure


###############################################################################
###############################################################################
###############################################################################
###############################################################################

  # classes here

###############################################################################
###############################################################################
###############################################################################
###############################################################################
def panda(args):
    panda = Pandaseq()
    panda.Parameters['-o'].on(args.olap)
    panda.Parameters['-l'].on(args.minlen)
    panda.Parameters['-f'].on(args.read_1)  
    panda.Parameters['-r'].on(args.read_2)
    return panda()
    
def doWork( args ):
    """ Main wrapper"""

    # generate fragment consensus sequences - pandaseq
    panda_out = panda(args)

    # dataset reduction - velvet or cd-hit

    # generate overlaps - phrap

    #read alignment - bwa or bowtie

    return 0

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    #parser.add_argument('positional_arg', help="Required")
    #parser.add_argument('positional_arg2', type=int, help="Integer argument")
    #parser.add_argument('positional_arg3', nargs='+', help="Multiple values")
    #parser.add_argument('-X', '--optional_X', action="store_true", default=False, help="flag")
    
    # parse the arguments
    args = parser.parse_args()        

    # do what we came here to do
    return doWork(args)

###############################################################################
###############################################################################
###############################################################################
###############################################################################
