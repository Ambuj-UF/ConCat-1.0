################################################################################################################
# Basic alignment Handling                                                                                     #
#                                                                                                              #
# Copyright (C) {2014}  {Ambuj Kumar, Kimball-Brain lab group, Biology Department, University of Florida}      #
#                                                                                                              #
# This program is free software: you can redistribute it and/or modify                                         #
# it under the terms of the GNU General Public License as published by                                         #
# the Free Software Foundation, either version 3 of the License, or                                            #
# (at your option) any later version.                                                                          #
#                                                                                                              #
# This program is distributed in the hope that it will be useful,                                              #
# but WITHOUT ANY WARRANTY; without even the implied warranty of                                               #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                                #
# GNU General Public License for more details.                                                                 #
#                                                                                                              #
# This program comes with ABSOLUTELY NO WARRANTY;                                                              #
# This is free software, and you are welcome to redistribute it                                                #
# under certain conditions;                                                                                    #
#                                                                                                              #
################################################################################################################

import sys
import os

fullpath = os.getcwd() + "/src/Utils"
sys.path.append(fullpath)



try:
    import src.Utils.Bio
except ImportError, e:
    sys.exit('BioPython not found on your system. Program Terminated')



try:
    import argparse
except ImportError, e:
    sys.exit('Requires python 2.7 and above to run ConCat-1.0. Program Terminated')


import textwrap
from src.nex import richNexusCall
from src.functions import Convert, ConvertAll



parser = argparse.ArgumentParser(prog='ConCat',
                                 version= 'ConCat-1.0',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=textwrap.dedent('''\
                                                            ----------------------------------------------------------------------------------------------------------
                                                            \t\t\t\t Welcome to ConCat-1.0 
                                                            \t\t\t A sequence alignment Concatenation Program
                                                            \t\t\tWritten by Ambuj Kumar, University of Florida
                                                            \t      ConCat sequence alignment concatenation program was designed at Kimball-Braun
                                                            \t\tlab group to handle large amount of sequence dataset for phylogeny study.
                                                            
                                                            ----------------------------------------------------------------------------------------------------------
                                                            
                                                            
                                                            '''))


group = parser.add_mutually_exclusive_group()

parser.add_argument('-i', type=str, default=None,
                    help='Enter the data folder name')

parser.add_argument('-log', type=str, default=None,
                    help='Enter the log file name')

parser.add_argument('-o', type=str, default=None,
                    help='Enter the output file name')

parser.add_argument('-ftype', type=str, default='nexus',
                    choices=['fasta', 'nexus', 'phylip', 'phylip-interleived', 'phylip-relaxed'],
                    help='Enter the input file format for Concatenation. Default is nexus.')

parser.add_argument('-otype', type=str, default='nexus',
                    choices=['fasta', 'nexus', 'phylip', 'phylip-interleived', 'phylip-relaxed'],
                    help='Enter the output file format for Concatenation. Default is nexus.')

parser.add_argument('-convert', action='store_true', default=False,
                    help='Converts fasta and phylip alignment files to nexus alignment')

parser.add_argument('-CA', action='store_true', default=False,
                    help='Use this function if you want to perform file conversion and analysis simultaneously')

parser.add_argument('-spell', action='store_true', default=False,
                    help='Include if you want to check spelling mistakes in alignment files')

parser.add_argument('-block', action='store_true', default=False,
                    help='Include if you have ConCat block defined in alignment files')

parser.add_argument('-pipe', action='store_true', default=False,
                    help='Include if you have IDs stored in taxon names in the alignment file')

parser.add_argument('-RNA', action='store_true', default=False,
                    help='Include if you want to run RNAfold structure prediction')

group.add_argument('-inc', nargs='?', type=argparse.FileType('r'), default=sys.stdin,
                    help='Enter the text file with taxon names to be included in the alignment file')

group.add_argument('-exc', nargs='?', type=argparse.FileType('r'), default=sys.stdin,
                    help='Enter the tet file with taxon names to be removed from the alignment file')

parser.add_argument('-shannon', action='store_true', default=False,
                    help='Include if you want to run Shannons Entropy Calculation for the alignments')

parser.add_argument('-rcv', action='store_true', default=False,
                    help='Include if you want to run RCV Calculation for the alignments')

parser.add_argument('-GC', action='store_true', default=False,
                    help='Include if you want to calculate GC content')

parser.add_argument('-addT', type=str, default=None,
                    choices=['Family', 'Order', 'Class', 'Phylum', 'Kingdom', 'Family-Order', 'Family-Order-Class', 'Family-Order-Class-Phylum', 'Family-Order-Class-Phylum-Kingdom'],
                    help='Enter the txaonomy classes to add in the final alignment taxa name. Sperate Multiple values using comma (-addT Class,Family,Order)')

parser.add_argument('-remT', type=int, default=None,
                    help='Enter the numer of txaonomy classes to remove from the end of taxon namees')

parser.add_argument('-RY', nargs='?', type=argparse.FileType('r'), default=sys.stdin,
                    help='Enter the text file with alignment file names and positions for RY coding')

parser.add_argument('-OV', type=float, default=None,
                    help='Enter the OV cutoff value for detecting fast evolving sites')

parser.add_argument('-rbin', type=str, default=None,
                    help='Enter the RCV range to bin data')

parser.add_argument('-ebin', type=str, default=None,
                    help='Enter the entropy range to bin data')

parser.add_argument('-gcbin', type=str, default=None,
                    help='Enter the GC range in percentage to bin data')

parser.add_argument('-pbin', action='store_true', default=False,
                    help='Use this function if you want to store RCV, GC and/or Entropy data in percentile bins')

parser.add_argument('-ugcbin', type=str, default=None,
                    help='Use this function if you want to create GC bin by supplied GC percentile partition.')


args = parser.parse_args()


if args.RNA == True and not args.block:
    parser.error('-block argument is required in "-RNA" mode.')

if args.rbin == True and not args.rcv:
    parser.error('-rcv argument is required in "-rbin" mode.')

if args.ebin == True and not args.shannon:
    parser.error('-shannon argument is required in "-ebin" mode.')

if args.gcbin == True and not args.GC:
    parser.error('-GC argument is required in "-gcbin" mode.')

if args.ugcbin != None and args.pbin == False:
    parser.error('-pbin argument is required with -ugcbin.')

if args.pbin == True:
    if args.GC == True or args.shannon == True or args.rcv == True:
        pass
    else:
        parser.error('-rcv,-shannon or -GC argument is required in "-pbin" mode.')


if args.log != None:
    logfile = open(args.log, "w")
    logfile.write("ConCat-build log file\n\n")
    logfile.write("Command used -\n")
    logfile.write("python " + " ".join(sys.argv) + "\n\n")
    logfile.write("Operation log\n\n")
    sys.stdout = open(args.log, 'a')
    logfile.close()


def main():
    if args.convert == True:
        ConvertAll(args.ftype)
    
    else:
        if args.ftype == 'nexus' and args.otype == 'nexus':
            richNexusCall(args.RNA,
                          args.inc,
                          args.exc,
                          args.shannon,
                          args.rcv,
                          args.addT,
                          args.remT,
                          args.pipe,
                          args.RY,
                          args.spell,
                          args.block,
                          args.OV,
                          args.rbin,
                          args.ebin,
                          args.GC,
                          args.gcbin,
                          args.pbin,
                          args.ugcbin,
                          args.i,
                          args.o
                          )
    
        else:
            if args.CA == True:
                ConvertAll(args.ftype)
            
            richNexusCall(args.RNA,
                          args.inc,
                          args.exc,
                          args.shannon,
                          args.rcv,
                          args.addT,
                          args.remT,
                          args.pipe,
                          args.RY,
                          args.spell,
                          args.block,
                          args.OV,
                          args.rbin,
                          args.ebin,
                          args.GC,
                          args.gcbin,
                          args.pbin,
                          args.ugcbin,
                          args.i,
                          args.o
                          )

            if args.o != None:
                Convert('nexus', args.otype, args.o)
            else:
                Convert('nexus', args.otype, 'Combined.nex')





if __name__ == "__main__":
    main()